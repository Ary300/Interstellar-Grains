from __future__ import annotations

import hashlib
import json
import pickle
import random
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np

from scientific_data import (
    K_B,
    K_B_ERG,
    M_H,
    surface_chemistry_data,
    uv_photon_flux,
)
from physical_rates import (
    adsorption_rate,
    h_desorption_rate,
    h_diffusion_rate,
    h2_formation_lh_rate,
    uv_h2_formation_rate,
    uv_h2_photofragmentation_rate,
    uv_photodesorption_rate,
)


def _thermal_rate(prefactor_s: float, barrier_eV: float, temperature_k: float) -> float:
    if temperature_k <= 0:
        return 0.0
    if barrier_eV <= 0:
        return float(prefactor_s)
    return float(prefactor_s) * float(np.exp(-float(barrier_eV) / (float(K_B) * float(temperature_k))))


class GrainCache:
    """
    Disk cache for reusing pre-generated grain topology + energetics across runs.

    Notes:
    - Cache contains only *structure* and energy maps. Occupancy (H/H2) is always fresh per run.
    - By default, the cache key intentionally excludes rng_seed so that repeated ensemble runs reuse
      a single grain realization (faster + reduces variance from structure re-rolls).
    """

    CACHE_VERSION = 2  # bump to invalidate old caches when layout/meaning changes
    _MEM_BYTES: Dict[str, bytes] = {}
    _MEM_DATA: Dict[str, Dict[str, Any]] = {}

    def __init__(self, cache_dir: str = "grain_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _key(self, params: Dict[str, Any], include_rng_seed: bool) -> str:
        key_fields = {
            "v": int(self.CACHE_VERSION),
            "grain_radius_um": float(params.get("grain_radius_um", 0.1)),
            "site_area_angstroms_sq": float(params.get("site_area_angstroms_sq", 9.0)),
            "use_3d_lattice": bool(params.get("use_3d_lattice", True)),
            "porosity_fraction": float(params.get("porosity_fraction", 0.2)),
            "chemisorption_fraction": float(params.get("chemisorption_fraction", 0.1)),
            "surface_defect_fraction": float(params.get("surface_defect_fraction", 0.15)),
            "E_phys_mean_meV": float(params.get("E_phys_mean_meV", 45.0)),
            "heterogeneity_E_bind_sigma_meV": float(params.get("heterogeneity_E_bind_sigma_meV", 5.0)),
            "E_chem_mean_eV": float(params.get("E_chem_mean_eV", 1.75)),
            "heterogeneity_E_chem_sigma_eV": float(params.get("heterogeneity_E_chem_sigma_eV", 0.25)),
        }
        if include_rng_seed:
            # When enabled, you get separate cached grains per RNG seed.
            key_fields["rng_seed"] = int(params.get("rng_seed", 0) or 0)
        payload = json.dumps(key_fields, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.md5(payload).hexdigest()

    def _path(self, key: str) -> Path:
        return self.cache_dir / f"grain_{key}.pkl"

    def load(self, key: str) -> Optional[Dict[str, Any]]:
        p = self._path(key)
        if not p.exists():
            return None
        cache_key = str(p)
        cached_obj = self._MEM_DATA.get(cache_key, None)
        if cached_obj is not None:
            return cached_obj
        b = self._MEM_BYTES.get(cache_key, None)
        if b is None:
            b = p.read_bytes()
            self._MEM_BYTES[cache_key] = b
        obj = pickle.loads(b)
        if isinstance(obj, dict):
            self._MEM_DATA[cache_key] = obj
        return obj

    def save(self, key: str, data: Dict[str, Any]) -> None:
        p = self._path(key)
        b = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
        p.write_bytes(b)
        self._MEM_BYTES[str(p)] = b
        self._MEM_DATA[str(p)] = data


class KineticMonteCarlo:
    def __init__(self, simulation_parameters: Dict[str, Any]):
        self.simulation_parameters = dict(simulation_parameters or {})

        seed = self.simulation_parameters.get("rng_seed", None)
        try:
            seed = int(seed) if seed is not None and str(seed).strip() != "" else None
        except (TypeError, ValueError):
            seed = None
        self.rng_seed: Optional[int] = seed
        if self.rng_seed is not None:
            random.seed(int(self.rng_seed))

        # Core time
        self.time = 0.0
        self.last_delta_t: Optional[float] = None

        # Event counters
        self.total_arrivals = 0  # total 'arrival' events executed (atoms + H2) in arrival mode
        self.total_impinging_h_atoms = 0
        self.total_impinging_h2_molecules = 0

        self.total_adsorbed_h_atoms = 0
        self.total_desorbed_h_atoms = 0

        self.h_atoms_on_surface = 0
        self.adjacent_h_pairs_count = 0

        self.h2_molecules_formed = 0
        self.h2_molecules_formed_LH = 0
        self.h2_molecules_formed_ER = 0
        self.h2_molecules_formed_UV = 0

        # Observable bookkeeping for Grieco-style ε: prompt desorption of *formed* H2 only.
        self.h2_molecules_desorbed = 0
        self.h2_molecules_desorbed_LH = 0
        self.h2_molecules_desorbed_ER = 0
        self.h2_molecules_desorbed_UV = 0

        # Beam-origin H2 bookkeeping (baseline channel)
        self.h2_molecules_desorbed_beam = 0
        self.h2_molecules_released_beam = 0
        self.h2_molecules_released_formed = 0

        # H2 blocking state
        self.h2_molecules_on_surface = 0
        self.h2_sites_formed: set[Tuple[int, int]] = set()
        self.h2_sites_beam: set[Tuple[int, int]] = set()

        # UV state
        self.uv_mode = str(self.simulation_parameters.get("uv_mode", "pulse") or "pulse").strip().lower()
        self.uv_pulse_enabled = bool(self.simulation_parameters.get("uv_pulse_enabled", True))
        self.uv_defect_creation_rate = float(self.simulation_parameters.get("uv_defect_creation_rate", 0.5))
        self.uv_pulse_duration = float(self.simulation_parameters.get("uv_pulse_duration", 1e-6))
        self.last_uv_pulse_time = 0.0
        self.uv_pulse_active = False

        # Energetics knobs
        self.E_bind_mean_meV = float(self.simulation_parameters.get("E_phys_mean_meV", 45.0))
        self.E_bind_sigma_meV = float(self.simulation_parameters.get("heterogeneity_E_bind_sigma_meV", 5.0))

        # Debug printing of Gillespie time steps
        self.debug_print_timestep = bool(self.simulation_parameters.get("debug_print_timestep", False))
        self.debug_timestep_every = int(self.simulation_parameters.get("debug_timestep_every", 1))

        # Initialize grain structure (+ optional cache)
        self.use_3d_lattice = bool(self.simulation_parameters.get("use_3d_lattice", True))
        self.porosity_fraction = float(self.simulation_parameters.get("porosity_fraction", 0.2))

        # Mechanism toggles / modeling choices
        # If enabled, LH (and UV-assisted) formation only uses non-chemisorption site pairs.
        # This preserves a chemisorbed reservoir for the high-T ER plateau during DED ramps.
        self.lh_exclude_chemisorption = bool(self.simulation_parameters.get("lh_exclude_chemisorption", False))
        self.initialize_3d_lattice()

        # Internal caches for rate-weighted micro-event selection (refreshed in calculate_rates).
        self._last_desorption_sites: list[Tuple[int, int, int]] = []
        self._last_desorption_weights: list[float] = []
        self._last_diffusion_sites: list[Tuple[int, int, int]] = []
        self._last_diffusion_weights: list[float] = []
        self._last_diffusion_empty_neighbors: list[list[Tuple[int, int, int]]] = []

    def _arrival_mode_enabled(self) -> bool:
        # Explicit arrival-rate knobs
        a = self.simulation_parameters.get("arrival_rate_s", None)
        b = self.simulation_parameters.get("arrival_rate_per_site_s", None)
        try:
            if a is not None and float(a) > 0:
                return True
        except (TypeError, ValueError):
            pass
        try:
            if b is not None and float(b) > 0:
                return True
        except (TypeError, ValueError):
            pass

        # Flux-based beam knob (converted to a total arrival rate using the current surface area).
        f = self.simulation_parameters.get("beam_flux_total_cm2_s", None)
        try:
            if f is not None and float(f) > 0:
                return True
        except (TypeError, ValueError):
            pass
        return False

    def _bool_param(self, name: str, default: bool) -> bool:
        v = self.simulation_parameters.get(name, default)
        if v is None:
            return default
        return bool(v)

    def _update_surface_temperature_from_ramp(self) -> None:
        ramp = self.simulation_parameters.get("temp_ramp", None)
        if not isinstance(ramp, dict) or not bool(ramp.get("enabled", False)):
            return
        T0 = float(ramp.get("T_start_K", self.simulation_parameters.get("surface_temperature_k", 10.0)))
        T1 = float(ramp.get("T_end_K", T0))
        t0 = float(ramp.get("t0_s", 0.0))
        rate = None
        if "rate_K_per_s" in ramp:
            rate = float(ramp["rate_K_per_s"])
        elif "rate_K_per_min" in ramp:
            rate = float(ramp["rate_K_per_min"]) / 60.0
        if rate is None:
            return
        dt = max(0.0, float(self.time) - t0)
        T = float(T0) + float(rate) * dt
        if T1 >= T0:
            T = min(float(T1), max(float(T0), T))
        else:
            T = max(float(T1), min(float(T0), T))
        self.simulation_parameters["surface_temperature_k"] = float(T)

    def initialize_3d_lattice(self) -> None:
        enable_cache = bool(self.simulation_parameters.get("enable_grain_cache", False))
        cache_dir = str(self.simulation_parameters.get("grain_cache_dir", "grain_cache"))
        include_seed = bool(self.simulation_parameters.get("grain_cache_include_rng_seed", False))

        cache = GrainCache(cache_dir) if enable_cache else None
        cache_key = cache._key(self.simulation_parameters, include_rng_seed=include_seed) if cache else None
        cached = cache.load(cache_key) if cache and cache_key else None

        if cached is None:
            grain_radius_um = float(self.simulation_parameters.get("grain_radius_um", 0.1))
            grain_radius_cm = grain_radius_um * 1e-4
            site_area_angstroms_sq = float(self.simulation_parameters.get("site_area_angstroms_sq", 9.0))
            site_area_cm2 = site_area_angstroms_sq * 1e-16
            grain_surface_area_cm2 = 4.0 * float(np.pi) * (grain_radius_cm**2)
            calculated_sites = max(1, int(grain_surface_area_cm2 / max(site_area_cm2, 1e-30)))
            surface_dimension = max(2, int(np.sqrt(calculated_sites)))
            depth_layers = 1 if not self.use_3d_lattice else max(3, surface_dimension // 10)

            # Deterministic structure RNG when caching is enabled; otherwise, tie to rng_seed (if provided).
            if enable_cache and cache_key is not None:
                seed_int = int(cache_key[:8], 16)
            else:
                seed_int = int(self.rng_seed) if self.rng_seed is not None else None
            rng = np.random.default_rng(seed_int)

            lattice = np.full((depth_layers, surface_dimension, surface_dimension), "C", dtype=object)
            porosity_mask = rng.random(lattice.shape) < float(self.porosity_fraction)
            lattice[porosity_mask] = None

            # Ensure at least a few accessible top-layer sites exist.
            if not np.any(lattice[0, :, :] != None):
                num_accessible = max(1, int(surface_dimension * surface_dimension * 0.1))
                idxs = rng.choice(surface_dimension * surface_dimension, num_accessible, replace=False)
                for idx in idxs:
                    r, c = int(idx) // surface_dimension, int(idx) % surface_dimension
                    lattice[0, r, c] = "C"

            site_types = np.where(lattice != None, 1, 0).astype(int)  # 0=void, 1=phys

            # Surface-only chemisorption/defect assignment.
            chem_frac = float(self.simulation_parameters.get("chemisorption_fraction", 0.1))
            def_frac = float(self.simulation_parameters.get("surface_defect_fraction", 0.15))
            surf_acc = np.argwhere(lattice[0] != None)
            n_surf = int(len(surf_acc))
            if n_surf > 0:
                n_chem = int(round(n_surf * max(0.0, min(1.0, chem_frac))))
                if n_chem > 0:
                    chem_idx = rng.choice(n_surf, n_chem, replace=False)
                    for i in chem_idx:
                        r, c = int(surf_acc[i][0]), int(surf_acc[i][1])
                        site_types[0, r, c] = 2

                # Defects assigned only on remaining physisorption sites.
                remaining = np.argwhere((lattice[0] != None) & (site_types[0] == 1))
                n_rem = int(len(remaining))
                n_def = int(round(n_surf * max(0.0, min(1.0, def_frac))))
                if n_def > 0 and n_rem > 0:
                    n_def = min(n_def, n_rem)
                    def_idx = rng.choice(n_rem, n_def, replace=False)
                    for i in def_idx:
                        r, c = int(remaining[i][0]), int(remaining[i][1])
                        site_types[0, r, c] = 3

            # Energetics: physisorption vs chemisorption vs defect.
            E_bind = np.zeros(lattice.shape, dtype=float)
            E_diff = np.zeros(lattice.shape, dtype=float)

            mean_phys = float(self.simulation_parameters.get("E_phys_mean_meV", self.E_bind_mean_meV)) / 1000.0
            sig_phys = float(self.simulation_parameters.get("heterogeneity_E_bind_sigma_meV", self.E_bind_sigma_meV)) / 1000.0
            mean_chem = float(self.simulation_parameters.get("E_chem_mean_eV", 1.75))
            sig_chem = float(self.simulation_parameters.get("heterogeneity_E_chem_sigma_eV", 0.25))
            mean_def = float(self.simulation_parameters.get("E_defect_mean_eV", 0.35))
            sig_def = float(self.simulation_parameters.get("heterogeneity_E_defect_sigma_eV", 0.05))

            phys_mask = (lattice != None) & (site_types == 1)
            chem_mask = (lattice != None) & (site_types == 2)
            def_mask = (lattice != None) & (site_types == 3)

            if np.any(phys_mask):
                E_bind[phys_mask] = rng.normal(mean_phys, sig_phys, int(np.sum(phys_mask)))
            if np.any(chem_mask):
                E_bind[chem_mask] = rng.normal(mean_chem, sig_chem, int(np.sum(chem_mask)))
            if np.any(def_mask):
                E_bind[def_mask] = rng.normal(mean_def, sig_def, int(np.sum(def_mask)))

            E_bind = np.clip(E_bind, 0.0, None)

            # Diffusion barrier map (not currently used for rates, but cached for completeness/tests).
            E_diff[phys_mask] = 0.3 * E_bind[phys_mask]
            E_diff[chem_mask] = 0.3 * E_bind[chem_mask]
            E_diff[def_mask] = 0.3 * E_bind[def_mask]

            cached = {
                "surface_dimension": int(surface_dimension),
                "depth_layers": int(depth_layers),
                "lattice_base": lattice,
                "site_types": site_types.astype(int),
                "E_bind_eV_map": E_bind.astype(float),
                "E_diff_eV_map": E_diff.astype(float),
            }
            if cache and cache_key:
                cache.save(cache_key, cached)

        self.surface_dimension = int(cached["surface_dimension"])
        self.depth_layers = int(cached["depth_layers"])
        self.lattice = np.array(cached["lattice_base"], copy=True)
        # These maps are treated as read-only during simulation. Avoid extra copies to speed up
        # large ensemble campaigns where we repeatedly instantiate KMC objects.
        self.site_types = np.array(cached["site_types"], copy=False)
        self.E_bind_eV_map = np.array(cached["E_bind_eV_map"], copy=False)
        self.E_diff_eV_map = np.array(cached["E_diff_eV_map"], copy=False)
        for arr in (self.site_types, self.E_bind_eV_map, self.E_diff_eV_map):
            try:
                arr.setflags(write=False)
            except Exception:
                pass

        self.total_accessible_surface_sites = int(np.sum(self.lattice[0, :, :] != None))

        # Fast occupancy structures
        self.occupied: set[Tuple[int, int, int]] = set()
        top_accessible = np.where(self.lattice[0, :, :] != None)
        self.empty_surface: set[Tuple[int, int]] = set(zip(top_accessible[0].tolist(), top_accessible[1].tolist()))
        self.occupied_chemisorption_surface: set[Tuple[int, int]] = set()

        # Seed initial H (fresh; not part of the cache)
        initial_h_coverage = float(self.simulation_parameters.get("initial_h_coverage", 0.0) or 0.0)
        initial_h_count = self.simulation_parameters.get("initial_h_count", None)
        try:
            initial_h_count_i = int(initial_h_count) if initial_h_count is not None else None
        except (TypeError, ValueError):
            initial_h_count_i = None
        self._initialize_h_atoms(initial_coverage=initial_h_coverage, initial_count=initial_h_count_i)

    def _initialize_h_atoms(self, initial_coverage: float = 0.0, initial_count: Optional[int] = None) -> None:
        if (initial_count is None or initial_count <= 0) and (initial_coverage is None or float(initial_coverage) <= 0):
            return

        # Candidate sites: empty accessible surface sites.
        candidates = list(self.empty_surface)
        if not candidates:
            return

        chem_only = self.simulation_parameters.get("initial_h_chemisorption_only", None)
        if chem_only:
            candidates = [(r, c) for (r, c) in candidates if int(self.site_types[0, r, c]) == 2]
            if not candidates:
                return

        if isinstance(initial_count, int) and initial_count > 0:
            n = min(int(initial_count), len(candidates))
        else:
            n = int(round(len(candidates) * float(initial_coverage)))
            n = max(0, min(n, len(candidates)))

        if n <= 0:
            return

        chosen = random.sample(candidates, n)
        for r, c in chosen:
            self._adsorb_h_at_surface(int(r), int(c))
        self._update_adjacent_h_pairs_count()

    def get_neighbors_3d(self, d: int, r: int, c: int):
        neighbors = []
        moves_2d = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dr, dc in moves_2d:
            nr, nc = (r + dr) % self.surface_dimension, (c + dc) % self.surface_dimension
            if self.lattice[d, nr, nc] is not None:
                neighbors.append((d, nr, nc))

        if d < self.depth_layers - 1 and self.lattice[d + 1, r, c] is not None:
            neighbors.append((d + 1, r, c))

        if d > 0 and self.lattice[d - 1, r, c] is not None:
            neighbors.append((d - 1, r, c))

        return neighbors

    def _update_adjacent_h_pairs_count(self) -> None:
        self.adjacent_h_pairs_count = 0
        for d, r, c in self.get_occupied_sites():
            for nd, nr, nc in self.get_neighbors_3d(d, r, c):
                if (
                    self.lattice[nd, nr, nc] == "H"
                    and (d, r, c) < (nd, nr, nc)
                    and self._lh_pair_allowed(int(d), int(r), int(c), int(nd), int(nr), int(nc))
                ):
                    self.adjacent_h_pairs_count += 1

    def update_adjacent_h_pairs_count(self, d: int, r: int, c: int, add_atom: bool) -> None:
        change = 1 if add_atom else -1
        for nd, nr, nc in self.get_neighbors_3d(d, r, c):
            if self.lattice[nd, nr, nc] == "H" and self._lh_pair_allowed(int(d), int(r), int(c), int(nd), int(nr), int(nc)):
                self.adjacent_h_pairs_count += change

    def get_accessible_surface_sites(self):
        if not self.empty_surface:
            return (np.array([], dtype=int), np.array([], dtype=int))
        rows, cols = zip(*self.empty_surface)
        return (np.fromiter(rows, dtype=int), np.fromiter(cols, dtype=int))

    def get_num_accessible_surface_sites(self) -> int:
        return int(len(self.empty_surface))

    def get_occupied_sites(self):
        return list(self.occupied) if self.occupied else []

    def _theta_h2(self) -> float:
        denom = float(self.total_accessible_surface_sites or 0)
        if denom <= 0:
            return 0.0
        return float(self.h2_molecules_on_surface) / denom

    def _h2_stick_probability(self, surface_temp_k: float) -> float:
        if not bool(self.simulation_parameters.get("enable_h2_blocking", False)):
            return 0.0
        t_transition = float(self.simulation_parameters.get("h2_stick_transition_K", 20.0))
        p_low = float(self.simulation_parameters.get("h2_stick_prob_lowT", 0.9))
        return float(p_low) if float(surface_temp_k) < float(t_transition) else 0.0

    def _diffusion_rate(self, site_type: int, surface_temp_k: float) -> float:
        """
        Diffusion is the fastest process at low T and can dominate Gillespie step counts.

        For the purposes of this project (and to keep DED ramps computationally tractable),
        we use a simple Arrhenius form with configurable barriers, rather than inheriting the
        very fast "defect" channel from `physical_rates.h_diffusion_rate`.
        """
        if surface_temp_k <= 0:
            return 0.0

        # Chemisorbed H is typically far less mobile than physisorbed H. Keep it configurable
        # because some unit tests / exploratory runs expect generic diffusion.
        if site_type == 2 and not bool(self.simulation_parameters.get("enable_chemisorption_diffusion", True)):
            return 0.0

        pref = float(self.simulation_parameters.get("diffusion_prefactor_s", 1e12))
        if site_type == 3:
            barrier = float(self.simulation_parameters.get("diffusion_barrier_defect_eV", 0.025))
        elif site_type == 2:
            barrier = float(self.simulation_parameters.get("diffusion_barrier_chem_eV", 0.03))
        else:
            barrier = float(self.simulation_parameters.get("diffusion_barrier_phys_eV", 0.025))
        rate = _thermal_rate(pref, barrier, float(surface_temp_k))
        cap = self.simulation_parameters.get("diffusion_rate_cap_s", None)
        if cap is not None:
            try:
                rate = min(float(rate), float(cap))
            except (TypeError, ValueError):
                pass
        return float(rate)

    def _lh_pair_allowed(self, d1: int, r1: int, c1: int, d2: int, r2: int, c2: int) -> bool:
        """
        Whether an adjacent H–H pair should contribute to LH/UV formation.

        When `lh_exclude_chemisorption` is enabled, pairs that include any chemisorption
        site are excluded. This prevents the LH channel from draining the chemisorbed
        reservoir during DED ramps, preserving the intended ER-driven high-T plateau.
        """
        if not bool(self.lh_exclude_chemisorption):
            return True
        try:
            st1 = int(self.site_types[int(d1), int(r1), int(c1)])
            st2 = int(self.site_types[int(d2), int(r2), int(c2)])
        except Exception:
            return True
        return st1 != 2 and st2 != 2

    def _handle_new_h2(self, place_site: Optional[Tuple[int, int]] = None, origin: str = "formed", mechanism: str = "") -> None:
        """
        Handle a newly formed (or beam-origin) H2 molecule.

        - If it sticks (low-T blocking), place it on the surface lattice at place_site.
        - Else, count it as prompt-desorbed (for formed) or prompt baseline (for beam).
        """
        T = float(self.simulation_parameters.get("surface_temperature_k", 10.0))
        p_stick = self._h2_stick_probability(T)

        if origin == "beam":
            base = float(self.simulation_parameters.get("h2_beam_stick_probability", 1.0))
            p_stick = min(1.0, max(0.0, base * p_stick))

        if place_site is not None and p_stick > 0 and random.random() < float(p_stick):
            r, c = int(place_site[0]), int(place_site[1])
            if self.lattice[0, r, c] == "C":
                self.lattice[0, r, c] = "H2"
                self.empty_surface.discard((r, c))
                self.h2_molecules_on_surface += 1
                if origin == "beam":
                    self.h2_sites_beam.add((r, c))
                else:
                    self.h2_sites_formed.add((r, c))
                return

        # Prompt release
        if origin == "beam":
            self.h2_molecules_desorbed_beam += 1
            self.h2_molecules_released_beam += 1
        else:
            self.h2_molecules_desorbed += 1
            if mechanism == "LH":
                self.h2_molecules_desorbed_LH += 1
            elif mechanism == "ER":
                self.h2_molecules_desorbed_ER += 1
            elif mechanism == "UV":
                self.h2_molecules_desorbed_UV += 1

    def _adsorb_h_at_surface(self, r: int, c: int) -> None:
        self.lattice[0, r, c] = "H"
        self.h_atoms_on_surface += 1
        self.total_adsorbed_h_atoms += 1
        self.occupied.add((0, r, c))
        self.empty_surface.discard((r, c))
        self.update_adjacent_h_pairs_count(0, r, c, True)
        if int(self.site_types[0, r, c]) == 2:
            self.occupied_chemisorption_surface.add((r, c))

    def _desorb_h_at(self, d: int, r: int, c: int) -> None:
        self.lattice[d, r, c] = "C"
        self.h_atoms_on_surface -= 1
        self.total_desorbed_h_atoms += 1
        self.update_adjacent_h_pairs_count(d, r, c, False)
        self.occupied.discard((d, r, c))
        if d == 0:
            if self.lattice[0, r, c] is not None:
                self.empty_surface.add((r, c))
            self.occupied_chemisorption_surface.discard((r, c))

    def calculate_rates(self) -> Dict[str, float]:
        # Keep surface temperature updated if a ramp is active.
        self._update_surface_temperature_from_ramp()

        surface_temp_k = float(self.simulation_parameters.get("surface_temperature_k", 10.0))
        gas_temp_k = float(self.simulation_parameters.get("gas_temperature_k", 100.0))
        h_gas_density = float(self.simulation_parameters.get("h_gas_density_cm3", 0.0))
        uv_flux_factor = float(self.simulation_parameters.get("uv_flux_factor", 0.0))

        site_area_cm2 = float(self.simulation_parameters.get("site_area_angstroms_sq", 9.0)) * 1e-16
        sticking_probability = float(self.simulation_parameters.get("sticking_probability", 0.3))

        rates: Dict[str, float] = {}

        arrival_mode = self._arrival_mode_enabled()

        # Arrival/beam mode takes precedence over gas adsorption mode.
        if arrival_mode:
            a = self.simulation_parameters.get("arrival_rate_s", None)
            b = self.simulation_parameters.get("arrival_rate_per_site_s", None)
            f = self.simulation_parameters.get("beam_flux_total_cm2_s", None)
            arrival_rate = 0.0
            try:
                if a is not None and float(a) > 0:
                    arrival_rate = float(a)
                elif b is not None and float(b) > 0:
                    arrival_rate = float(b) * float(self.total_accessible_surface_sites)
                elif f is not None and float(f) > 0:
                    # Convert beam flux (per cm^2 per s) into a total arrival rate on this surface.
                    flux = float(f)
                    angle = self.simulation_parameters.get("beam_incidence_angle_deg", None)
                    if angle is not None:
                        try:
                            flux *= float(np.cos(np.deg2rad(float(angle))))
                        except (TypeError, ValueError):
                            pass
                    surface_area_cm2 = float(self.total_accessible_surface_sites) * float(site_area_cm2)
                    arrival_rate = float(flux) * float(surface_area_cm2)
            except (TypeError, ValueError):
                arrival_rate = 0.0
            if arrival_rate > 0:
                rates["arrival"] = float(arrival_rate)
        else:
            # Gas adsorption rate: depends on empty accessible area.
            num_accessible_sites = self.get_num_accessible_surface_sites()
            accessible_area_cm2 = float(num_accessible_sites) * float(site_area_cm2)
            rates["adsorption"] = float(adsorption_rate(h_gas_density, gas_temp_k, sticking_probability, accessible_area_cm2))

        # H2 desorption (blocking channel)
        if self.h2_molecules_on_surface > 0:
            E_h2 = float(self.simulation_parameters.get("E_h2_bind_eV", 0.03))
            pref = float(self.simulation_parameters.get("h2_desorption_prefactor_s", 1e12))
            k = _thermal_rate(pref, E_h2, surface_temp_k)
            rates["h2_desorption"] = float(k) * float(self.h2_molecules_on_surface)

        # ER (gas-mode only). Arrival-mode ER is handled at the moment of arrival.
        if (not arrival_mode) and self.h_atoms_on_surface > 0:
            v_th = float(np.sqrt(8.0 * float(K_B_ERG) * float(gas_temp_k) / (float(np.pi) * float(M_H))))
            gas_flux = 0.25 * float(h_gas_density) * v_th  # atoms cm^-2 s^-1
            sigma = float(self.simulation_parameters.get("er_cross_section_cm2", surface_chemistry_data.get("er_cross_section_cm2", 1e-15)))
            p = float(self.simulation_parameters.get("er_reaction_probability", 0.1))
            rates["h2_formation_ER"] = float(gas_flux) * float(sigma) * float(p) * float(self.h_atoms_on_surface)

        # H desorption and diffusion
        self._last_desorption_sites = []
        self._last_desorption_weights = []
        self._last_diffusion_sites = []
        self._last_diffusion_weights = []
        self._last_diffusion_empty_neighbors = []

        diffusion_mode = str(self.simulation_parameters.get("diffusion_mode", "explicit") or "explicit").strip().lower()
        # LH mode options:
        # - "pairs" (default): use explicit adjacent pairs (requires diffusion to create adjacency)
        # - "diffusion_limited": compute LH from diffusion in a fast-mixing approximation (no explicit adjacency needed)
        lh_mode = str(self.simulation_parameters.get("lh_formation_mode", "pairs") or "pairs").strip().lower()

        # In diffusion-limited LH mode, keep track of an eligible surface subset.
        lh_surface_h_sites: list[Tuple[int, int, int]] = []
        lh_surface_sites_total = None
        total_diff_surface_lh = 0.0

        if self.h_atoms_on_surface > 0 and self.occupied:
            total_des = 0.0
            total_diff = 0.0
            # Default: diffusion is ON in gas-kinetic mode, OFF in arrival-mode unless explicitly enabled.
            enable_diffusion = self._bool_param("enable_diffusion", False if arrival_mode else True)

            for d, r, c in self.occupied:
                d_i, r_i, c_i = int(d), int(r), int(c)
                st = int(self.site_types[d_i, r_i, c_i])
                bind_e = float(self.E_bind_eV_map[d_i, r_i, c_i])
                # Chemisorbed H is effectively non-desorbing over 10–250 K in this model; skipping
                # the exp() call here is a large performance win for high-T plateau runs where the
                # chemisorbed reservoir can contain hundreds of sites.
                des_k = 0.0 if st == 2 else float(h_desorption_rate(bind_e, surface_temp_k))
                if des_k > 0:
                    total_des += des_k
                    self._last_desorption_sites.append((d_i, r_i, c_i))
                    self._last_desorption_weights.append(des_k)

                if enable_diffusion:
                    # Only mobile atoms contribute to diffusion.
                    empties = [nb for nb in self.get_neighbors_3d(d_i, r_i, c_i) if self.lattice[nb] == "C"]
                    if empties:
                        diff_k = float(self._diffusion_rate(st, surface_temp_k))
                        if diff_k > 0:
                            total_diff += diff_k
                            if diffusion_mode == "explicit":
                                self._last_diffusion_sites.append((d_i, r_i, c_i))
                                self._last_diffusion_weights.append(diff_k)
                                self._last_diffusion_empty_neighbors.append(empties)

                            if lh_mode == "diffusion_limited" and d_i == 0:
                                if not self.lh_exclude_chemisorption or int(st) != 2:
                                    lh_surface_h_sites.append((d_i, r_i, c_i))
                                    total_diff_surface_lh += diff_k

            if total_des > 0:
                rates["desorption"] = float(total_des)
            if enable_diffusion and diffusion_mode == "explicit" and total_diff > 0:
                rates["diffusion"] = float(total_diff)

        # LH formation
        # Default: LH is ON in gas-kinetic mode, OFF in arrival-mode unless explicitly enabled.
        enable_LH = self._bool_param("enable_LH", False if arrival_mode else True)
        if enable_LH:
            if lh_mode == "diffusion_limited":
                # Fast-mixing approximation: diffusion creates encounter pairs and reaction is effectively prompt.
                if lh_surface_sites_total is None:
                    # Eligible surface sites for LH (exclude chemisorption sites if requested).
                    if self.lh_exclude_chemisorption:
                        lh_surface_sites_total = int(
                            np.sum((self.lattice[0] != None) & (self.site_types[0] != 2))
                        )
                    else:
                        lh_surface_sites_total = int(np.sum(self.lattice[0] != None))
                n_sites = int(lh_surface_sites_total or 0)
                n_h = int(len(lh_surface_h_sites))
                if n_sites > 0 and n_h >= 2 and total_diff_surface_lh > 0:
                    # Approximate encounter probability per hop based on coverage.
                    coverage = float(n_h) / float(max(1, n_sites))
                    z = float(self.simulation_parameters.get("lh_encounter_neighbors", 3.0))
                    p_enc = min(1.0, max(0.0, float(z) * float(coverage)))
                    fac = float(self.simulation_parameters.get("lh_diffusion_factor", 0.5))
                    rates["h2_formation_LH"] = float(fac) * float(total_diff_surface_lh) * float(p_enc)
            else:
                if self.adjacent_h_pairs_count > 0:
                    rates["h2_formation_LH"] = float(h2_formation_lh_rate(surface_temp_k, int(self.adjacent_h_pairs_count)))

        # UV processes (optional)
        if uv_flux_factor > 0:
            uv_photon_flux_total = float(uv_photon_flux["integrated_fuv_photon_flux_photons_cm2_s"]) * float(uv_flux_factor)
            uv_active = self.uv_mode == "continuous"
            uv_h2_mode = str(self.simulation_parameters.get("uv_h2_mode", "adjacent_pair") or "adjacent_pair").strip().lower()
            if self.uv_mode != "continuous":
                base_uv_rate = self.simulation_parameters.get("uv_pulse_start_rate_s", 5.0e-8)
                if self.uv_pulse_enabled:
                    rates["uv_pulse_start"] = float(base_uv_rate) * float(uv_flux_factor)
                uv_active = self.uv_pulse_active

            if uv_active:
                if self.h_atoms_on_surface > 0:
                    rates["uv_photodesorption"] = float(uv_photodesorption_rate(uv_photon_flux_total, self.h_atoms_on_surface))
                if uv_h2_mode == "chemisorption_photofrag":
                    n_chem = int(len(self.occupied_chemisorption_surface))
                    min_h = int(self.simulation_parameters.get("uv_photofrag_min_chemisorbed_h", 2) or 2)
                    h_per_event = float(self.simulation_parameters.get("uv_photofrag_h_per_event", 2.0) or 2.0)
                    sigma_frag = self.simulation_parameters.get("uv_photofrag_cross_section_cm2", None)
                    branch_frag = self.simulation_parameters.get("uv_photofrag_branching_ratio", None)
                    if n_chem >= min_h:
                        active_motifs = max(0.0, float(n_chem - min_h + 1) / max(h_per_event, 1.0))
                        rates["h2_formation_UV"] = float(
                            uv_h2_photofragmentation_rate(
                                uv_photon_flux_total,
                                n_chem,
                                active_motifs,
                                absorption_cross_section_cm2=sigma_frag,
                                branching_ratio=branch_frag,
                            )
                        )
                elif self.adjacent_h_pairs_count > 0:
                    rates["h2_formation_UV"] = float(uv_h2_formation_rate(uv_photon_flux_total, self.adjacent_h_pairs_count))
                uv_diff = float(self.simulation_parameters.get("uv_stimulated_diffusion_factor", 1.0))
                if uv_diff > 1.0 and "diffusion" in rates:
                    rates["uv_stimulated_diffusion"] = float(uv_diff - 1.0) * float(rates["diffusion"])

        return {k: float(v) for k, v in rates.items() if float(v) > 0.0}

    def _execute_arrival(self) -> None:
        self.total_arrivals += 1

        tau = float(self.simulation_parameters.get("beam_dissociation_fraction", 1.0))
        tau = max(0.0, min(1.0, tau))
        is_atom = random.random() < float(tau)

        if is_atom:
            self.total_impinging_h_atoms += 1
        else:
            self.total_impinging_h2_molecules += 1

        # Blocking: coverage-dependent reduction factors.
        theta = self._theta_h2()
        stick_block = float(self.simulation_parameters.get("sticking_blocking_strength", 0.0) or 0.0)
        er_block = float(self.simulation_parameters.get("er_blocking_strength", 0.0) or 0.0)

        if not is_atom:
            # Undissociated H2 arrival: may stick (low T) or bounce (baseline channel).
            if not self.empty_surface:
                # No room: treat as immediate release.
                self.h2_molecules_desorbed_beam += 1
                self.h2_molecules_released_beam += 1
                return

            # Try to stick onto a random empty surface site.
            r, c = random.choice(tuple(self.empty_surface))
            self._handle_new_h2(place_site=(r, c), origin="beam", mechanism="")
            return

        # Atom arrival: attempt ER/abstraction against the chemisorbed reservoir.
        surface_area_cm2 = float(self.total_accessible_surface_sites) * float(
            float(self.simulation_parameters.get("site_area_angstroms_sq", 9.0)) * 1e-16
        )
        sigma = float(self.simulation_parameters.get("er_cross_section_cm2", surface_chemistry_data.get("er_cross_section_cm2", 1e-15)))
        p_er = float(self.simulation_parameters.get("er_reaction_probability", 0.1))
        p_er_eff = max(0.0, min(1.0, float(p_er) * max(0.0, 1.0 - float(er_block) * float(theta))))

        n_targets = float(len(self.occupied_chemisorption_surface))
        p_react = 0.0
        if surface_area_cm2 > 0 and n_targets > 0 and sigma > 0 and p_er_eff > 0:
            p_react = min(1.0, float(sigma) * float(p_er_eff) * float(n_targets) / float(surface_area_cm2))

        if p_react > 0 and random.random() < float(p_react) and self.occupied_chemisorption_surface:
            r, c = random.choice(tuple(self.occupied_chemisorption_surface))
            # Consume one chemisorbed H (incoming atom + surface H -> H2)
            self._desorb_h_at(0, int(r), int(c))
            self.h2_molecules_formed += 1
            self.h2_molecules_formed_ER += 1
            self._handle_new_h2(place_site=(int(r), int(c)), origin="formed", mechanism="ER")
            return

        # Otherwise, try sticking to an empty site.
        if not self.empty_surface:
            return

        # Sticking probability model for arrival mode.
        stick_model = str(self.simulation_parameters.get("sticking_temp_model", "constant") or "constant").strip().lower()
        p_stick = float(self.simulation_parameters.get("sticking_probability", 0.3))
        if stick_model not in {"constant", "none"}:
            # Simple empirical falloff with gas temperature (matching adsorption_rate's convention).
            gas_temp_k = float(self.simulation_parameters.get("gas_temperature_k", 100.0))
            p_stick = float(p_stick) * float(np.exp(-float(gas_temp_k) / 100.0))
        p_stick *= max(0.0, 1.0 - float(stick_block) * float(theta))
        p_stick = max(0.0, min(1.0, float(p_stick)))

        if random.random() < float(p_stick):
            r, c = random.choice(tuple(self.empty_surface))
            self._adsorb_h_at_surface(int(r), int(c))

    def execute_event(self, event_type: str) -> None:
        if event_type == "arrival":
            self._execute_arrival()
            return

        if event_type == "adsorption":
            if not self.empty_surface:
                return
            r, c = random.choice(tuple(self.empty_surface))
            self._adsorb_h_at_surface(int(r), int(c))
            return

        if event_type in {"desorption", "uv_photodesorption"}:
            if not self.occupied:
                return

            if event_type == "uv_photodesorption":
                # Photodesorption is treated as uniform per H atom (rate already scales with count).
                d, r, c = random.choice(tuple(self.occupied))
                self._desorb_h_at(int(d), int(r), int(c))
                return

            # Thermal desorption: choose a site weighted by its per-site desorption rate.
            if self._last_desorption_sites and self._last_desorption_weights:
                d, r, c = random.choices(self._last_desorption_sites, weights=self._last_desorption_weights, k=1)[0]
                self._desorb_h_at(int(d), int(r), int(c))
                return

            # Fallback (e.g., direct unit-test call without calculate_rates): compute weights on the fly.
            T = float(self.simulation_parameters.get("surface_temperature_k", 10.0))
            sites = list(self.occupied)
            weights = [float(h_desorption_rate(float(self.E_bind_eV_map[d, r, c]), T)) for (d, r, c) in sites]
            if not any(w > 0 for w in weights):
                return
            d, r, c = random.choices(sites, weights=weights, k=1)[0]
            self._desorb_h_at(int(d), int(r), int(c))
            return

        if event_type in {"diffusion", "uv_stimulated_diffusion"}:
            if not self.occupied:
                return

            # Diffusion: choose a mobile site weighted by its diffusion rate.
            if self._last_diffusion_sites and self._last_diffusion_weights and self._last_diffusion_empty_neighbors:
                idx = random.choices(
                    range(len(self._last_diffusion_sites)), weights=self._last_diffusion_weights, k=1
                )[0]
                d, r, c = self._last_diffusion_sites[idx]
                empties = self._last_diffusion_empty_neighbors[idx]
                if not empties:
                    return
                nd, nr, nc = random.choice(empties)
            else:
                # Fallback: compute candidates on the fly (e.g., direct unit-test call).
                T = float(self.simulation_parameters.get("surface_temperature_k", 10.0))
                candidates = []
                weights = []
                empty_lists = []
                for d, r, c in self.occupied:
                    empties = [nb for nb in self.get_neighbors_3d(int(d), int(r), int(c)) if self.lattice[nb] == "C"]
                    if not empties:
                        continue
                    rate = float(self._diffusion_rate(int(self.site_types[int(d), int(r), int(c)]), T))
                    if rate <= 0:
                        continue
                    candidates.append((int(d), int(r), int(c)))
                    weights.append(rate)
                    empty_lists.append(empties)
                if not candidates:
                    return
                idx = random.choices(range(len(candidates)), weights=weights, k=1)[0]
                d, r, c = candidates[idx]
                nd, nr, nc = random.choice(empty_lists[idx])

            # Move H
            self.lattice[d, r, c] = "C"
            self.update_adjacent_h_pairs_count(int(d), int(r), int(c), False)
            self.occupied.discard((int(d), int(r), int(c)))
            if int(d) == 0 and self.lattice[0, int(r), int(c)] is not None:
                self.empty_surface.add((int(r), int(c)))
                self.occupied_chemisorption_surface.discard((int(r), int(c)))

            self.lattice[nd, nr, nc] = "H"
            self.update_adjacent_h_pairs_count(int(nd), int(nr), int(nc), True)
            self.occupied.add((int(nd), int(nr), int(nc)))
            if int(nd) == 0:
                self.empty_surface.discard((int(nr), int(nc)))
                if int(self.site_types[0, int(nr), int(nc)]) == 2:
                    self.occupied_chemisorption_surface.add((int(nr), int(nc)))
            return

        if event_type in {"h2_formation_LH", "h2_formation_UV"}:
            lh_mode = str(self.simulation_parameters.get("lh_formation_mode", "pairs") or "pairs").strip().lower()
            uv_h2_mode = str(self.simulation_parameters.get("uv_h2_mode", "adjacent_pair") or "adjacent_pair").strip().lower()
            if event_type == "h2_formation_UV" and uv_h2_mode == "chemisorption_photofrag":
                candidates = [(0, int(r), int(c)) for (r, c) in self.occupied_chemisorption_surface if self.lattice[0, int(r), int(c)] == "H"]
                if len(candidates) < 2:
                    return
                (d1, r1, c1), (d2, r2, c2) = random.sample(candidates, 2)
            elif event_type == "h2_formation_LH" and lh_mode == "diffusion_limited":
                # Pick any two eligible surface H atoms (fast-diffusion/mixing approximation).
                candidates = []
                for d, r, c in self.get_occupied_sites():
                    if int(d) != 0:
                        continue
                    if self.lh_exclude_chemisorption and int(self.site_types[0, int(r), int(c)]) == 2:
                        continue
                    candidates.append((int(d), int(r), int(c)))
                if len(candidates) < 2:
                    return
                (d1, r1, c1), (d2, r2, c2) = random.sample(candidates, 2)
            else:
                if self.adjacent_h_pairs_count <= 0:
                    return
                pairs = []
                for d, r, c in self.get_occupied_sites():
                    for nd, nr, nc in self.get_neighbors_3d(int(d), int(r), int(c)):
                        if (
                            self.lattice[nd, nr, nc] == "H"
                            and (int(d), int(r), int(c)) < (int(nd), int(nr), int(nc))
                            and self._lh_pair_allowed(int(d), int(r), int(c), int(nd), int(nr), int(nc))
                        ):
                            pairs.append(((int(d), int(r), int(c)), (int(nd), int(nr), int(nc))))
                if not pairs:
                    return
                (d1, r1, c1), (d2, r2, c2) = random.choice(pairs)

            # Remove the two H atoms.
            self.lattice[d1, r1, c1] = "C"
            self.update_adjacent_h_pairs_count(d1, r1, c1, False)
            self.occupied.discard((d1, r1, c1))
            if d1 == 0 and self.lattice[0, r1, c1] is not None:
                self.empty_surface.add((r1, c1))
                self.occupied_chemisorption_surface.discard((r1, c1))

            self.lattice[d2, r2, c2] = "C"
            self.update_adjacent_h_pairs_count(d2, r2, c2, False)
            self.occupied.discard((d2, r2, c2))
            if d2 == 0 and self.lattice[0, r2, c2] is not None:
                self.empty_surface.add((r2, c2))
                self.occupied_chemisorption_surface.discard((r2, c2))

            self.h_atoms_on_surface -= 2
            self.h2_molecules_formed += 1
            if event_type == "h2_formation_LH":
                self.h2_molecules_formed_LH += 1
                mech = "LH"
            else:
                self.h2_molecules_formed_UV += 1
                mech = "UV"

            # Place product on surface (if possible) for blocking, otherwise count prompt.
            place = (r1, c1) if d1 == 0 else ((r2, c2) if d2 == 0 else None)
            self._handle_new_h2(place_site=place, origin="formed", mechanism=mech)
            return

        if event_type == "h2_formation_ER":
            if not self.occupied:
                return
            d, r, c = random.choice(tuple(self.occupied))
            self._desorb_h_at(int(d), int(r), int(c))
            self.h2_molecules_formed += 1
            self.h2_molecules_formed_ER += 1
            place = (int(r), int(c)) if int(d) == 0 else None
            self._handle_new_h2(place_site=place, origin="formed", mechanism="ER")
            return

        if event_type == "h2_desorption":
            if self.h2_molecules_on_surface <= 0:
                return
            all_sites = list(self.h2_sites_beam | self.h2_sites_formed)
            if not all_sites:
                return
            r, c = random.choice(all_sites)
            if self.lattice[0, r, c] != "H2":
                # Stale bookkeeping; repair.
                self.h2_sites_beam.discard((r, c))
                self.h2_sites_formed.discard((r, c))
                return

            self.lattice[0, r, c] = "C"
            self.empty_surface.add((r, c))
            self.h2_molecules_on_surface -= 1
            if (r, c) in self.h2_sites_formed:
                self.h2_sites_formed.discard((r, c))
                self.h2_molecules_released_formed += 1
            else:
                self.h2_sites_beam.discard((r, c))
                self.h2_molecules_released_beam += 1
            return

        if event_type == "uv_pulse_start":
            self.uv_pulse_active = True
            self.last_uv_pulse_time = float(self.time)
            return

    def run_gillespie(
        self,
        max_time: float,
        max_steps: Optional[int] = None,
        callback: Optional[Callable[["KineticMonteCarlo", str], None]] = None,
    ):
        step_count = 0

        # Stop criteria that should be interpreted relative to the start of this run call.
        start_arrivals = int(self.total_arrivals)
        max_arrivals = self.simulation_parameters.get("max_arrivals", None)
        try:
            max_arrivals = int(max_arrivals) if max_arrivals is not None else None
        except (TypeError, ValueError):
            max_arrivals = None

        target_exposure = self.simulation_parameters.get("target_exposure_atoms_cm2", None)
        try:
            target_exposure = float(target_exposure) if target_exposure is not None else None
        except (TypeError, ValueError):
            target_exposure = None
        start_impinging_atoms = int(self.total_impinging_h_atoms)

        while float(self.time) < float(max_time):
            if max_steps and step_count >= int(max_steps):
                break

            # UV pulse end condition (simple timer)
            if self.uv_pulse_active and (float(self.time) - float(self.last_uv_pulse_time)) >= float(self.uv_pulse_duration):
                self.uv_pulse_active = False

            # Arrival-limited runs (Grieco harness)
            if max_arrivals is not None and self._arrival_mode_enabled():
                if (int(self.total_arrivals) - int(start_arrivals)) >= int(max_arrivals):
                    break

            # Exposure-limited runs (optional; matches "few × 1e15 atoms/cm2" protocol statement)
            if target_exposure is not None and self._arrival_mode_enabled():
                area_cm2 = float(self.total_accessible_surface_sites) * float(
                    float(self.simulation_parameters.get("site_area_angstroms_sq", 9.0)) * 1e-16
                )
                if area_cm2 > 0:
                    imp = int(self.total_impinging_h_atoms) - int(start_impinging_atoms)
                    exposure = float(imp) / float(area_cm2)
                    if exposure >= float(target_exposure):
                        break

            rates = self.calculate_rates()
            if not rates:
                break
            total_rate = float(sum(rates.values()))
            if total_rate <= 0:
                break

            delta_t = random.expovariate(float(total_rate))
            self.last_delta_t = float(delta_t)

            if float(self.time) + float(delta_t) > float(max_time):
                self.time = float(max_time)
                self._update_surface_temperature_from_ramp()
                break

            self.time = float(self.time) + float(delta_t)
            self._update_surface_temperature_from_ramp()

            if self.debug_print_timestep and (step_count % max(1, int(self.debug_timestep_every)) == 0):
                try:
                    print(
                        f"[KMC] step={step_count} time={self.time:.6e} total_rate={total_rate:.6e} "
                        f"delta_t={delta_t:.6e} next_time={self.time:.6e}"
                    )
                except Exception:
                    pass

            chosen_event = random.choices(list(rates.keys()), weights=list(rates.values()), k=1)[0]
            self.execute_event(str(chosen_event))
            if callback is not None:
                callback(self, str(chosen_event))
            step_count += 1
        return []
