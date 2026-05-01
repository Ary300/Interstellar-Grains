"""
Cazaux & Tielens analytic H2 formation efficiency prescriptions.

This module originally implemented the reduced CT02 helper from:
  - Cazaux & Tielens (2002), ApJL 575, L29, Eqs. (15)–(18)

It now also exposes an erratum-inspired comparison variant based on:
  - Cazaux & Tielens (2010), ApJ 715, 698, erratum to the 2004 paper

Important note:
The 2010 erratum corrects the transmission coefficients / α_pc. The exact
barrier-level corrected model requires additional potential-curve parameters
from the 2004 treatment that are not represented in this reduced CT02 helper.
Accordingly, `erratum2010_approx` is an explicit approximation calibrated to
the erratum statement that efficiencies above ~25 K are about 3.5× larger.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional


@dataclass(frozen=True)
class CT02Params:
    # Table 1 (silicate surface) in the letter.
    E_H2_K: float = 320.0
    E_HP_K: float = 600.0
    E_HC_K: float = 10000.0
    E_S_K: float = 200.0
    mu: float = 0.005
    nu_H2_s: float = 3.0e12
    nu_HC_s: float = 1.3e13
    erratum_gain_hot: float = 3.5
    erratum_transition_K: float = 25.0


def _beta(nu_s: float, E_K: float, T_K: float) -> float:
    # β = ν exp(-E/kT), with E already in Kelvin units (E/k_B).
    if T_K <= 0:
        return 0.0
    return float(nu_s) * math.exp(-float(E_K) / float(T_K))


def _term_ratio(params: CT02Params) -> float:
    # Common geometric factor in Eqs. (16)–(17):
    #   (1 + sqrt((E_HC - E_S)/(E_HP - E_S)))^2
    num = float(params.E_HC_K) - float(params.E_S_K)
    den = float(params.E_HP_K) - float(params.E_S_K)
    if den <= 0 or num <= 0:
        return float("nan")
    return (1.0 + math.sqrt(num / den)) ** 2


def beta_hp_over_alpha_pc(T_K: float, params: Optional[CT02Params] = None) -> float:
    """
    Eq. (17): β_HP / α_pc (dimensionless).
    """
    params = params or CT02Params()
    g2 = _term_ratio(params)
    if not math.isfinite(g2) or T_K <= 0:
        return float("nan")
    return 0.25 * g2 * math.exp(-float(params.E_S_K) / float(T_K))


def xi_correction(T_K: float, F_ML_s: float, params: Optional[CT02Params] = None) -> float:
    """
    Eq. (16): high-temperature correction factor ξ(T, F).
    """
    params = params or CT02Params()
    if T_K <= 0 or F_ML_s <= 0:
        return 0.0
    g2 = _term_ratio(params)
    if not math.isfinite(g2):
        return float("nan")
    num = float(params.nu_HC_s) * math.exp(-1.5 * float(params.E_HC_K) / float(T_K)) * g2
    return 1.0 / (1.0 + (num / (2.0 * float(F_ML_s))))


def epsilon_ct02_original(T_K: float, F_ML_s: float, params: Optional[CT02Params] = None) -> float:
    """
    Eq. (15): ε_H2(T, F).

    Parameters
    ----------
    T_K : float
        Grain surface temperature (K).
    F_ML_s : float
        Impingement flux in monolayers per second (arrivals per site per second).
    params : CT02Params
        Parameter set (defaults to Table 1 silicate values).
    """
    params = params or CT02Params()
    if T_K <= 0 or F_ML_s <= 0:
        return 0.0

    beta_h2 = _beta(params.nu_H2_s, params.E_H2_K, T_K)
    if beta_h2 <= 0:
        return 0.0

    ratio = beta_hp_over_alpha_pc(T_K, params=params)
    xi = xi_correction(T_K, F_ML_s, params=params)

    denom = 1.0 + (float(params.mu) * float(F_ML_s)) / (2.0 * float(beta_h2)) + float(ratio)
    if denom <= 0:
        return 0.0
    eps = (1.0 / denom) * float(xi)
    return max(0.0, min(1.0, float(eps)))


def epsilon_ct10_erratum_approx(T_K: float, F_ML_s: float, params: Optional[CT02Params] = None) -> float:
    """
    Approximate 2010-erratum correction to the reduced CT02 helper.

    The erratum states that for grain temperatures above ~25 K the corrected
    efficiencies are about 3.5 times higher. We preserve the CT02 functional
    form and apply that published high-temperature gain as an explicit
    approximation.
    """
    params = params or CT02Params()
    if T_K <= 0 or F_ML_s <= 0:
        return 0.0

    beta_h2 = _beta(params.nu_H2_s, params.E_H2_K, T_K)
    if beta_h2 <= 0:
        return 0.0

    gain = 1.0 if float(T_K) <= float(params.erratum_transition_K) else float(params.erratum_gain_hot)
    ratio = beta_hp_over_alpha_pc(T_K, params=params) / float(gain)
    xi = xi_correction(T_K, F_ML_s, params=params)

    denom = 1.0 + (float(params.mu) * float(F_ML_s)) / (2.0 * float(beta_h2)) + float(ratio)
    if denom <= 0:
        return 0.0
    eps = (1.0 / denom) * float(xi)
    return max(0.0, min(1.0, float(eps)))


def epsilon_ct02(
    T_K: float,
    F_ML_s: float,
    params: Optional[CT02Params] = None,
    variant: str = "erratum2010_approx",
) -> float:
    """
    Convenience wrapper for analytic comparison variants.
    """
    key = str(variant or "erratum2010_approx").strip().lower()
    if key in {"original", "ct02", "ct2002"}:
        return epsilon_ct02_original(T_K, F_ML_s, params=params)
    if key in {"erratum2010_approx", "ct10", "corrected", "erratum"}:
        return epsilon_ct10_erratum_approx(T_K, F_ML_s, params=params)
    raise ValueError(f"Unknown CT variant: {variant}")


def rate_ct02_per_area_cm2_s(
    *,
    gas_flux_cm2_s: float,
    sticking: float,
    T_surf_K: float,
    site_area_angstrom2: float,
    params: Optional[CT02Params] = None,
    variant: str = "erratum2010_approx",
) -> float:
    """
    CT02 Eq. (18) re-expressed as a per-grain *surface area* rate (cm^-2 s^-1).

    We use:
      F_ML_s = gas_flux_cm2_s * site_area_cm2
    and:
      R_area = 0.5 * gas_flux_cm2_s * sticking * ε_H2(T, F)

    Returns H2 molecules released to gas per cm^2 of grain surface per second.
    """
    params = params or CT02Params()
    site_area_cm2 = float(site_area_angstrom2) * 1.0e-16
    F = float(gas_flux_cm2_s) * float(site_area_cm2)
    eps = epsilon_ct02(T_surf_K, F, params=params, variant=variant)
    return 0.5 * float(gas_flux_cm2_s) * float(sticking) * float(eps)
