# Reanalysis Closeout (Analysis Only)

This note closes the remaining analysis-side review questions using the current repository outputs on 2026-05-21.

## D1: `t_H2 / t_ff` inconsistency

Status: `resolved in current repo outputs`

The representative solar-metallicity timescale table now gives

- `T = 60 K`
- `n_H = 10^3 cm^-3`
- `Z' = 1`
- `t_H2 / t_ff = 0.9331976868624089`

from:

- `results/tables/table_timescales_representative.csv`
- `results/tables/table_timescales_representative_metallicity.csv`

So, in the current generated analysis products, the `0.933` value is the live value and the stale `0.38` value is not supported.

## G1: free-fall time factor-of-two complaint

Status: `not reproduced`

The current repo uses

- `t_ff = sqrt(3 pi / (32 G rho))`
- `rho = mu m_H n_H`
- `mu = 1.4`

and the resulting free-fall times are internally consistent with the generated CSVs:

- `n_H = 100 cm^-3`: `t_ff = 4.3488 Myr`
- `n_H = 10^3 cm^-3`: `t_ff = 1.3752 Myr`
- `n_H = 10^4 cm^-3`: `t_ff = 0.43488 Myr`

## D6: `k_f` at 20 K

Status: `resolved in favor of the repo-backed simulation conversion`

The stale manuscript-side values `1.0 x 10^-17` and `1.5 x 10^-17 cm^3 s^-1` are both inconsistent with the current repo outputs.

Using the simulation output row at

- `T = 20 K`
- `n_H = 100 cm^-3`
- `epsilon = 0.28515`
- `R_area = 3.198596959998044e5 cm^-2 s^-1`

the repo convention gives

- `k_eff = 4 sigma_H R_area / n_H`
- with `sigma_H = 1e-21 cm^2`
- therefore `k_eff = 1.2794387839992175e-17 cm^3 s^-1`

This is the correct code-backed value for the current repository convention.

The reviewer's `1.03 x 10^-17` value is reproduced only if the factor of 4 is omitted. Under the current repo convention that factor is required because the stored release rate is per unit grain surface area while `sigma_H` is the total geometric cross-section per H nucleus.

## D9: factor-of-4 convention

Status: `resolved analytically`

The current repo uses:

- `R_vol(H2) = 4 sigma_H n_H R_area`
- `k_eff = R_vol / n_H^2 = 4 sigma_H R_area / n_H`
- `t_H2 ~= 1 / (8 sigma_H R_area)`

So the factor of 4 is not an error; it is part of the documented surface-area conversion convention used by:

- `compute_ism_timescales.py`
- `make_representative_timescale_table.py`

## D16: density-dependent enhancement significance

Status: `supported strongly by saved aggregate outputs`

Using the four saved density means at `T = 100 K`, `G0 = 0`:

- `n_H = 10`: `epsilon = 0.215665`, `CI95 = 0.001574`
- `n_H = 100`: `epsilon = 0.226280`, `CI95 = 0.001416`
- `n_H = 10^3`: `epsilon = 0.239630`, `CI95 = 0.000994`
- `n_H = 10^4`: `epsilon = 0.250165`, `CI95 = 0.001320`

Key effect-size checks:

- absolute enhancement from `10` to `10^4 cm^-3`: `0.0345`
- ratio to the largest reported `CI95`: `21.9`

A weighted linear trend of `epsilon` versus `log10(n_H)` over those four means gives:

- slope: `0.0117423` per density decade
- slope SE: `0.0003246`
- `z`-like statistic: `36.18`
- `R^2 = 0.9976`

This is not a substitute for a raw-realisation nonparametric test across all four densities, because the archived deep raw campaign is only available for two densities. But it does show that the reported 16 per cent effect is far larger than the saved statistical uncertainty floor.

## D4 / D10: mechanism of the 100 K enhancement

Status: `still open`

The current saved aggregate outputs include:

- total surface H count
- LH and ER formed counts

but do **not** include separate saved `N_phys` and `N_chem` inventories for the four-density `100 K` comparison. Because of that, the detailed mechanistic attribution remains provisional.

What can be said safely from the saved outputs:

- the density enhancement itself is robust
- the total surface H inventory decreases with increasing density
- therefore the mechanism should not be phrased as simply "more total surface H at higher density"

## Bottom line

Analysis-side closeout status:

- `D1`: resolved
- `G1`: resolved
- `D6`: resolved
- `D9`: resolved
- `D16`: effectively supported from saved outputs
- `D4/D10`: still mechanistically open unless separate physisorbed / chemisorbed inventories are extracted from trajectory-level state
