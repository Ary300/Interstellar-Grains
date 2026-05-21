# Review Triage And Verified Fixes

This note records what was checked directly against the current repository state on 2026-05-21, with emphasis on the highest-impact review concerns.

## Verified in current repo

- `Calibration` is the correct term for the Grieco comparison, not `validation`.
  - The current figure/data pipeline calibrates against the Grieco benchmark and then evaluates goodness-of-fit.
  - The paper should describe `\chi^2_\mathrm{red}` as a calibration-set fit statistic.

- The warm ER plateau is analytic by construction.
  - Above the warm crossover, the model converges to the intensive limit `\epsilon \approx f_\mathrm{chem} P_\mathrm{ER} S`.
  - Any abstract or conclusion statement should explicitly frame the `\epsilon \approx 0.19` plateau as consistency with this limit, not as an independent stochastic discovery.

- The metallicity correction is real and important.
  - Using the representative solar-metallicity timescale table:
    - `T = 60 K, n_H = 10^3 cm^-3`: `t_H2 / t_ff = 0.933` at `Z' = 1`
    - the same scenario becomes `9.33` at `Z' = 0.1`
    - and `93.3` at `Z' = 0.01`
  - Any high-redshift discussion tied to `Z' = 0.01–0.1` must therefore be reframed.

## Not reproduced from current repo outputs

- The free-fall-time complaint in review items `D1/G1` does **not** match the current repository tables/scripts.
  - The current code uses
    - `t_ff = sqrt(3 pi / (32 G rho))`
    - `rho = mu m_H n_H`
    - `mu = 1.4`
  - This gives:
    - `n_H = 100 cm^-3`: `t_ff = 4.35 Myr`
    - `n_H = 10^3 cm^-3`: `t_ff = 1.38 Myr`
    - `n_H = 10^4 cm^-3`: `t_ff = 0.435 Myr`
  - These values are exactly what appear in the current `results/tables/table_timescales_representative.csv`.
  - The reviewer numbers that imply `13.8 Myr` at `n_H = 100 cm^-3` do not match the code path now present in this repository.

- Because the manuscript source is not present in this repository, it is not possible here to determine whether the current paper text/equations still contain stale or inconsistent timescale values copied from an earlier calculation path.

## In-repo fixes applied

- Added metallicity-aware representative-timescale generation in `make_representative_timescale_table.py`.
  - The script now accepts repeatable `--Z` metallicity scaling values.
  - Dust abundance scales linearly with `Z'`, so `sigma_H,eff = Z' sigma_H`, `k_eff \propto Z'`, and `t_H2 \propto 1/Z'`.

## Recommended manuscript wording changes

- Calibration framing:
  - "The model reproduces the Grieco calibration dataset with `\chi^2_\mathrm{red} = ...`; we interpret this as a calibration-set goodness-of-fit statistic rather than independent validation."

- Warm plateau framing:
  - "At `T \gtrsim 150 K`, the simulation converges to the analytic ER limit `\epsilon_\mathrm{high-T} = f_\mathrm{chem} P_\mathrm{ER} S`, giving an efficiency plateau near 0.19."

- High-redshift metallicity framing:
  - "At solar metallicity the representative `60 K, 10^3 cm^-3` case gives `t_H2 / t_ff \approx 0.93`; however, applying the paper's own linear dust-metallicity scaling increases this to `\sim 9.3` at `Z' = 0.1` and `\sim 93` at `Z' = 0.01`, so dust-surface formation alone cannot keep pace with collapse in low-metallicity high-redshift gas."

## Next paper-level fixes once manuscript source is available

- Replace all remaining uses of `validation` with `calibration` where the Grieco dataset is the fitting target.
- Update the abstract/conclusion high-redshift claim to include the metallicity-scaled timescales.
- Audit Equation (13), Appendix E, and any manually entered timescale tables against the current generated CSVs.
- Add a one-sentence note that only the product `f_\mathrm{chem} P_\mathrm{ER}` is constrained by the warm plateau.
