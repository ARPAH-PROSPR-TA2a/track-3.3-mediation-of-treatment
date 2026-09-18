# CALERIE DNAm beta mediation of 24-month MetS

`run_DNAm_MetS_24mo.R` is a top-level R script, following the existing CALERIE
runners. Edit its four paths (`repo`, `treatment_path`, `outcome_path`, `out_dir`)
and source it, or run it section by section. It sets the working directory to
the local GitHub checkout and sources `mediation/main.R` from that checkout.

```r
# From the runner's folder, after editing the paths:
source("run_DNAm_MetS_24mo.R")
```

The script visibly loads the two RDS files, releases their unused coefficient
tables to save memory, checks the intended FUs/strata, calls `FAST_mediation()`,
and saves the results. `treatment_results`, `outcome_results`,
`mediation_results`, and `mediation_summary` remain in the R session to inspect.

## Input paths

The treatment path matches the confirmed production output directory:

```text
~/FAST/Outputs/1.1.1/DNAm_Bvals_1.1.1/DNAm_betas_1.1.1.rds
```

The outcome input is:

```text
~/FAST/Outputs/3.3/DNAm_betas_3.3OWAS_MetS/DNAm_betas_3.3OWAS_MetS_results.rds
```

Both paths must use DNAm beta values and compatible treatment coding and
covariates. FU1 is 12-month DNAm and FU2 is 24-month DNAm; MetS is fixed at
24 months. The MetS analysis requires observed outcome, so its participants can
differ from those in the treatment analysis. Effect summaries alone cannot
verify participant overlap or measurement scale.

The runner requires both analyses (`analysis_change`, `analysis_level`), all
three strata (`all`, `male`, `female`), both FUs, and nonempty filtered probe
sets. The core validates finite estimates, positive SEs, unique keys, and
matching filtered-probe membership. Partial joins stop for review. No upstream
significance screening is applied before mediation; BH is calculated separately
within each analysis/stratum/FU and within the filtered probes.

## Outputs and resources

Default output folder: `~/FAST/Outputs/3.3/DNAm_betas_3.3Med_MetS/`.

- `DNAm_betas_3.3Med_MetS_results.rds`: nested indirect effects and join diagnostics.
- `summary.tsv`: tested/filtered counts and significance counts for all 12
  analysis/stratum/FU combinations.
- `provenance.rds`: input and repository paths, timing, FU definitions, method,
  and R session information.
- `run.log`: timestamped input paths, load times and object sizes, per-FU/stratum
  input and filtered counts, join diagnostics, mediation/save times, the full
  significance summary, and completion.

Progress and the summary print to both the R console and `run.log`. Core
mediation errors are logged; errors elsewhere are shown by R on the console,
with the last log entry identifying the stage. Existing outputs require a new
`out_dir` to avoid overwriting them. The script is tracked under
`mediation/CALERIE/DNAm/MetS_24mo/` in the GitHub repository.

Use one R process: approximately 2 vCPUs and 32 GiB RAM, with 64 GiB for additional
headroom. The full RDS files must each fit in memory while loading; unused
coefficient tables are removed before the next file is loaded. Allow roughly
5 GiB of additional output space. A 15–30 minute initial time allowance is an
estimate; production deserialization has not been measured.

Focused synthetic validation: `Rscript tests/mediation/test_dnam_mediation_runners.R`.
