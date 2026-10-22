# CALERIE proteomic-clock OutcomeWAS

Use [the INF runner](INF_24mo/run_Proteomic_Clocks_INF_24mo.R) or
[the MetS runner](MetS_24mo/run_Proteomic_Clock_MetS_24mo.R) for the corresponding
fixed 24-month outcome. Both use adjusted clock values produced by Track 1.1.1,
with `omics_type = "other"`.

## Visits and sample matching

The outcome is measured at 24 months. Clock measurements at 3, 6, 12, and
24 months are each paired with the same subject's baseline measurement for
separate change and level analyses. A missing visit does not remove a subject
from the other visits for which a baseline pair is available.

| Clock visit (months) | FU |
| --- | --- |
| 0 (baseline) | 0 |
| 3 | 1 |
| 6 | 2 |
| 12 | 3 |
| 24 | 4 |

These codes match Track 1.1.1's proteomic-clock treatment results. The old
24-month-only runs encoded 24 months as `FU1`; their effects must not be joined
directly to the current Track 1.1.1 results.

The clock CSV is long: `Organ` and `proteomics_barcode` together identify a row,
with `Adjusted_Predicted_Age` holding the value. Multiple organs per barcode are
expected. The runners validate that composite key and align the pivoted sample
columns to the phenotype sample IDs.

## Configured paths

Check the explicit path block in each runner for the target machine:

- Clock values: `~/CALERIE/Proteomic_Clock/data/omics_Proteomic_Clock.csv`.
- Phenotypes: `~/CALERIE/Proteomic_Clock/data/CALERIE_clinical_w_omics_crosswalks.csv`.
- Genetic PCs: `~/CALERIE/data/Proteomics/CALERIE_GeneticPCs_20200727.xlsx`.
- Outcomes: `~/CALERIE/Proteomics/3.3/data/INF_CALERIE_24mo_Outcome.csv` or
  `MetS_CALERIE_24mo_Outcome.csv` in the same directory.
- Pipeline checkout: `~/FAST/GitHub/track-3.3`.

The clock and phenotype defaults match the current Track 1.1.1 clock producer.
They describe configured inputs, not a verification of files on the server.

## Outputs and checkpoints

The runners write under `~/CALERIE/Proteomic_Clock/3.3/`, in
`INF_24mo_allFU_v1` or `MetS_24mo_allFU_v1`. Each contains `results.rds`,
`reports.rds`, `run.log`, `checkpoints/`, and plots under `Figures/change` and
`Figures/level`.

Use these separate directories for the all-visit analyses. Do not copy the old
24-month-only checkpoints into them: their `FU1` has a different meaning.
Choose a fresh directory again if inputs, model settings, or result schema
change. The `OLD_` scripts are retained for reference.

Synthetic runner validation, including manual regression comparisons:

```bash
Rscript tests/outcomewas/test_proteomic_clock_runners.R
```
