# Mediation Example

This example uses the compact, balanced Track 1.1.1 proteomics fixture as the
shared input to both upstream pipelines:

- 48 subjects;
- follow-ups 0, 1, and 2 for every subject;
- 12 subjects in each treatment-by-sex cell;
- 60 proteomic analytes;
- no missing phenotype or omics values.

`pheno_base.rds` and `proteomics_log2.rds` are copied from the tracked Track
1.1.1 example data. `generate_example_data.R` adds a deterministic continuous
outcome based on the standardized FU1 change in `PROT_00001`, treatment, BMI,
and seeded Gaussian noise. The outcome is only intended to create a stable
pipeline test signal.

From the Track 3.3 repository root:

```bash
Rscript mediation/Examples/generate_example_data.R
Rscript mediation/Examples/generate_example_inputs.R
Rscript mediation/Examples/run_mediation_example.R
```

`generate_example_inputs.R` uses the sibling `../Track1.1.1` checkout by
default. Set `TRACK111_DIR` when it lives elsewhere:

```bash
TRACK111_DIR=/path/to/Track1.1.1 \
  Rscript mediation/Examples/generate_example_inputs.R
```

The committed input result objects let the final mediation example run without
the Track 1.1.1 checkout. The live input-generation test regenerates them in a
temporary directory and checks current upstream compatibility when that
checkout is available.

## Synthetic multi-arm example

The standalone example creates three treatment arms, fits joint alpha models,
then runs the common-beta OutcomeWAS model and passes its results into
mediation:

```bash
Rscript mediation/Examples/run_multiarm_example.R /tmp/fast_multiarm_mediation
```

The optional argument selects an output directory; without it the script uses
R's temporary directory. All data are simulated, and the script needs no
Track 1.1.1 checkout or production data. Its alpha fit demonstrates the same
joint active-versus-control contract rather than splitting the study into
separate two-group fits.

Inspect the preserved treatment contrasts and labels, the common beta reused
for each analyte's indirect effects, and the saved provenance. The products
`alpha_1 * beta` and `alpha_3 * beta` differ because the interventions have
different alpha effects. No treatment interactions or model-selection option
are needed; the outcome model is the same as in the binary example.
