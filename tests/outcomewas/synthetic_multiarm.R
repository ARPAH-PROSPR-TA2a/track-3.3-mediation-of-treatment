# Entirely simulated subjects, measurements and outcomes; no example RDS inputs.
.make_multiarm_outcome_fixture <- function(outcome_type = "continuous",
                                          arms = c("0", "2", "7"),
                                          n_per_cell = 16L,
                                          n_analytes = 6L) {
  stopifnot(outcome_type %in% c("continuous", "tte"))
  set.seed(20261015)
  subjects <- expand.grid(
    replicate = seq_len(n_per_cell), FEMALE = 0:1,
    TREATMENT_GROUP = arms, stringsAsFactors = FALSE
  )
  n <- nrow(subjects)
  subjects$SUBJECT_ID <- sprintf("synthetic_subject_%03d", seq_len(n))
  subjects$AGE <- rnorm(n, 50, 8)
  latent <- rnorm(n)
  treatment_shift <- c(0, seq(0.2, 0.5, length.out = length(arms) - 1L))[
    match(subjects$TREATMENT_GROUP, arms)
  ]
  baseline <- matrix(rnorm(n_analytes * n), nrow = n_analytes)
  visits <- lapply(0:2, function(fu) {
    frame <- subjects[c("SUBJECT_ID", "FEMALE", "TREATMENT_GROUP", "AGE")]
    frame$FU <- fu
    frame$SAMPLE_ID <- paste0(frame$SUBJECT_ID, "_FU", fu)
    frame
  })
  measurements <- lapply(0:2, function(fu) {
    if (fu == 0L) return(baseline)
    baseline + matrix(rnorm(n_analytes * n, sd = 0.8), nrow = n_analytes) +
      rep(0.4 * latent + fu * treatment_shift, each = n_analytes)
  })
  signal <- measurements[[2L]][1L, ] - baseline[1L, ]
  linear_predictor <- 0.4 * signal + treatment_shift +
    0.15 * subjects$FEMALE + 0.02 * (subjects$AGE - 50)
  continuous <- linear_predictor + rnorm(n, sd = 1)
  event_time <- rexp(n, rate = 0.25 * exp(0.4 * linear_predictor))
  censor_time <- rexp(n, rate = 0.12)
  pheno <- do.call(rbind, visits)
  subject_index <- match(pheno$SUBJECT_ID, subjects$SUBJECT_ID)
  if (outcome_type == "continuous") {
    pheno$OUTCOME <- continuous[subject_index]
  } else {
    pheno$OUTCOME_TIME <- pmin(event_time, censor_time)[subject_index]
    pheno$OUTCOME_STATUS <- as.integer(event_time <= censor_time)[subject_index]
  }
  pheno$TREATMENT_GROUP <- factor(pheno$TREATMENT_GROUP, levels = rev(arms))
  pheno$FEMALE <- factor(pheno$FEMALE, levels = 0:1)
  pheno$FU <- factor(pheno$FU, levels = 0:2)
  omics <- data.frame(
    ANALYTE_NAME = sprintf("synthetic_analyte_%02d", seq_len(n_analytes)),
    do.call(cbind, measurements), check.names = FALSE
  )
  names(omics)[-1L] <- pheno$SAMPLE_ID
  labels <- setNames(c("Usual care", paste("Intervention", arms[-1L])), arms)
  list(pheno = pheno, omics = omics, treatment_labels = labels)
}

.quiet_multiarm <- function(expression) suppressMessages(suppressWarnings(expression))

# Reference models prepare data directly from the simulated input, independently
# of the pipeline's validation, preparation and coefficient extraction helpers.
.manual_multiarm_outcome <- function(fixture, response, stratum, fu, analyte,
                                     outcome_type = "continuous") {
  pheno <- fixture$pheno
  rows <- as.character(pheno$FU) == as.character(fu)
  if (stratum != "all") {
    rows <- rows & as.character(pheno$FEMALE) == if (stratum == "male") "0" else "1"
  }
  md <- pheno[rows, , drop = FALSE]
  baseline <- pheno[as.character(pheno$FU) == "0", , drop = FALSE]
  baseline <- baseline[match(md$SUBJECT_ID, baseline$SUBJECT_ID), , drop = FALSE]
  values <- fixture$omics[match(analyte, fixture$omics$ANALYTE_NAME), ]
  md$analyte_baseline <- as.numeric(values[, baseline$SAMPLE_ID, drop = TRUE])
  followup <- as.numeric(values[, md$SAMPLE_ID, drop = TRUE])
  md$analyte <- if (response == "change") followup - md$analyte_baseline else followup
  outcome_columns <- if (outcome_type == "continuous") "OUTCOME" else c("OUTCOME_TIME", "OUTCOME_STATUS")
  md <- md[complete.cases(md[c("analyte", "analyte_baseline", "AGE", "TREATMENT_GROUP", "FEMALE", outcome_columns)]), ]
  md$TREATMENT_GROUP <- factor(as.character(md$TREATMENT_GROUP),
                              levels = as.character(sort(unique(as.integer(as.character(md$TREATMENT_GROUP))))))
  md$FEMALE <- factor(md$FEMALE)
  terms <- c("analyte", "analyte_baseline")
  if (nlevels(md$TREATMENT_GROUP) > 1L) terms <- c(terms, "TREATMENT_GROUP")
  if (length(unique(md$FEMALE)) > 1L) terms <- c(terms, "FEMALE")
  terms <- c(terms, "AGE")
  response_term <- if (outcome_type == "continuous") "OUTCOME" else "survival::Surv(OUTCOME_TIME, OUTCOME_STATUS)"
  formula <- as.formula(paste(response_term, "~", paste(terms, collapse = " + ")))
  fit <- if (outcome_type == "continuous") {
    lm(formula, data = md)
  } else {
    survival::coxph(formula, data = md, ties = "efron")
  }
  coefficients <- coef(summary(fit))
  columns <- if (outcome_type == "continuous") {
    c("Estimate", "Std. Error", "Pr(>|t|)")
  } else {
    c("coef", "se(coef)", "Pr(>|z|)")
  }
  setNames(unname(coefficients["analyte", columns]), c("EFFECT_SIZE", "SE", "P_VALUE"))
}
