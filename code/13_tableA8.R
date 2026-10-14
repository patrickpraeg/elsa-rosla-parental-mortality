suppressPackageStartupMessages({library(tidyverse); library(survival); library(rdrobust)})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir  <- get_script_dir()
pkg_root    <- dirname(script_dir)
tempdata    <- file.path(pkg_root, "data/tempdata")
output_path <- file.path(pkg_root, "Outputs")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
source(file.path(script_dir, "helpers.R"))

add_reduced_dummies <- function(df) {
  df %>%
    mutate(
      par_class2_int = as.integer(as.character(par_class2)),
      par_class2_d1  = as.numeric(par_class2_int == 1),
      ngor_d1        = as.numeric(ngor == 1)
    )
}

nc <- add_reduced_dummies(prepare_df(readRDS(file.path(tempdata, "ELSA_bmonth-year.rds"))))
sn <- add_reduced_dummies(prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds"))))

reduced_controls <- c("fem", "nonwhite", "par_class2_d1", "ngor_d1", "nsibs")

# First step: Cox model on school-leaving age, predicted linear predictor
fit_step1 <- function(data, time_var, status_var) {
  f <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ factor(schlage)"))
  coxph(f, data = data, id = idauniq, x = TRUE, model = TRUE, na.action = na.omit)
}

cox_pred_full <- function(data, fit) {
  used_rows <- as.integer(rownames(model.frame(fit)))
  pred <- rep(NA_real_, nrow(data))
  pred[used_rows] <- predict(fit, type = "lp")
  pred
}

nc$yhat_f <- cox_pred_full(nc, fit_step1(nc, "flongevity2", "fdead"))
nc$yhat_m <- cox_pred_full(nc, fit_step1(nc, "mlongevity2", "mdead"))
sn$yhat_f <- cox_pred_full(sn, fit_step1(sn, "flongevity2", "fdead"))
sn$yhat_m <- cox_pred_full(sn, fit_step1(sn, "mlongevity2", "mdead"))

# Second step: fuzzy regression discontinuity on the predicted hazard
run_rd <- function(data, y, x, fuzzy, covs) {
  vars <- unique(c(y, x, fuzzy, covs))
  keep <- complete.cases(data[, vars, drop = FALSE])
  rdrobust(y = data[[y]][keep], x = data[[x]][keep], c = 0, covs_drop = TRUE,
           fuzzy = data[[fuzzy]][keep], covs = data[keep, covs, drop = FALSE])
}

datasets <- list(
  list(name = "NatCen",       df = nc),
  list(name = "SN 5050/8346", df = sn)
)

parents <- list(
  list(parent = "Father", y = "yhat_f", xcov = "nmalive_d"),
  list(parent = "Mother", y = "yhat_m", xcov = "nfalive_d")
)

reforms <- list(
  list(reform = "1947", x = "R_months47", fuzzy = "schlage15"),
  list(reform = "1972", x = "R_months72", fuzzy = "schlage16")
)

rows <- list()
for (dd in datasets) {
  for (pp in parents) {
    for (rr in reforms) {
      res <- run_rd(dd$df, y = pp$y, x = rr$x, fuzzy = rr$fuzzy,
                    covs = c(reduced_controls, pp$xcov))
      rows[[paste(dd$name, pp$parent, rr$reform)]] <- tibble(
        dataset     = dd$name,
        parent      = pp$parent,
        reform      = rr$reform,
        coef_conv   = round(res$coef["Conventional", 1], 4),
        se_conv     = round(res$se["Conventional", 1], 4),
        p_conv      = signif(res$pv["Conventional", 1], 3),
        p_robust    = signif(res$pv["Robust", 1], 3),
        ci_lo_rob   = round(res$ci["Robust", 1], 4),
        ci_hi_rob   = round(res$ci["Robust", 2], 4),
        HR          = round(exp(res$coef["Conventional", 1]), 3),
        HR_lo_rob   = round(exp(res$ci["Robust", 1]), 3),
        HR_hi_rob   = round(exp(res$ci["Robust", 2]), 3),
        bw          = round(res$bws["h", 1], 2),
        n_eff_left  = res$N_h[1],
        n_eff_right = res$N_h[2]
      )
    }
  }
}

rd_tab <- bind_rows(rows)
cat("\nTable A8, predicted-hazard RDD column\n")
print(as.data.frame(rd_tab))
write_csv(rd_tab, file.path(output_path, "tableA8.csv"))
