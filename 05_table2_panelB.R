suppressPackageStartupMessages({library(tidyverse); library(survival)})

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

nc <- prepare_df(readRDS(file.path(tempdata, "ELSA_bmonth-year.rds")))
sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

bandwidth <- 24

# Reduced-form Cox model at a 24-month bandwidth
rows <- list()
for (nm in c("sn", "nc")) {
  ds       <- get(nm)
  ds_label <- if (nm == "sn") "SN 5050/8346" else "NatCen"

  for (cl in cells) {
    d <- ds %>%
      filter(.data[[cl$cohort]]) %>%
      filter(abs(.data[[cl$rv]]) <= bandwidth)

    f <- as.formula(paste0("Surv(", cl$time, ", ", cl$event, ") ~ ", cl$iv,
                           " + ", cl$rv, " + fem + par_class2 + nsibs_cat"))
    m <- coxph(f, data = d)

    s  <- summary(m)
    cf <- s$coefficients[cl$iv, ]
    ci <- s$conf.int[cl$iv, ]

    rows[[paste(nm, cl$parent, cl$reform)]] <- tibble(
      dataset   = ds_label,
      parent    = cl$parent,
      reform    = cl$reform,
      bandwidth = bandwidth,
      HR        = round(unname(cf["exp(coef)"]), 4),
      lo        = round(unname(ci["lower .95"]), 4),
      hi        = round(unname(ci["upper .95"]), 4),
      p         = signif(unname(cf["Pr(>|z|)"]), 3),
      N         = m$n,
      events    = m$nevent,
      AIC       = round(AIC(m), 1),
      BIC       = round(BIC(m), 1)
    )
  }
}

itt <- bind_rows(rows)

cat("\nTable 2, Panel B (SN 5050/8346); Table A5, Panel B (both datasets)\n")
print(as.data.frame(itt))
write_csv(itt, file.path(output_path, "table2_panelB.csv"))
