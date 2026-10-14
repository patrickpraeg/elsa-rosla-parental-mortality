suppressPackageStartupMessages({library(tidyverse); library(survival)})

set.seed(321)
B_SUB <- 1000

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

nc_cohort1 <- fit_first_stage(nc %>% filter(cohort1), first_stage_formula("post47", "R_months47"))$data
nc_cohort2 <- fit_first_stage(nc %>% filter(cohort2), first_stage_formula("post72", "R_months72"))$data

rows <- list()
for (cl in cells) {
  d <- get(paste0("nc_", cl$cohort))
  for (sg in subgroups) {
    sub <- d[!is.na(d[[sg$var]]) &
               as.character(d[[sg$var]]) == as.character(sg$value), , drop = FALSE]
    rhs  <- setdiff(c("schlage", "presid", cl$rv, "fem", "par_class2", "nsibs_cat"), sg$var)
    f_sg <- as.formula(paste0("Surv(", cl$time, ", ", cl$event, ") ~ ",
                              paste(rhs, collapse = " + ")))
    m <- tryCatch(coxph(f_sg, data = sub), error = function(e) NULL)
    if (is.null(m)) next

    bs <- run_bootstrap(d, cl$time, cl$event, cl$iv, cl$rv, B = B_SUB,
                        subgroup_var = sg$var, subgroup_value = sg$value)

    rows[[paste(cl$parent, cl$reform, sg$label)]] <- tibble(
      panel = sg$panel,
      outcome = paste0(cl$parent, "s, ", cl$reform, " reform"),
      subgroup = sg$label,
      HR = round(unname(summary(m)$coefficients["schlage", "exp(coef)"]), 4),
      lo = unname(bs["lo"]), hi = unname(bs["hi"]),
      N = m$n, events = m$nevent, n_ok = unname(bs["n_ok"])
    )
  }
}

res <- bind_rows(rows)
cat("\nTable A6, NatCen columns\n")
print(as.data.frame(res), row.names = FALSE)
write_csv(res, file.path(output_path, "tableA6_natcen.csv"))
