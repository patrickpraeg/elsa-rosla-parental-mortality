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

sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

cohort1 <- fit_first_stage(sn %>% filter(cohort1), first_stage_formula("post47", "R_months47"))$data
cohort2 <- fit_first_stage(sn %>% filter(cohort2), first_stage_formula("post72", "R_months72"))$data

specs <- list(
  list(lab = "Father, 1947", d = cohort1, t = "flongevity2", e = "fdead", rv = "R_months47"),
  list(lab = "Father, 1972", d = cohort2, t = "flongevity2", e = "fdead", rv = "R_months72"),
  list(lab = "Mother, 1947", d = cohort1, t = "mlongevity2", e = "mdead", rv = "R_months47"),
  list(lab = "Mother, 1972", d = cohort2, t = "mlongevity2", e = "mdead", rv = "R_months72"))

# Supplement C.1: time-varying coefficient
cat("\nSupplement C.1: time-varying coefficient, beta(t) = b0 + b1*log(t)\n")
tv <- list()
for (s in specs) {
  f <- as.formula(paste0("Surv(", s$t, ", ", s$e, ") ~ schlage + tt(schlage) + presid + ",
                         s$rv, " + fem + par_class2 + nsibs_cat"))
  m <- coxph(f, data = s$d, tt = function(x, t, ...) x * log(t))
  cf <- summary(m)$coefficients
  tv[[s$lab]] <- tibble(
    model      = s$lab,
    b_main     = round(cf["schlage", "coef"], 4),
    p_main     = signif(cf["schlage", "Pr(>|z|)"], 3),
    b_tt       = round(cf["tt(schlage)", "coef"], 4),
    p_tt       = signif(cf["tt(schlage)", "Pr(>|z|)"], 3),
    HR_at_70   = round(exp(cf["schlage","coef"] + cf["tt(schlage)","coef"]*log(70)), 3),
    HR_at_85   = round(exp(cf["schlage","coef"] + cf["tt(schlage)","coef"]*log(85)), 3))
}
tv_df <- bind_rows(tv)
print(as.data.frame(tv_df))

# Supplement C.1: IV estimates within age-at-death strata
cat("\nSupplement C.1: IV estimates within age-at-death strata\n")
st <- list()
for (s in specs) {
  for (br in list(c(20, 75), c(75, 101))) {
    d <- s$d %>% filter(.data[[s$t]] >= br[1], .data[[s$t]] < br[2])
    if (sum(d[[s$e]], na.rm = TRUE) < 50) next
    m <- tryCatch(coxph(as.formula(paste0("Surv(", s$t, ", ", s$e, ") ~ schlage + presid + ",
                                          s$rv, " + fem + par_class2 + nsibs_cat")),
                        data = d), error = function(e) NULL)
    if (is.null(m)) next
    ci <- summary(m)$conf.int["schlage", ]
    st[[paste(s$lab, br[1])]] <- tibble(
      model = s$lab, stratum = paste0("age ", br[1], "-", br[2] - 1),
      HR = round(unname(ci[1]), 3), lo = round(unname(ci[3]), 3),
      hi = round(unname(ci[4]), 3), N = m$n, events = m$nevent)
  }
}
strata_df <- bind_rows(st)
print(as.data.frame(strata_df))

# Supplement C.1: piecewise-constant coefficient by age-at-death band
cat("\nSupplement C.1: IV estimates by age-at-death band, likelihood-ratio test of a constant coefficient\n")
bands <- list()
for (s in specs) {
  d <- s$d %>%
    filter(!is.na(.data[[s$t]]), !is.na(.data[[s$e]])) %>%
    transmute(time = .data[[s$t]], ev = .data[[s$e]], schlage, presid,
              rv = .data[[s$rv]], fem, par_class2, nsibs_cat) %>%
    filter(complete.cases(.), time > 0)
  sp <- survSplit(Surv(time, ev) ~ ., data = d, cut = c(65, 85), episode = "band")
  sp$band <- factor(sp$band, levels = 1:3, labels = c("<65", "65-84", "85+"))
  rhs <- "presid + rv + fem + par_class2 + nsibs_cat"
  m0  <- coxph(as.formula(paste("Surv(tstart, time, ev) ~ schlage +", rhs)), data = sp)
  m1  <- coxph(as.formula(paste("Surv(tstart, time, ev) ~ schlage:band +", rhs)), data = sp)
  lr  <- anova(m0, m1)
  cf  <- summary(m1)$conf.int
  rw  <- grep("^schlage:band", rownames(cf))
  bn  <- sub("schlage:band", "", rownames(cf)[rw])
  bands[[s$lab]] <- tibble(
    model = s$lab, band = bn,
    deaths = as.integer(tapply(sp$ev, sp$band, sum)[bn]),
    HR = round(cf[rw, "exp(coef)"], 3), lo = round(cf[rw, "lower .95"], 3),
    hi = round(cf[rw, "upper .95"], 3),
    p_LR_constant = signif(lr$`Pr(>|Chi|)`[2], 3))
}
bands_df <- bind_rows(bands)
print(as.data.frame(bands_df))

write_csv(tv_df,     file.path(output_path, "ph_timevarying.csv"))
write_csv(strata_df, file.path(output_path, "ph_strata.csv"))
write_csv(bands_df,  file.path(output_path, "ph_bands.csv"))
