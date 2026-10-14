suppressPackageStartupMessages(library(tidyverse))

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

# Table A2
father_analytic <- make_analytical(sn) %>% filter(!is.na(flongevity2), !is.na(fdead))
mother_analytic <- make_analytical(sn) %>% filter(!is.na(mlongevity2), !is.na(mdead))

row_vars <- tribble(
  ~label,                        ~var,          ~show_father, ~show_mother,
  "Father longevity",            "flongevity2",  TRUE,  FALSE,
  "Mother longevity",            "mlongevity2",  FALSE, TRUE,
  "Mother's still alive",        "malive",       FALSE, TRUE,
  "Father's still alive",        "falive",       TRUE,  FALSE,
  "Age finished education",      "schlage",      TRUE,  TRUE,
  "Left school at age 16",       "schlage16",    TRUE,  TRUE,
  "Left school at age 15",       "schlage15",    TRUE,  TRUE,
  "Woman",                       "fem",          TRUE,  TRUE,
  "Number of living siblings",   "nsibs",        TRUE,  TRUE,
  "Ethnic origin: non-white",    "nonwhite",     TRUE,  TRUE,
  "Year of birth",               "byear",        TRUE,  TRUE,
  "Upper class",                 "upper_class",  TRUE,  TRUE,
  "Working class",               "working_class", TRUE, TRUE
)

stats3 <- function(x) {
  c(mean = mean(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
}

class_dummy <- function(df, code) as.numeric(as.character(df$par_class2) == as.character(code))

build_row <- function(label, var, show_father, show_mother) {
  get_col <- function(df, show) {
    if (!show) return(c(mean = NA_real_, min = NA_real_, max = NA_real_))
    if (var == "upper_class")  return(stats3(class_dummy(df, 1)))
    if (var == "working_class") return(stats3(class_dummy(df, 3)))
    stats3(df[[var]])
  }
  initial <- get_col(sn, TRUE)
  father  <- get_col(father_analytic, show_father)
  mother  <- get_col(mother_analytic, show_mother)
  tibble(
    variable      = label,
    initial_mean  = initial["mean"], initial_min = initial["min"], initial_max = initial["max"],
    father_mean   = father["mean"],  father_min  = father["min"],  father_max  = father["max"],
    mother_mean   = mother["mean"],  mother_min  = mother["min"],  mother_max  = mother["max"]
  )
}

tableA2 <- purrr::pmap_dfr(row_vars, build_row)

n_row <- tibble(
  variable = "N",
  initial_mean = nrow(sn), initial_min = NA_real_, initial_max = NA_real_,
  father_mean  = nrow(father_analytic), father_min = NA_real_, father_max = NA_real_,
  mother_mean  = nrow(mother_analytic), mother_min = NA_real_, mother_max = NA_real_
)

tableA2 <- bind_rows(tableA2, n_row)

cat("\nTable A2\n")
print(as.data.frame(tableA2), digits = 4)
write_csv(tableA2, file.path(output_path, "tableA2.csv"))

# Table A3
a2_column <- function(ds) {
  m <- make_analytical(ds) %>% filter(!is.na(mlongevity2), !is.na(mdead))
  pc <- as.character(m$par_class2)
  c(`Mother longevity`               = mean(m$mlongevity2),
    `Mother's still alive`           = mean(m$malive, na.rm = TRUE),
    `Father's still alive`           = mean(m$falive, na.rm = TRUE),
    `Age finished education`         = mean(m$schlage),
    `Left school at age 16`          = mean(m$schlage16),
    `Left school at age 15`          = mean(m$schlage15),
    `Woman`                          = mean(m$fem),
    `Upper class`                    = mean(pc == "1"),
    `Working class`                  = mean(pc == "3"),
    `Ethnic origin: non-white`       = mean(m$nonwhite),
    `Number of siblings`             = mean(m$nsibs),
    `Year of birth`                  = mean(m$byear),
    `In 1947 cohort, +/-60 months`   = mean(m$cohort1),
    `In 1972 cohort, +/-60 months`   = mean(m$cohort2),
    N                                = nrow(m))
}

tableA3 <- tibble(variable = names(a2_column(nc)),
                  NatCen = a2_column(nc),
                  UK_Data_Service = a2_column(sn))

cat("\nTable A3\n")
print(as.data.frame(tableA3 %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
write_csv(tableA3, file.path(output_path, "tableA3.csv"))
