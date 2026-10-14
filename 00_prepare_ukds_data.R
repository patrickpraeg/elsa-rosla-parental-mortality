suppressPackageStartupMessages({library(dplyr); library(tidyr); library(haven)})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
pkg_root   <- dirname(script_dir)

elsa_8346 <- file.path(pkg_root, "data/elsa sn 8346/UKDA-8346-stata/stata/stata13_se")
elsa_5050 <- file.path(pkg_root, "data/elsa sn 5050/UKDA-5050-stata/stata/stata13_se")
tempdata  <- file.path(pkg_root, "data/tempdata")
dir.create(tempdata, showWarnings = FALSE, recursive = TRUE)

core_files <- c(
  file.path(elsa_5050, c("wave_1_core_data_v3.dta",
                         "wave_2_core_data_v4.dta",
                         "wave_3_elsa_data_v4.dta")),
  file.path(elsa_8346, c("wave_4_elsa_data_sl_protect.dta",
                         "wave_5_elsa_data_sl_protect.dta",
                         "wave_6_elsa_data_sl_protect.dta",
                         "wave_7_elsa_data_sl_protect.dta",
                         "wave_8_elsa_data_sl_v3_protect.dta",
                         "wave_9_elsa_data_sl_v2_protect.dta",
                         "wave_10_elsa_data_sl_v4_protect.dta",
                         "wave_11_elsa_data_sl_protect.dta"))
)
ifs_file <- function(i) file.path(elsa_5050, sprintf("wave_%d_ifs_derived_variables.dta", i))

read_lower <- function(path) {
  df <- read_dta(path)
  names(df) <- tolower(names(df))
  df %>% mutate(across(where(haven::is.labelled), haven::zap_labels))
}

# Wave-level extraction
ifs_vars <- c("idauniq", "malive", "mothage", "magedied", "falive", "fathage",
              "fagedied", "edend", "sex", "difjobm", "nonwhite", "nsibs", "wgt")

wave_long <- list()
dob       <- list()
for (i in 1:11) {
  core <- read_lower(core_files[i])

  d <- read_lower(ifs_file(i)) %>%
    select(any_of(ifs_vars)) %>%
    rename(mage = mothage, fage = fathage) %>%
    mutate(wave = i)

  core_extra <- core %>% select(idauniq)
  core_extra$gor <- if ("gor" %in% names(core)) as.character(core$gor) else NA_character_
  wave_long[[i]] <- d %>% left_join(core_extra, by = "idauniq")

  if (i >= 4) dob[[i]] <- core %>% select(idauniq, dobmonth, dobyear)
}

dob <- bind_rows(dob) %>%
  filter(!(is.na(dobmonth) & is.na(dobyear))) %>%
  distinct(idauniq, .keep_all = TRUE)

value_vars <- c("malive", "mage", "magedied", "falive", "fage", "fagedied",
                "edend", "sex", "difjobm", "gor", "nonwhite", "nsibs", "wgt")

wide <- bind_rows(wave_long) %>%
  pivot_wider(id_cols = idauniq, names_from = wave, values_from = all_of(value_vars),
              names_glue = "{.value}{wave}")

df <- wide %>% left_join(dob, by = "idauniq")

# Region
ons_to_letter <- c(E12000001 = "A", E12000002 = "B", E12000003 = "D", E12000004 = "E",
                   E12000005 = "F", E12000006 = "G", E12000007 = "H", E12000008 = "J",
                   E12000009 = "K")

for (i in 1:11) {
  v <- paste0("gor", i)
  x <- as.character(df[[v]])
  recoded <- ons_to_letter[x]
  x <- ifelse(!is.na(recoded), unname(recoded), x)
  x[x %in% c("W99999999", "NA")] <- NA
  df[[v]] <- x
}

gor_cols <- paste0("gor", 1:11)
df$gor <- apply(df[gor_cols], 1, function(row) {
  v <- row[!is.na(row) & row != ""]
  if (length(v)) v[[1]] else NA_character_
})

region_lookup <- c(A = 2, B = 3, D = 4, E = 5, F = 6, G = 7, H = 8, J = 9, K = 10)
df$ngor <- unname(region_lookup[df$gor])

# Parents' survival
first_valid <- function(cols, invalid = c(-1, -2, -8, -9)) {
  m <- as.matrix(cols)
  storage.mode(m) <- "double"
  m[m %in% invalid] <- NA
  n <- nrow(m); out <- rep(NA_real_, n); resolved <- rep(FALSE, n)
  for (j in seq_len(ncol(m))) {
    take <- !resolved & !is.na(m[, j])
    out[take] <- m[take, j]
    resolved <- resolved | take
  }
  out
}

fwd  <- function(prefix) df[paste0(prefix, 1:11)]
bwd  <- function(prefix) df[paste0(prefix, 11:1)]

df$fagedied <- first_valid(fwd("fagedied"))
df$magedied <- first_valid(fwd("magedied"))

df$falive <- first_valid(bwd("falive"))
df$malive <- first_valid(bwd("malive"))

age_with_gate <- function(prefix, alive) {
  primary  <- first_valid(bwd(paste0(prefix, "age")))
  primary[is.na(alive) | alive != 1] <- NA
  fallback <- first_valid(bwd(paste0(prefix, "age")), invalid = numeric(0))
  coalesce(primary, fallback)
}
df$fage <- age_with_gate("f", df$falive)
df$mage <- age_with_gate("m", df$malive)

for (P in c("m", "f")) {
  for (i in 1:11) {
    v <- df[[paste0(P, "alive", i)]]
    v <- ifelse(v == 2, 0, ifelse(v < 0, NA, v))
    df[[paste0(P, "lve", i)]] <- v
  }
}
df$nmalive_d <- as.integer(rowSums(df[paste0("mlve", 1:11)], na.rm = TRUE) > 0)
df$nfalive_d <- as.integer(rowSums(df[paste0("flve", 1:11)], na.rm = TRUE) > 0)

# Education and covariates
df$edend <- first_valid(fwd("edend"))

df$sex      <- first_valid(fwd("sex"), invalid = numeric(0))
df$difjobm  <- first_valid(fwd("difjobm"), invalid = c(-1, -8, -9))
df$nonwhite <- first_valid(fwd("nonwhite"), invalid = c(-2, -8, -9))
df$nsibs    <- first_valid(fwd("nsibs"), invalid = c(-1, -2, -8, -9))

df$par_class2 <- case_when(
  df$difjobm %in% c(2, 3, 4)            ~ 1,
  df$difjobm %in% c(1, 5, 7, 8)         ~ 2,
  df$difjobm %in% c(6, 9, 10, 11, 12)   ~ 3,
  df$difjobm %in% c(13, 14, 15, 16)     ~ 4,
  TRUE                                   ~ NA_real_
)

df <- df %>% rename(bmonth = dobmonth, byear = dobyear)

final <- df %>%
  select(idauniq, bmonth, byear, fagedied, magedied, fage, mage, falive, malive,
         nmalive_d, nfalive_d, edend, sex, par_class2, ngor, nonwhite, nsibs,
         all_of(paste0("wgt", 1:11))) %>%
  arrange(idauniq)

saveRDS(final, file.path(tempdata, "ELSA_base_dataset.rds"))
cat("Saved", nrow(final), "respondents to ELSA_base_dataset.rds\n")
