suppressPackageStartupMessages({library(tidyverse); library(haven)})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
pkg_root   <- dirname(script_dir)

data_natcen <- file.path(pkg_root, "data/elsa from natcen")
tempdata    <- file.path(pkg_root, "data/tempdata")
dir.create(tempdata, showWarnings = FALSE, recursive = TRUE)

# Merge waves 1-8
waves <- lapply(1:8, function(i) {
  df <- read_sav(file.path(data_natcen, paste0("Wave", i, "_final.sav")))
  names(df) <- tolower(names(df))
  df$idauniq <- zap_labels(df$idauniq)
  df$wave <- i
  names(df) <- ifelse(names(df) == "idauniq", "idauniq",
                      paste0(names(df), "_w", i))
  df
})
elsa <- Reduce(function(a, b) full_join(a, b, by = "idauniq"), waves)

coalesce_invalid <- function(..., invalid = c(-1, -2, -8, -9)) {
  vals <- list(...)
  result <- rep(NA_real_, length(vals[[1]]))
  for (v in vals) {
    v_clean <- if_else(v %in% invalid, NA_real_, as.numeric(v))
    result <- coalesce(result, v_clean)
  }
  result
}

# Parents' survival
for (i in 1:8) {
  elsa <- elsa %>% rename(
    !!paste0("mage_w", i) := !!paste0("mothage_w", i),
    !!paste0("fage_w", i) := !!paste0("fathage_w", i)
  )
}

elsa <- elsa %>%
  mutate(
    fagedied = coalesce_invalid(!!!syms(paste0("fagedied_w", 1:8))),
    magedied = coalesce_invalid(!!!syms(paste0("magedied_w", 1:8))),
    falive   = coalesce_invalid(!!!syms(paste0("falive_w", 8:1))),
    malive   = coalesce_invalid(!!!syms(paste0("malive_w", 8:1))),
    fage     = coalesce_invalid(!!!syms(paste0("fage_w", 8:1))),
    mage     = coalesce_invalid(!!!syms(paste0("mage_w", 8:1)))
  )

# Birth date and reform variables
elsa <- elsa %>%
  rename(
    dob_month_w1 = indobm_w1,
    dob_month_w2 = indobm_w2
  ) %>%
  mutate(
    across(c(dob_month_w7, dob_month_w8), as.numeric)
  )

elsa <- elsa %>%
  mutate(bmonth = coalesce_invalid(!!!syms(paste0("dob_month_w", 1:8)),
                                   invalid = -1))

elsa <- elsa %>%
  rename(
    dob_year_w1 = dhdobyr_w1,
    dob_year_w2 = dhdobyr_w2,
    dob_year_w3 = dhdobyr_w3,
    dob_year_w4 = indobyr_w4,
    dob_year_w5 = indobyr_w5,
    dob_year_w6 = indobyr_w6,
    dob_year_w7 = indobyr_w7,
    dob_year_w8 = indobyr_w8
  )

elsa <- elsa %>%
  mutate(
    byear = coalesce_invalid(!!!syms(paste0("dob_year_w", 1:8)),
                             invalid = c(-8, -7))
  ) %>%
  mutate(byear = if_else(byear %in% c(-7, 1908, 1909, 1911), 1912, byear))

valid_dob <- !is.na(elsa$byear) & !is.na(elsa$bmonth)
mob_str <- rep(NA_character_, nrow(elsa))
mob_str[valid_dob] <- sprintf("%d-%02d-01", elsa$byear[valid_dob], elsa$bmonth[valid_dob])
elsa$MoB <- as.Date(mob_str)

elsa <- elsa %>%
  mutate(
    post47     = as.integer(MoB >= as.Date("1933-04-01")),
    R_months47 = ((byear - 1933) * 12) + (bmonth - 4),
    post72     = as.integer(MoB >= as.Date("1957-09-01")),
    R_months72 = ((byear - 1957) * 12) + (bmonth - 9)
  )

# Education and covariates
elsa <- elsa %>%
  mutate(
    edend    = coalesce_invalid(!!!syms(paste0("edend_w", 1:8))),
    sex      = zap_labels(coalesce(!!!syms(paste0("sex_w", 1:8)))),
    difjobm  = coalesce_invalid(!!!syms(paste0("difjobm_w", 1:8))),
    nonwhite = coalesce_invalid(!!!syms(paste0("nonwhite_w", 1:8)),
                                invalid = c(-2, -8, -9)),
    nsibs    = coalesce_invalid(!!!syms(paste0("nsibs_w", 1:8)))
  )

elsa <- elsa %>%
  mutate(
    par_class2 = case_when(
      difjobm %in% c(2, 3, 4)          ~ 1L,
      difjobm %in% c(1, 5, 7, 8)       ~ 2L,
      difjobm %in% c(6, 9, 10, 11, 12) ~ 3L,
      difjobm %in% c(13, 14, 15, 16)   ~ 4L,
      TRUE                             ~ NA_integer_
    )
  )

# Region
elsa <- elsa %>% mutate(gor = NA_character_)

for (i in 1:8) {
  col <- paste0("gor_w", i)
  if (col %in% names(elsa)) {
    elsa <- elsa %>%
      mutate(gor = if_else(is.na(gor) | gor == "W99999999",
                           as.character(.data[[col]]), gor))
  }
}

elsa <- elsa %>%
  mutate(gor = case_when(
    gor == "E12000001" ~ "A",
    gor == "E12000002" ~ "B",
    gor == "E12000003" ~ "D",
    gor == "E12000004" ~ "E",
    gor == "E12000005" ~ "F",
    gor == "E12000006" ~ "G",
    gor == "E12000007" ~ "H",
    gor == "E12000008" ~ "J",
    gor == "E12000009" ~ "K",
    gor == "W99999999" ~ NA_character_,
    TRUE               ~ gor
  )) %>%
  mutate(ngor = factor(gor,
                       levels = c("A","B","D","E","F","G","H","J","K"),
                       labels = c("North East","North West",
                                  "Yorkshire and The Humber",
                                  "East Midlands","West Midlands",
                                  "East of England","London",
                                  "South East","South West")))

# Parents observed alive
for (i in 1:8) {
  for (P in c("m", "f")) {
    col_alive <- paste0(P, "alive_w", i)
    col_lve   <- paste0(P, "lve", i)
    elsa <- elsa %>%
      mutate(!!col_lve := case_when(
        .data[[col_alive]] == 1  ~ 1,
        .data[[col_alive]] == 2  ~ 0,
        .data[[col_alive]] <  0  ~ NA_real_,
        TRUE                     ~ NA_real_
      ))
  }
}

elsa <- elsa %>%
  mutate(
    nmalive_d = as.integer(rowSums(select(., paste0("mlve", 1:8)), na.rm = TRUE) > 0),
    nfalive_d = as.integer(rowSums(select(., paste0("flve", 1:8)), na.rm = TRUE) > 0)
  )

elsa_reduced <- elsa %>%
  select(idauniq, R_months72, R_months47, post72, post47, MoB, bmonth, byear,
         fagedied, magedied, fage, mage, falive, malive, nmalive_d, nfalive_d,
         edend, sex, par_class2, ngor, nonwhite, nsibs)

saveRDS(elsa_reduced, file.path(tempdata, "ELSA_bmonth-year.rds"))
cat("Saved", nrow(elsa_reduced), "respondents to ELSA_bmonth-year.rds\n")
