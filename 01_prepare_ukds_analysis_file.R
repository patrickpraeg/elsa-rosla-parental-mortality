suppressPackageStartupMessages(library(tidyverse))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
tempdata   <- file.path(dirname(script_dir), "data/tempdata")

df <- readRDS(file.path(tempdata, "ELSA_base_dataset.rds"))

df <- df %>%
  mutate(
    ngor       = as.integer(ngor),
    par_class2 = as.integer(par_class2),
    falive     = as.integer(falive),
    malive     = as.integer(malive)
  )

# Birth date and reform variables
invalid_to_na <- function(x) {
  x <- as.numeric(x)
  x[x %in% c(-1, -2, -8, -9)] <- NA_real_
  x
}
df$byear <- invalid_to_na(df$byear)
df <- df %>%
  mutate(byear  = if_else(byear %in% c(-7, 1908, 1909, 1911), 1912, byear),
         bmonth = if_else(!is.na(bmonth) & bmonth <= 0, NA_real_, as.numeric(bmonth)))

valid_dob <- !is.na(df$byear) & !is.na(df$bmonth) & df$byear > 0 & df$bmonth > 0
mob_str <- rep(NA_character_, nrow(df))
mob_str[valid_dob] <- sprintf("%d-%02d-01", df$byear[valid_dob], df$bmonth[valid_dob])
df$MoB <- as.Date(mob_str)

df <- df %>%
  mutate(
    post47     = if_else(is.na(MoB), NA_integer_,
                         as.integer(MoB >= as.Date("1933-04-01"))),
    R_months47 = ((byear - 1933) * 12) + (bmonth - 4),
    post72     = if_else(is.na(MoB), NA_integer_,
                         as.integer(MoB >= as.Date("1957-09-01"))),
    R_months72 = ((byear - 1957) * 12) + (bmonth - 9)
  )

saveRDS(df, file.path(tempdata, "ELSA_sn8346_bmonth-year.rds"))
cat("Saved", nrow(df), "respondents to ELSA_sn8346_bmonth-year.rds\n")
