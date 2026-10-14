suppressPackageStartupMessages({library(tidyverse); library(haven)})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir  <- get_script_dir()
data_path   <- file.path(dirname(script_dir), "data/hse")
output_path <- file.path(dirname(script_dir), "figures")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Health Survey for England waves
read_hse <- function(file, year, vars = c("age", "educend", "sex")) {
  df <- zap_labels(read_dta(file.path(data_path, file)))
  names(df) <- tolower(names(df))
  df <- df[, intersect(vars, names(df)), drop = FALSE]
  df$year <- year
  df
}

waves <- list(
  list("hse2013ai.dta",            2013),
  list("hse2012ai.dta",            2012),
  list("hse2011ai.dta",            2011),
  list("hse10ai.dta",              2010),
  list("hse09ai.dta",              2009),
  list("hse08ai.dta",              2008),
  list("hse07ai.dta",              2007),
  list("hse06ai.dta",              2006),
  list("hse05ai.dta",              2005),
  list("hse04etha.dta",            2004),
  list("hse03ai.dta",              2003),
  list("hse02ai.dta",              2002),
  list("hse01ai.dta",              2001),
  list("hse00ai.dta",              2000),
  list("hse99et3.dta",             1999),
  list("hse98ai.dta",              1998),
  list("hse97ai_v5_seg.dta",       1997),
  list("hse96ai_v5_with_seg.dta",  1996),
  list("hse95ai_v6_seg.dta",       1995),
  list("hse94ai_v4_seg.dta",       1994),
  list("3316.dta",                 1993)
)

df <- map_dfr(waves, ~ read_hse(.x[[1]], .x[[2]]))

df_9192 <- read_hse("3238.dta", NA, vars = c("age", "educend", "sex", "year"))
df_9192 <- df_9192 %>%
  mutate(year = case_when(
    year == 91 ~ 1991,
    year == 92 ~ 1992,
    TRUE       ~ as.numeric(year)
  ))

df <- bind_rows(df, df_9192)

# Sample restrictions and school-leaving age categories
df <- df %>% filter(age != -1, !is.na(age))
df <- df %>% mutate(yob = year - age)

invalid_educend <- c(-9, -8, -1, 1, 2, 9, NA)
df <- df %>% filter(!educend %in% invalid_educend)

df <- df %>%
  mutate(educendx = case_when(
    educend == 3       ~ "14 y. or under",
    educend == 4       ~ "15 y.",
    educend == 5       ~ "16 y.",
    educend %in% 6:8   ~ "17 y. or over",
    TRUE               ~ NA_character_
  )) %>%
  filter(!is.na(educendx)) %>%
  mutate(educendx = factor(educendx,
                           levels = c("14 y. or under", "15 y.",
                                      "16 y.",           "17 y. or over")))

collapsed <- df %>%
  filter(!is.na(yob)) %>%
  group_by(yob, educendx) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(yob) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  filter(yob >= 1920, yob <= 1980)

# Figure 1
line_colours <- c(
  "17 y. or over"  = "#1A476F",
  "16 y."          = "#5293FF",
  "15 y."          = "#FF58AC",
  "14 y. or under" = "#9DB0C1"
)
line_types <- c(
  "17 y. or over"  = "dotted",
  "16 y."          = "solid",
  "15 y."          = "solid",
  "14 y. or under" = "dashed"
)
line_widths <- c(
  "17 y. or over"  = 0.6,
  "16 y."          = 1.0,
  "15 y."          = 1.0,
  "14 y. or under" = 0.6
)

fig1 <- ggplot(collapsed, aes(x = yob, y = pct, colour = educendx,
                              linetype = educendx, linewidth = educendx)) +
  annotate("rect", xmin = 1956.4, xmax = 1957.6,
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) +
  annotate("rect", xmin = 1932.4, xmax = 1933.6,
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) +
  geom_line() +
  annotate("text", x = 1938, y = 90, label = "1947 reform",  size = 3.2) +
  annotate("text", x = 1962, y = 90, label = "1972 reform",  size = 3.2) +
  scale_colour_manual(values = line_colours,   name = "School-leaving age") +
  scale_linetype_manual(values = line_types,   name = "School-leaving age") +
  scale_linewidth_manual(values = line_widths, name = "School-leaving age") +
  scale_x_continuous(breaks = seq(1920, 1980, by = 10),
                     limits = c(1920, 1980), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     limits = c(0, 100), expand = expansion(mult = 0.02)) +
  labs(
    x = "Year of birth",
    y = "% of birth cohort by age of leaving education"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#E9E9E9"),
    panel.grid.major.x = element_blank(),
    axis.line.x        = element_line(colour = "#4C4C4C", linewidth = 0.4),
    axis.line.y        = element_line(colour = "#4C4C4C", linewidth = 0.4),
    axis.ticks         = element_blank(),
    legend.position           = "inside",
    legend.position.inside    = c(0.985, 0.985),
    legend.justification.inside = c(1, 1),
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 9),
    legend.key.width   = unit(0.55, "cm"),
    legend.key.height  = unit(0.40, "cm"),
    legend.margin      = margin(3, 5, 3, 5),
    legend.background  = element_rect(fill = "white", colour = "#5E5E5E",
                                      linewidth = 0.3)
  ) +
  guides(colour    = guide_legend(nrow = 4, reverse = TRUE),
         linetype  = guide_legend(nrow = 4, reverse = TRUE),
         linewidth = guide_legend(nrow = 4, reverse = TRUE))

ggsave(file.path(output_path, "figure1.pdf"), fig1,
       width = 6, height = 4, device = cairo_pdf)
