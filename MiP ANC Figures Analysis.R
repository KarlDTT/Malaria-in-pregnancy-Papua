# ============================================================================
# Malaria in Pregnancy Analysis - Papua, Indonesia
# Target Journal: Lancet Regional Health - Western Pacific
# Figure Analysis
# ============================================================================


library(pacman)
pacman::p_load(tidyverse, sf, viridis, patchwork, scales, lubridate, 
               ggrepel, RColorBrewer, rio, corrr, GGally, stats)

province_order <- c("Southwest Papua", "West Papua", "Central Papua", "Papua",
                    "Highland Papua", "South Papua" )


# ---- Data ------------------------------------------------------------------
Yearly_Data_Complete  <- import("data/Yearly Data Complete.xlsx",  sheet = "Sheet1")
Monthly_Data_Complete <- import("data/Monthly Data Complete.xlsx", sheet = "Sheet1")
Regmal1_MiP           <- import("data/Regmal1 MiP.xlsx", sheet = "Sheet1")

Papua_map    <- st_read("data/pulau_papua/papua_admin2.shp")
Papua_border <- st_read("data/pulau_papua/papua_admin1.shp")

province_order <- c("West Papua", "Southwest Papua", "Central Papua", "Papua",
                    "Highland Papua", "South Papua")


# Join yearly data with district shapefile (Fig 1C, Fig S2B use this pattern)
Papua_map_TPR <- Papua_map %>%
  left_join(
    Yearly_Data_Complete %>%
      filter(Malaria_Testing_ANC >= 30) %>%
      group_by(Regency) %>%
      summarise(
        `Total MiPs`            = sum(MiPs, na.rm = TRUE),
        `Total Malaria Testing` = sum(Malaria_Testing_ANC, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        `Test Positivity Rate` = (`Total MiPs` / `Total Malaria Testing`) * 100,
        admin2 = toupper(Regency)
      ),
    by = c("admin2" = "admin2")
  ) %>%
  ungroup()


# Shared theme components (avoids repeating the same theme() block everywhere)
base_theme <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.4, colour = "black"),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )


# ============================================================================
#Figure 1A
# ============================================================================


# 1A: Mean ANC1 / testing coverage by province (bar + range)
Fig1A <- Yearly_Data_Complete %>%
  filter(Malaria_Testing_ANC >= 30) %>%
  group_by(Province, Regency) %>%
  mutate(
    `ANC1 Coverage`    = (ANC1_imputed / Pregnancies_imputed) * 100,
    `Testing Coverage` = (Malaria_Testing_ANC / Pregnancies_imputed) * 100,
    Province = factor(Province, levels = province_order)
  ) %>%
  summarise(
    `Mean ANC1 Coverage`    = mean(`ANC1 Coverage`, na.rm = TRUE),
    `Mean Testing Coverage` = mean(`Testing Coverage`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Province) %>%
  summarise(
    Mean_ANC1_Coverage    = mean(`Mean ANC1 Coverage`, na.rm = TRUE),
    Min_ANC1_Coverage     = min(`Mean ANC1 Coverage`, na.rm = TRUE),
    Max_ANC1_Coverage     = max(`Mean ANC1 Coverage`, na.rm = TRUE),
    Mean_Testing_Coverage = mean(`Mean Testing Coverage`, na.rm = TRUE),
    Min_Testing_Coverage  = min(`Mean Testing Coverage`, na.rm = TRUE),
    Max_Testing_coverage  = max(`Mean Testing Coverage`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Mean_ANC1_Coverage, Mean_Testing_Coverage),
    names_to = "Metric", values_to = "Coverage"
  ) %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("Mean_ANC1_Coverage", "Mean_Testing_Coverage"),
                    labels = c("ANC1 Coverage", "Testing Coverage")),
    ymin = ifelse(Metric == "Testing Coverage", Min_Testing_Coverage, Min_ANC1_Coverage),
    ymax = ifelse(Metric == "Testing Coverage", Max_Testing_coverage, Max_ANC1_Coverage)
  ) %>%
  ggplot(aes(x = Province, y = Coverage, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.8), width = 0.3) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02"), name = "") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(tag = "A", x = "", y = "Coverage (%)") +
  base_theme +
  theme(
    legend.position          = "inside",
    legend.position.inside   = c(0.7, 0.975),
    legend.background        = element_rect(fill = "transparent", colour = NA),
    legend.key               = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(fill = guide_legend(nrow = 1))


# 1B: ANC1 test positivity over time, by province
Fig1B <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  summarise(`Total ANC1 Testing` = sum(Malaria_Testing_ANC, na.rm = TRUE), .groups = "drop") %>%
  filter(`Total ANC1 Testing` >= 30) %>%
  left_join(Monthly_Data_Complete, by = c("Province", "Regency", "Year")) %>%
  group_by(Province, Month, Year) %>%
  summarise(
    `Total ANC1 Testing` = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total MiPs`         = sum(MiPs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate` = replace_na((`Total MiPs` / `Total ANC1 Testing`) * 100, 0),
    date     = ymd(paste0(Year, "-", Month, "-01")),
    Province = factor(Province, levels = province_order)
  ) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = `Test Positivity Rate`, group = Province, col = Province),
            linewidth = 0.8, na.rm = TRUE) +
  scale_x_date(
    breaks = seq(ymd("2019-07-01"), ymd("2023-07-01"), by = "1 year"),
    date_labels = "%Y", expand = c(0.01, 0.01), minor_breaks = NULL
  ) +
  geom_vline(
    xintercept = seq(ymd("2019-01-01"), ymd("2023-12-31"), by = "1 year"),
    linetype = "dashed", colour = "grey60", linewidth = 0.3
  ) +
  labs(tag = "B", x = "", y = "ANC1 test positivity (%)", colour = NULL) +
  base_theme +
  theme(
    legend.position         = "inside",
    legend.position.inside  = c(0.7, 0.9),
    legend.background       = element_rect(fill = "transparent", colour = NA),
    legend.key              = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(color = guide_legend(nrow = 3))


# 1C: District-level test positivity map
Fig1C <- ggplot(Papua_map_TPR) +
  geom_sf(aes(fill = `Test Positivity Rate`), colour = "white", size = 0.2) +
  geom_sf(data = Papua_border, fill = NA, colour = "black", linewidth = 0.8) +
  scale_fill_viridis(
    option = "plasma", name = "ANC test\npositivity (%)",
    breaks = seq(0, 100, 20), labels = seq(0, 100, 20),
    limits = c(0, 40), oob = scales::squish,
    guide = guide_colourbar(barwidth = unit(2, "cm"), barheight = unit(0.4, "cm"))
  ) +
  labs(tag = "C") +
  theme_void(base_size = 10) +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.3, 0.1),
    legend.direction       = "horizontal"
  )


#Making Figure 1
Fig1 <- free(Fig1A) / free(Fig1B) / Fig1C +
  plot_layout(heights = c(0.8, 0.8, 1))
ggsave("Fig1.jpg", Fig1, width = 25, height = 30, dpi = 600, units = "cm")


# ============================================================================
# Supplementary Figure 1: heatmaps of coverage / positivity by year
# ============================================================================


FigS1A <- Yearly_Data_Complete %>%
  mutate(
    `Testing Ratio` = (Malaria_Testing_ANC / Pregnancies_imputed) * 100,
    Province = factor(Province, levels = province_order),
    Year = factor(Year)
  ) %>%
  ggplot(aes(x = Year, y = Regency, fill = `Testing Ratio`)) +
  geom_tile(width = 0.9) +
  scale_fill_viridis(
    option = "magma", name = "Testing\ncoverage (%)",
    limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = c("0", "25", "50", "75", ">100"), oob = scales::squish
  ) +
  facet_grid(Province ~ ., scales = "free_y", space = "free_y") +
  labs(tag = "A", x = "", y = "District") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  base_theme +
  theme(strip.text.y = element_blank())



FigS1B <- Yearly_Data_Complete %>%
  mutate(
    MiPs = ifelse(is.na(MiPs), 0, MiPs),
    `Test Positivity Rate` = MiPs / Malaria_Testing_ANC * 100,
    Province = factor(Province, levels = province_order),
    Year = factor(Year)
  ) %>%
  ggplot(aes(x = Year, y = Regency, fill = `Test Positivity Rate`)) +
  geom_tile(width = 0.9) +
  scale_fill_viridis(
    option = "magma", name = "ANC1 test\npositivity (%)",
    limits = c(0, 100), oob = scales::squish
  ) +
  facet_grid(Province ~ ., scales = "free_y", space = "free_y") +
  labs(tag = "B", x = "") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  base_theme +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())


#Making Supplementary Figure 1
FigS1 <- FigS1A | FigS1B
ggsave("FigS1.png", FigS1, width = 20, height = 25, dpi = 600, units = "cm")


# ============================================================================
# Figure 2: species / maternal age trends
# ============================================================================


#Figure 2A
Fig2A <- Regmal1_MiP %>%
  group_by(Province, Month, Year) %>%
  summarise(
    P.falciparum = sum(Parasite == "P.falciparum", na.rm = TRUE),
    P.vivax      = sum(Parasite == "P.vivax", na.rm = TRUE),
    P.mix        = sum(Parasite == "P.Mix", na.rm = TRUE),
    P.malariae   = sum(Parasite == "P.malariae", na.rm = TRUE),
    P.ovalae     = sum(Parasite == "P.ovalae", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total     = P.falciparum + P.vivax + P.mix + P.malariae + P.ovalae,
    `Vivax %` = replace_na((P.vivax / Total) * 100, 0),
    Date      = ymd(paste0(Year, "-", Month, "-01")),
    Province  = factor(Province, levels = province_order)
  ) %>%
  arrange(Province, Date) %>%
  group_by(Province) %>%
  mutate(prop_vivax_smooth = zoo::rollmean(`Vivax %`, k = 3, fill = NA, align = "center")) %>%
  ungroup() %>%
  ggplot(aes(x = Date, y = prop_vivax_smooth, colour = Province)) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  scale_x_date(
    breaks = seq(ymd("2019-07-01"), ymd("2023-07-01"), by = "1 year"),
    date_labels = "%Y", expand = c(0.01, 0.01), minor_breaks = NULL
  ) +
  geom_vline(
    xintercept = seq(ymd("2019-01-01"), ymd("2023-12-31"), by = "1 year"),
    linetype = "dashed", colour = "grey60", linewidth = 0.3
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(tag = "A", x = NULL, y = "P. vivax (%)", colour = NULL) +
  base_theme +
  theme(
    legend.position         = "inside",
    legend.position.inside  = c(0.7, 0.9),
    legend.background       = element_rect(fill = "transparent", colour = NA),
    legend.key              = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(color = guide_legend(nrow = 2))


#Figure 2B
Fig2B <- Regmal1_MiP %>%
  mutate(
    Age_Group = case_when(
      Age_Group %in% c("10-14", "15-19")         ~ "10-19",
      Age_Group %in% c("20-24", "25-29", "30-34") ~ "20-34",
      Age_Group == "35+"                          ~ "35+",
      TRUE ~ Age_Group
    ),
    Age_Group = factor(Age_Group, levels = c("10-19", "20-34", "35+")),
    Date      = ymd(paste0(Year, "-", Month, "-01")),
    Quarter   = quarter(Date),
    year_quarter = paste0(Year, "-Q", Quarter)
  ) %>%
  group_by(Year, Quarter, year_quarter) %>%
  summarise(
    `10-19` = sum(Age_Group == "10-19", na.rm = TRUE),
    `20-34` = sum(Age_Group == "20-34", na.rm = TRUE),
    `35+`   = sum(Age_Group == "35+", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total = `10-19` + `20-34` + `35+`,
    `Adolescent Pregnancy %` = (`10-19` / Total) * 100,
    `Optimal Pregnancy %`    = (`20-34` / Total) * 100,
    `Advanced Pregnancy %`   = (`35+` / Total) * 100,
    Date = as.Date(paste0(Year, "-", Quarter * 3 - 2, "-01"))
  ) %>%
  pivot_longer(
    cols = `Adolescent Pregnancy %`:`Advanced Pregnancy %`,
    names_to = "Age Group", values_to = "Proportion"
  ) %>%
  mutate(
    `Age Group` = factor(`Age Group`,
                         levels = c("Adolescent Pregnancy %", "Optimal Pregnancy %", "Advanced Pregnancy %"),
                         labels = c("Adolescent (10-19y)", "Optimal (20-34y)", "Advanced (\u226535y)"))
  ) %>%
  ggplot(aes(x = Date, y = Proportion, colour = `Age Group`)) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Maternal age") +
  scale_x_date(
    breaks = seq(ymd("2019-07-01"), ymd("2023-07-01"), by = "1 year"),
    date_labels = "%Y", expand = c(0.01, 0.01), minor_breaks = NULL
  ) +
  geom_vline(
    xintercept = seq(ymd("2019-01-01"), ymd("2023-12-31"), by = "1 year"),
    linetype = "dashed", colour = "grey60", linewidth = 0.3
  ) +
  labs(tag = "B", x = "", y = "Proportion of cases (%)") +
  base_theme +
  theme(legend.margin = margin(t = -2), legend.box.margin = margin(t = -4))


#Figure 2C
Fig2C <- Regmal1_MiP %>%
  mutate(
    Age_Group = case_when(
      Age_Group %in% c("10-14", "15-19")         ~ "10-19",
      Age_Group %in% c("20-24", "25-29", "30-34") ~ "20-34",
      Age_Group == "35+"                          ~ "35+",
      TRUE ~ Age_Group
    ),
    Age_Group = factor(Age_Group, levels = c("10-19", "20-34", "35+"))
  ) %>%
  group_by(Province) %>%
  summarise(
    `10-19` = sum(Age_Group == "10-19", na.rm = TRUE),
    `20-34` = sum(Age_Group == "20-34", na.rm = TRUE),
    `35+`   = sum(Age_Group == "35+", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total = `10-19` + `20-34` + `35+`,
    `Adolescent Pregnancy %` = (`10-19` / Total) * 100,
    `Optimal Pregnancy %`    = (`20-34` / Total) * 100,
    `Advanced Pregnancy %`   = (`35+` / Total) * 100,
    Province = factor(Province, levels = province_order)
  ) %>%
  pivot_longer(
    cols = `Adolescent Pregnancy %`:`Advanced Pregnancy %`,
    names_to = "Age Group", values_to = "Proportion"
  ) %>%
  mutate(
    `Age Group` = factor(`Age Group`,
                         levels = c("Adolescent Pregnancy %", "Optimal Pregnancy %", "Advanced Pregnancy %"),
                         labels = c("Adolescent (10-19y)", "Optimal (20-34y)", "Advanced (\u226535y)"))
  ) %>%
  ggplot(aes(x = Province, y = Proportion, fill = `Age Group`)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set1", name = "Maternal age") +
  labs(tag = "C", x = "", y = "Proportion of cases (%)") +
  base_theme +
  theme(legend.margin = margin(t = -2), legend.box.margin = margin(t = -4))


#Making Figure 2
Fig2 <- Fig2A / Fig2B / Fig2C
ggsave("Fig2.jpg", Fig2, width = 25, height = 30, dpi = 600, units = "cm")


# ============================================================================
# Supplementary Figure 2: species composition
# ============================================================================


#Figure S2A
FigS2A <- Regmal1_MiP %>%
  mutate(Date = ymd(paste0(Year, "-", Month, "-01"))) %>%
  group_by(Date) %>%
  summarise(
    P.falciparum      = sum(Parasite == "P.falciparum", na.rm = TRUE),
    P.vivax           = sum(Parasite == "P.vivax", na.rm = TRUE),
    `Mix (Pf + Pv)`   = sum(Parasite == "P.Mix", na.rm = TRUE),
    P.malariae        = sum(Parasite == "P.malariae", na.rm = TRUE),
    P.ovale           = sum(Parasite == "P.ovale", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = P.falciparum:P.ovale, names_to = "Species", values_to = "Count") %>%
  mutate(Species = factor(Species,
                          levels = c("P.falciparum", "P.vivax", "Mix (Pf + Pv)", "P.malariae", "P.ovale"))) %>%
  ggplot(aes(x = Date, y = Count, fill = Species)) +
  geom_area(alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  scale_x_date(
    breaks = seq(ymd("2019-07-01"), ymd("2023-07-01"), by = "1 year"),
    date_labels = "%Y", expand = c(0.01, 0.01), minor_breaks = NULL
  ) +
  geom_vline(
    xintercept = seq(ymd("2019-01-01"), ymd("2023-12-31"), by = "1 year"),
    linetype = "dashed", colour = "grey60", linewidth = 0.3
  ) +
  labs(tag = "A", x = "", y = "Number of positive cases", fill = NULL) +
  base_theme


#Figure S2B
Papua_map_PV <- Papua_map %>%
  left_join(
    Regmal1_MiP %>%
      group_by(Regency) %>%
      summarise(
        P.falciparum = sum(Parasite == "P.falciparum", na.rm = TRUE),
        P.vivax      = sum(Parasite == "P.vivax", na.rm = TRUE),
        P.mix        = sum(Parasite == "P.Mix", na.rm = TRUE),
        P.malariae   = sum(Parasite == "P.malariae", na.rm = TRUE),
        P.ovalae     = sum(Parasite == "P.ovalae", na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Total     = P.falciparum + P.vivax + P.mix + P.malariae + P.ovalae,
        `Vivax %` = (P.vivax / Total) * 100,
        admin2    = toupper(Regency)
      ),
    by = c("admin2" = "admin2")
  ) %>%
  ungroup()


FigS2B <- ggplot(Papua_map_PV) +
  geom_sf(aes(fill = `Vivax %`), colour = "white", size = 0.2) +
  geom_sf(data = Papua_border, fill = NA, colour = "black", linewidth = 0.8) +
  scale_fill_viridis(option = "plasma", name = "P. vivax (%)") +
  labs(tag = "B") +
  theme_void(base_size = 10) +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.3, 0.1),
    legend.direction       = "horizontal"
  )


#Making Supplementary Figure 2
FigS2 <- free(FigS2A) / FigS2B
ggsave("FigS2.jpg", FigS2, width = 25, height = 30, dpi = 600, units = "cm")


# ============================================================================
# Figure 3: coverage vs. ANC1 volume, quarterly heatmap, district variation
# ============================================================================


#Figure 3A
Fig3A <- Yearly_Data_Complete %>%
  filter(Malaria_Testing_ANC >= 30) %>%
  group_by(Province) %>%
  summarise(
    `Total Pregnancies` = sum(Pregnancies_imputed, na.rm = TRUE),
    `Total ANC1`        = sum(ANC1_imputed, na.rm = TRUE),
    `Total Testing`     = sum(Malaria_Testing_ANC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `ANC1 Coverage`    = pmin((`Total ANC1` / `Total Pregnancies`) * 100, 100),
    `Testing Coverage` = (`Total Testing` / `Total Pregnancies`) * 100,
    `ANC1 Group` = cut(`Total ANC1`, breaks = c(0, 40000, 80000, Inf),
                       labels = c("<40,000", "40,000-80,000", ">80,000"), right = FALSE)
  ) %>%
  ggplot(aes(x = `ANC1 Coverage`, y = `Testing Coverage`, size = `ANC1 Group`, label = Province)) +
  geom_point(alpha = 0.7, color = "#e41a1c") +
  geom_text_repel(size = 3, max.overlaps = 10) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 97, linetype = "dashed", color = "gray50") +
  scale_size_manual(
    name = "Total\nANC1 visits",
    values = c("<40,000" = 4, "40,000-80,000" = 8, ">80,000" = 12)
  ) +
  labs(tag = "A", x = "ANC1 coverage (%)", y = "Malaria testing coverage (%)") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  base_theme +
  theme(legend.direction = "horizontal")


#Figure 3B 
Fig3B <- Monthly_Data_Complete %>%
  mutate(Quarter = quarter(Date), year_quarter = paste0(Year, "-Q", Quarter)) %>%
  group_by(Province, Year, Quarter, year_quarter) %>%
  summarise(
    `Total Pregnancies` = sum(Pregnancies_imputed),
    `Total Testing`     = sum(Malaria_Testing_ANC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Testing Coverage` = (`Total Testing` / `Total Pregnancies`) * 100,
    Province = factor(Province, levels = rev(province_order))
  ) %>%
  arrange(Year, Quarter) %>%
  mutate(year_quarter = factor(year_quarter, levels = unique(year_quarter))) %>%
  ggplot(aes(x = year_quarter, y = Province, fill = `Testing Coverage`)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis(option = "magma", name = "Testing\ncoverage (%)", limits = c(0, 10)) +
  scale_x_discrete(labels = rep(c("1", "2", "3", "4"), length(unique(Monthly_Data_Complete$Year)))) +
  labs(tag = "B", x = "Quarter", y = "") +
  facet_grid(. ~ Year, scales = "free_x", space = "free_x", switch = "x") +
  base_theme +
  theme(
    panel.spacing.x  = unit(2, "pt"),
    legend.margin    = margin(t = -2),
    legend.box.margin = margin(t = -4)
  )


#Figure 3C
Fig3C <- Yearly_Data_Complete %>%
  filter(Malaria_Testing_ANC >= 30) %>%
  group_by(Province, Regency) %>%
  summarise(
    `Total Pregnancies` = sum(Pregnancies_imputed),
    `Total Testing`     = sum(Malaria_Testing_ANC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Testing Coverage` = (`Total Testing` / `Total Pregnancies`) * 100,
    Province = factor(Province, levels = province_order)
  ) %>%
  ggplot(aes(x = Province, y = `Testing Coverage`, fill = Province)) +
  geom_boxplot(show.legend = FALSE, alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(tag = "C", x = NULL, y = "Testing coverage (%)") +
  base_theme


#Making Figure 3
Fig3 <- (Fig3A | Fig3B) / Fig3C
ggsave("Fig3.jpg", Fig3, width = 19, height = 20, dpi = 600, units = "cm")


# ============================================================================
# Supplementary Figure 3: district-level correlation (rho + 95% CI)
# ============================================================================


#Figure S3
FigS3 <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  mutate(
    MiPs = ifelse(is.na(MiPs), 0, MiPs),
    `Test Positivity Rate` = MiPs / Malaria_Testing_ANC * 100
  ) %>%
  ungroup() %>%
  mutate(across(c(`Test Positivity Rate`, `Population Test Positivity Rate`), ~replace_na(., 0))) %>%
  group_by(Province, Regency) %>%
  summarise(
    rho     = cor(`Test Positivity Rate`, `Population Test Positivity Rate`,
                  method = "spearman", use = "complete.obs"),
    n       = n(),
    sd_test = sd(`Test Positivity Rate`, na.rm = TRUE),
    sd_pop  = sd(`Population Test Positivity Rate`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    no_variation = sd_test == 0 | sd_pop == 0,
    se        = sqrt((1 - rho^2) / (n - 2)),
    ci_lower  = tanh(atanh(rho) - 1.96 * se),
    ci_upper  = tanh(atanh(rho) + 1.96 * se),
    rho_plot  = ifelse(no_variation, 0, rho)   # districts with no variance plotted at 0, flagged
  ) %>%
  arrange(Province, rho_plot) %>%
  ggplot(aes(x = rho_plot, y = reorder(Regency, rho_plot), colour = Province)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(data = ~ filter(.x, !no_variation),
                 aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, alpha = 0.6) +
  geom_point(aes(shape = no_variation), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
  scale_color_brewer(palette = "Set2", name = NULL) +
  labs(
    x = "Spearman's correlation coefficient (\u03c1)", y = "District", col = NULL,
    caption = "\u00d7 = no correlation estimable (no variation in test positivity)"
  ) +
  base_theme +
  theme(plot.caption = element_text(size = 7, hjust = 0)) +
  guides(color = guide_legend(nrow = 1))


#Making Supplementary Figure 3
ggsave("FigS3.jpg", FigS3, width = 25, height = 30, dpi = 600, units = "cm")


# ============================================================================
# Figure 4: ANC1 vs. population-level test positivity relationship
# ============================================================================


#Figure 4A
Fig4A <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  ungroup() %>%
  group_by(Province, Year, Month, Date) %>%
  summarise(
    `Total MiPs`               = sum(MiPs, na.rm = TRUE),
    `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
    `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100,
    across(c(`Test Positivity Rate`, `Population Test Positivity Rate`), ~replace_na(., 0)),
    Province = factor(Province, levels = province_order)
  ) %>%
  {
    cor_data  <- cor.test(.$`Test Positivity Rate`, .$`Population Test Positivity Rate`,
                          method = "spearman")
    cor_label <- sprintf("Spearman rho = %.3f\np < 0.001", cor_data$estimate)
    
    ggplot(., aes(x = `Test Positivity Rate`, y = `Population Test Positivity Rate`, colour = Province)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "black",
                  linetype = "dashed", linewidth = 0.8, alpha = 0.2) +
      scale_color_brewer(palette = "Set2", name = NULL) +
      annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.5, size = 3.5) +
      labs(tag = "A", x = "ANC1 test positivity (%)",
           y = "Population-level test positivity (%)", col = NULL) +
      base_theme +
      theme(
        legend.margin        = margin(t = -2),
        legend.box.margin    = margin(t = -4),
        legend.key.height    = unit(8, "pt"),
        legend.key.spacing.y = unit(1, "pt"),
        legend.text          = element_text(margin = margin(t = 0, b = 0))
      )
  }


# Alternative version with per-province + overall regression lines (not used in main figure)
Fig4A_alt <- Monthly_Data_Complete %>%
  group_by(Province, Year, Month, Date) %>%
  summarise(
    `Total MiPs`               = sum(MiPs, na.rm = TRUE),
    `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
    `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100,
    across(c(`Test Positivity Rate`, `Population Test Positivity Rate`), ~replace_na(., 0)),
    Province = factor(Province, levels = province_order)
  ) %>%
  {
    cor_data  <- cor.test(.$`Test Positivity Rate`, .$`Population Test Positivity Rate`,
                          method = "spearman")
    cor_label <- sprintf("Spearman rho = %.3f\np < 0.001", cor_data$estimate)
    
    ggplot(., aes(x = `Test Positivity Rate`, y = `Population Test Positivity Rate`, colour = Province)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, aes(group = Province, colour = Province),
                  linetype = "dashed", linewidth = 0.8, alpha = 0.2) +
      geom_smooth(method = "lm", se = TRUE, colour = "black",
                  linetype = "dashed", linewidth = 1, alpha = 0.2) +
      scale_color_brewer(palette = "Set2", name = "Province") +
      annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.5, size = 3.5) +
      labs(x = "ANC1 test positivity (%)", y = "Population test positivity (%)") +
      base_theme +
      theme(legend.position = "right", plot.title = element_text(size = 10, face = "bold"))
  }


#Figure 4B
Fig4B <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  ungroup() %>%
  group_by(Year, Month, Date) %>%
  summarise(
    `Total MiPs`               = sum(MiPs, na.rm = TRUE),
    `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
    `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100
  ) %>%
  {
    y1_range <- range(.$`Test Positivity Rate`, na.rm = TRUE)
    y2_range <- range(.$`Population Test Positivity Rate`, na.rm = TRUE)
    scale_factor  <- diff(y2_range) / diff(y1_range)
    y2_shift      <- y2_range[1] - y1_range[1] * scale_factor
    se_positivity <- sd(.$`Test Positivity Rate`, na.rm = TRUE) / sqrt(nrow(.))
    
    ggplot(., aes(x = Date)) +
      geom_ribbon(aes(ymin = `Test Positivity Rate` - se_positivity,
                      ymax = `Test Positivity Rate` + se_positivity),
                  alpha = 0.2, fill = "#e41a1c") +
      geom_line(aes(y = `Test Positivity Rate`, color = "ANC1 test positivity"), linewidth = 1) +
      geom_line(aes(y = (`Population Test Positivity Rate` - y2_shift) / scale_factor,
                    color = "Population-level test positivity"), linewidth = 1) +
      scale_y_continuous(
        name = "ANC1 test positivity (%)",
        sec.axis = sec_axis(~ . * scale_factor + y2_shift, name = "Population-level test positivity (%)")
      ) +
      scale_color_manual(name = "", values = c(
        "ANC1 test positivity" = "#e41a1c", "Population-level test positivity" = "#377eb8"
      )) +
      scale_x_date(
        breaks = seq(ymd("2019-07-01"), ymd("2023-07-01"), by = "1 year"),
        date_labels = "%Y", expand = c(0.01, 0.01), minor_breaks = NULL
      ) +
      geom_vline(
        xintercept = seq(ymd("2019-01-01"), ymd("2023-12-01"), by = "1 year"),
        linetype = "dashed", colour = "grey60", linewidth = 0.3
      ) +
      labs(tag = "B", x = "") +
      base_theme +
      theme(
        legend.margin     = margin(t = -2),
        legend.box.margin = margin(t = -4),
        axis.title.y.left  = element_text(color = "#e41a1c"),
        axis.title.y.right = element_text(color = "#377eb8"),
        axis.text.y.left   = element_text(color = "#e41a1c"),
        axis.text.y.right  = element_text(color = "#377eb8")
      )
  }


#Figure 4C
plot_df <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  ungroup() %>%
  group_by(Province, Year, Month, Date) %>%
  summarise(
    `Total MiPs`               = sum(MiPs, na.rm = TRUE),
    `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
    `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100,
    across(c(`Test Positivity Rate`, `Population Test Positivity Rate`), ~replace_na(., 0)),
    Province = factor(Province, levels = province_order)
  )

correlation    <- plot_df %>% group_by(Province) %>%
  summarise(rho = cor(`Test Positivity Rate`, `Population Test Positivity Rate`,
                      method = "spearman", use = "complete.obs"), .groups = "drop")

scale_factor_df <- plot_df %>% group_by(Province) %>%
  summarise(sf = max(`Test Positivity Rate`, na.rm = TRUE) /
              max(`Population Test Positivity Rate`, na.rm = TRUE), .groups = "drop")

plot_df <- plot_df %>%
  left_join(correlation, by = "Province") %>%
  left_join(scale_factor_df, by = "Province") %>%
  mutate(facet_label = sprintf("%s\n\u03c1 = %.2f", Province, rho))

make_panel <- function(prov_data, position = c("left", "middle", "right")) {
  position <- match.arg(position)
  sf_i <- unique(prov_data$sf)
  
  p <- ggplot(prov_data, aes(x = Date)) +
    geom_line(aes(y = `Test Positivity Rate`, colour = "ANC1 prevalence"), linewidth = 0.7) +
    geom_line(aes(y = `Population Test Positivity Rate` * sf_i, colour = "Population-level prevalence"),
              linewidth = 0.7) +
    scale_color_manual(name = "", values = c(
      "ANC1 prevalence" = "#e41a1c", "Population-level prevalence" = "#377eb8"
    )) +
    scale_y_continuous(name = "ANC1 (%)", sec.axis = sec_axis(~ . / sf_i, name = "Population-level (%)")) +
    scale_x_date(
      breaks = seq(ymd("2019-07-01"), ymd("2023-07-01"), by = "1 year"),
      date_labels = "%Y", expand = c(0.01, 0.01), minor_breaks = NULL
    ) +
    geom_vline(
      xintercept = seq(ymd("2019-01-01"), ymd("2023-12-01"), by = "1 year"),
      linetype = "dashed", colour = "grey60", linewidth = 0.3
    ) +
    labs(title = unique(prov_data$facet_label), x = NULL) +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.line   = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.4, colour = "black"),
      panel.grid  = element_blank(),
      legend.position = "none",
      plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title.y.left  = element_text(color = "#e41a1c", size = 9),
      axis.title.y.right = element_text(color = "#377eb8", size = 9),
      axis.text.y = element_text(size = 7)
    )
  
  if (position == "left") {
    p <- p + theme(axis.title.y.right = element_blank())
  } else if (position == "right") {
    p <- p + theme(axis.title.y.left = element_blank())
  } else {
    p <- p + theme(axis.title.y.left = element_blank(), axis.title.y.right = element_blank())
  }
  p
}

prov_list  <- plot_df %>% group_split(Province)
ncol_grid  <- 3
positions  <- sapply(seq_along(prov_list), function(i) {
  col_idx <- ((i - 1) %% ncol_grid) + 1
  if (col_idx == 1) "left" else if (col_idx == ncol_grid) "right" else "middle"
})

panels <- map2(prov_list, positions, make_panel)

Fig4C <- wrap_plots(panels, ncol = ncol_grid) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


#Making Figure 4
Fig4 <- (Fig4A | Fig4B) / Fig4C +
  plot_annotation(tag_levels = list(c("A", "B", "C")))
ggsave("Fig4.jpg", Fig4, width = 20, height = 20, dpi = 600, units = "cm")


# ============================================================================
# Supplementary Figure 4: cross-correlation lag + seasonal decomposition
# ============================================================================


#Figure S4A
FigS4A <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  ungroup() %>%
  group_by(Year, Month, Date) %>%
  summarise(
    `Total MiPs`               = sum(MiPs, na.rm = TRUE),
    `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
    `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100
  ) %>%
  {
    max_lag <- 6  # months
    ccf_results <- data.frame(lag = -max_lag:max_lag, correlation = numeric(2 * max_lag + 1))
    
    for (i in seq_len(nrow(ccf_results))) {
      lag_val <- ccf_results$lag[i]
      if (lag_val < 0) {
        mip_data <- head(.$`Test Positivity Rate`, -abs(lag_val))
        pop_data <- tail(.$`Population Test Positivity Rate`, -abs(lag_val))
      } else if (lag_val > 0) {
        mip_data <- tail(.$`Test Positivity Rate`, -lag_val)
        pop_data <- head(.$`Population Test Positivity Rate`, -lag_val)
      } else {
        mip_data <- .$`Test Positivity Rate`
        pop_data <- .$`Population Test Positivity Rate`
      }
      ccf_results$correlation[i] <- cor(mip_data, pop_data, method = "spearman", use = "complete.obs")
    }
    
    ggplot(ccf_results, aes(x = lag, y = correlation)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linewidth = 0.8, col = "red", linetype = 2) +
      geom_line(linewidth = 1, color = "#377eb8") +
      geom_point(size = 3, color = "#377eb8") +
      scale_x_continuous(breaks = -max_lag:max_lag) +
      labs(
        tag = "A",
        x = "Lag (months)\nNegative: MiP leads population | Positive: Population leads MiP",
        y = "Spearman's correlation coefficient"
      ) +
      base_theme
  }


#Figure S4B
FigS4B <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  ungroup() %>%
  group_by(Year, Month, Date) %>%
  summarise(
    `Total MiPs`               = sum(MiPs, na.rm = TRUE),
    `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
    `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100
  ) %>%
  {
    ts_mip <- ts(.$`Test Positivity Rate`, start = c(2019, 1), frequency = 12)
    ts_pop <- ts(.$`Population Test Positivity Rate`, start = c(2019, 1), frequency = 12)
    
    decomp_mip <- stl(ts_mip, s.window = "periodic")
    decomp_pop <- stl(ts_pop, s.window = "periodic")
    
    decomp_data <- bind_rows(
      data.frame(date = .$Date, metric = "ANC1 prevalence",
                 observed = as.numeric(ts_mip),
                 trend    = as.numeric(decomp_mip$time.series[, "trend"]),
                 seasonal = as.numeric(decomp_mip$time.series[, "seasonal"]),
                 remainder = as.numeric(decomp_mip$time.series[, "remainder"])),
      data.frame(date = .$Date, metric = "Population-level prevalence",
                 observed = as.numeric(ts_pop),
                 trend    = as.numeric(decomp_pop$time.series[, "trend"]),
                 seasonal = as.numeric(decomp_pop$time.series[, "seasonal"]),
                 remainder = as.numeric(decomp_pop$time.series[, "remainder"]))
    ) %>%
      pivot_longer(cols = c(observed, trend, seasonal, remainder),
                   names_to = "component", values_to = "value") %>%
      mutate(component = factor(component,
                                levels = c("observed", "trend", "seasonal", "remainder"),
                                labels = c("Observed", "Trend", "Seasonal", "Remainder")))
    
    ggplot(decomp_data, aes(x = date, y = value, color = metric)) +
      geom_line(linewidth = 0.7) +
      scale_color_manual(values = c("#e41a1c", "#377eb8"), name = "") +
      facet_wrap(~ component, scales = "free_y", ncol = 1) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(tag = "B", x = "Year", y = "Value") +
      base_theme
  }


#Making Supplementary Figure 4
FigS4 <- FigS4A / FigS4B
ggsave("FigS4.jpg", FigS4, width = 25, height = 30, dpi = 600, units = "cm")


# ============================================================================
# Summary Statistics and Tables
# ============================================================================


#Table S1: District-level summary statistics
TableS1 <- Yearly_Data_Complete %>%
  filter(Malaria_Testing_ANC >= 30) %>%
  group_by(Province, Regency) %>%
  summarise(
    `Total Pregnancies`         = sum(Pregnancies_imputed, na.rm = TRUE),
    `Total ANC1 Visits`         = sum(ANC1_imputed, na.rm = TRUE),
    `Total ANC1 Testing`        = sum(Malaria_Testing_ANC, na.rm = TRUE),
    `Total MiPs`                = sum(MiPs, na.rm = TRUE),
    `Total Population Testing`  = sum(Malaria_Population_Testing, na.rm = TRUE),
    `Total Malaria Positives`   = sum(`Malaria Positives`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `ANC1 Coverage`                    = (`Total ANC1 Visits` / `Total Pregnancies`) * 100,
    `ANC1 Testing Coverage`            = (`Total ANC1 Testing` / `Total Pregnancies`) * 100,
    `MiP Test Positivity Rate`         = (`Total MiPs` / `Total ANC1 Testing`) * 100,
    `Population Test Positivity Rate`  = (`Total Malaria Positives` / `Total Population Testing`) * 100
  ) %>%
  arrange(Province, Regency) %>%
  select(Province, Regency, `ANC1 Coverage`, `ANC1 Testing Coverage`,
         `MiP Test Positivity Rate`, `Population Test Positivity Rate`) %>%
  mutate(across(where(is.numeric), ~round(., 2)))


#Table S2: Species and age breakdown by province
TableS2 <- Regmal1_MiP %>%
  group_by(Province) %>%
  summarise(
    `Total MiPs`  = n(),
    P.falciparum  = sum(Parasite == "P.falciparum", na.rm = TRUE),
    P.vivax       = sum(Parasite == "P.vivax", na.rm = TRUE),
    P.mix         = sum(Parasite == "P.Mix", na.rm = TRUE),
    P.malariae    = sum(Parasite == "P.malariae", na.rm = TRUE),
    P.ovale       = sum(Parasite == "P.ovale", na.rm = TRUE),
    `10-19`      = sum(Age_Group == "10-19", na.rm = TRUE),
    `20-34`      = sum(Age_Group == "20-34", na.rm = TRUE),
    `35+`        = sum(Age_Group == "35+", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Pf %`   = P.falciparum / `Total MiPs` * 100,
    `Pv %`   = P.vivax / `Total MiPs` * 100,
    `PMix %` = P.mix / `Total MiPs` * 100,
    `Pm %`  = P.malariae / `Total MiPs` * 100,
    `Po %`  = P.ovale / `Total MiPs` * 100,
    `Adolescent Pregnancy %` = `10-19` / `Total MiPs` * 100,
    `Optimal Pregnancy %`    = `20-34` / `Total MiPs` * 100,
    `Advanced Pregnancy %`   = `35+` / `Total MiPs` * 100,
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  ) %>%
  mutate(Province = factor(Province, levels = province_order)) %>%
  pivot_longer(
    cols = -Province,
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Province,
    values_from = Value,
    names_sort = FALSE
  ) %>%
  select(
    Variable,
    all_of(province_order)
  ) %>% 
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Total MiPs", "P.falciparum", "P.vivax", "P.mix",
        "P.malariae", "P.ovale",
        "Pf %", "Pv %", "PMix %", "Pm %", "Po %",
        "10-19", "20-34", "35+",
        "Adolescent Pregnancy %", "Optimal Pregnancy %",
        "Advanced Pregnancy %"
      )
    )) %>% 
  arrange(Variable) %>%
  rename(Province = Variable)
  

#Table S3: Correlation analysis results
TableS3 <- Monthly_Data_Complete %>%
  group_by(Province, Regency, Year) %>%
  filter(
    sum(Malaria_Testing_ANC, na.rm = TRUE) >= 30,
    sum(Malaria_Population_Testing, na.rm = TRUE) >= 30
  ) %>%
  ungroup() %>%
  {
    prov_monthly <- (.) %>%
      group_by(Province, Year, Month, Date) %>%
      summarise(
        `Total MiPs`               = sum(MiPs, na.rm = TRUE),
        `Total ANC1 Testing`       = sum(Malaria_Testing_ANC, na.rm = TRUE),
        `Total Positives`          = sum(`Malaria Positives`, na.rm = TRUE),
        `Total Population Testing` = sum(Malaria_Population_Testing, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        `Test Positivity Rate`            = (`Total MiPs` / `Total ANC1 Testing`) * 100,
        `Population Test Positivity Rate` = (`Total Positives` / `Total Population Testing`) * 100,
        across(c(`Test Positivity Rate`, `Population Test Positivity Rate`), ~replace_na(., 0))
      )
    
    provincial <- prov_monthly %>%
      group_by(Province) %>%
      summarise(
        n_months = n(),
        rho = cor(`Test Positivity Rate`, `Population Test Positivity Rate`,
                  method = "spearman", use = "complete.obs"),
        .groups = "drop"
      ) %>%
      mutate(
        level = "Provincial",
        se = sqrt((1 - rho^2) / (n_months - 2)),
        ci_lower = tanh(atanh(rho) - 1.96 * se),
        ci_upper = tanh(atanh(rho) + 1.96 * se),
        p_value  = 2 * pt(-abs(rho * sqrt(n_months - 2) / sqrt(1 - rho^2)), df = n_months - 2)
      ) %>%
      select(level, Province, n = n_months, rho, ci_lower, ci_upper, p_value)
    
    district <- (.) %>%
      mutate(across(c(`Test Positivity Rate`, `Population Test Positivity Rate`), ~replace_na(., 0))) %>%
      group_by(Province, Regency) %>%
      summarise(
        n = n(),
        rho = cor(`Test Positivity Rate`, `Population Test Positivity Rate`,
                  method = "spearman", use = "complete.obs"),
        .groups = "drop"
      ) %>%
      mutate(
        level = "District",
        se = sqrt((1 - rho^2) / (n - 2)),
        ci_lower = tanh(atanh(rho) - 1.96 * se),
        ci_upper = tanh(atanh(rho) + 1.96 * se)
      ) %>%
      select(level, Province, Regency, n, rho, ci_lower, ci_upper) %>%
      arrange(Province, rho)
    
    bind_rows(provincial, district)
  } %>% 
  mutate(across(where(is.numeric), ~round(., 2)))


#Making the supplementary tables
export(list(
  "Table S1" = TableS1,
  "Table S2" = TableS2,
  "Table S3" = TableS3
), file = "Supplementary Tables.xlsx")


# ============================================================================
# SESSION INFO
# ============================================================================


# Print session info for reproducibility
cat("\n=== Session Information ===\n")
sessionInfo()


cat("\n=== Analysis Complete ===\n")
cat("\nMain figures saved:\n")
cat("  - Figure1_Geographic_Temporal_Overview.png\n")
cat("  - Figure2_Case_Characteristics.png\n")
cat("  - Figure3_Surveillance_Performance.png\n")
cat("  - Figure4_MiP_Population_Comparison.png\n")


cat("\nSupplementary figures saved:\n")
cat("  - FigureS1_District_Temporal_Details.png\n")
cat("  - FigureS2_Species_Distribution.png\n")
cat("  - FigureS3a_District_Correlations.png\n")
cat("  - FigureS3b_Scatter_Matrix.png\n")
cat("  - FigureS4_Temporal_Sensitivity.png\n")


cat("\nSupplementary tables saved:\n")
cat("  - Table S1: District-level summary statistics\n")
cat("  - Table S2: Species and age breakdown by province\n")
cat("  - Table S3: Correlation analysis results\n")
cat("  - Supplementary Tables.xlsx\n")
cat("\nOutput files generated successfully.\n")


# ============================================================================
# END OF SCRIPT
# ============================================================================