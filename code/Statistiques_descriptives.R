library(readxl)
library(tidyverse)
library(ggplot2)
library(patchwork)

#######################################
# Variables cliniques
#######################################

#setwd("")
data_DEXA <- read_excel("data/caractéristiques cliniques et sanguines/AGBRESA_BCD_DEXA_TotalBodyComp_2_corrected180320202.xlsx")

data_DEXA <- data_DEXA %>%
  select(-Study, -gBMC, -gLeanBMC, -PercFat) %>%
  filter(Region == "Total") %>%
  filter(Sujet %in% c("F", "J", "W", "D", "T")) %>%
  mutate(
    PercFat  = gFat * 100 / gTotal,
    PercLean = gLean * 100 / gTotal
  )

# attribution échelle
data_DEXA <- data_DEXA %>%
  mutate(Jour = case_when(
    Examination == "BDC-13" ~ 0,
    Examination == "BDC-4"  ~ 9,
    Examination == "HDT15"  ~ 28,
    Examination == "HDT30"  ~ 43,
    Examination == "HDT45"  ~ 58,
    Examination == "HDT60"  ~ 73,
    Examination == "R+11"   ~ 84
  ))

# plot masse totale
p_masse_totale <- ggplot(data_DEXA, aes(x = Jour, y = gTotal/1000, group = Sujet, 
                      color = Sujet, shape = Sujet, linetype = Sujet)) +
  geom_vline(xintercept = c(13, 73), 
             color = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0, 9, 28, 43, 58, 73, 84),
                     labels = c("BDC-13", "BDC-4", "HDT15", "HDT30", "HDT45", "HDT60", "R+11")) +
  scale_color_manual(values = c(
    "D" = "grey55",
    "F" = "grey20",
    "J" = "grey20",
    "T" = "grey35",
    "W" = "grey10",
    "Z" = "grey35"
  )) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Evolution de la masse des sujets au cours du protocole",
       x = "Temps (condition)", y = "Masse (kg)", 
       color = "Sujet", shape = "Sujet", linetype = "Sujet")

p_masse_totale

# plot masse maigre
p_masse_maigre <- ggplot(data_DEXA, aes(x = Jour, y = PercLean, group = Sujet, 
                      color = Sujet, shape = Sujet, linetype = Sujet)) +
  geom_vline(xintercept = c(13, 73), 
             color = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0, 9, 28, 43, 58, 73, 84),
                     labels = c("BDC-13", "BDC-4", "HDT15", "HDT30", "HDT45", "HDT60", "R+11")) +
  scale_color_manual(values = c(
    "D" = "grey55",
    "F" = "grey20",
    "J" = "grey20",
    "T" = "grey35",
    "W" = "grey10",
    "Z" = "grey35"
  )) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Evolution de la masse maigre des sujets au cours du protocole",
       x = "Temps (condition)", y = "Masse maigre (% masse totale)", 
       color = "Sujet", shape = "Sujet", linetype = "Sujet")

# plot masse grasse
p_masse_grasse <- ggplot(data_DEXA, aes(x = Jour, y = PercFat, group = Sujet, 
                      color = Sujet, shape = Sujet, linetype = Sujet)) +
  geom_vline(xintercept = c(13, 73), 
             color = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0, 9, 28, 43, 58, 73, 84),
                     labels = c("BDC-13", "BDC-4", "HDT15", "HDT30", "HDT45", "HDT60", "R+11")) +
  scale_color_manual(values = c(
    "D" = "grey55",
    "F" = "grey20",
    "J" = "grey20",
    "T" = "grey35",
    "W" = "grey10",
    "Z" = "grey35"
  )) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Evolution de la masse grasse des sujets au cours du protocole",
       x = "Temps (condition)", y = "Masse grasse (% masse totale)", 
       color = "Sujet", shape = "Sujet", linetype = "Sujet")

# Par sujet
data_DEXA %>%
  filter(Examination %in% c("BDC-4", "HDT60")) %>%
  group_by(Sujet, Examination) %>%
  summarise(gTotal = mean(gTotal), PercFat = mean(PercFat), PercLean = mean(PercLean))

# Tous sujets confondus
data_DEXA %>%
  filter(Examination %in% c("BDC-4", "HDT60")) %>%
  group_by(Examination) %>%
  summarise(
    gTotal_mean = mean(gTotal), gTotal_sd = sd(gTotal),
    PercFat_mean = mean(PercFat), PercFat_sd = sd(PercFat),
    PercLean_mean = mean(PercLean), PercLean_sd = sd(PercLean)
  )

#######################################
# Variables systémiques
#######################################

#setwd("")
data_dosage <- read_excel("data/caractéristiques cliniques et sanguines/2020 AVRILMAI Dosages bedrest Rendus-hsCRP.xlsx")

data_dosage <- data_dosage %>%
  select(Sujet, Condition, `Ferritine\n (µg/l)`, `CRP\n (mg/l)`, `Haptoglobine \n(g/l)`) %>%
  filter(Sujet %in% c("F", "J", "W", "D", "T"))

# attribution échelle
data_dosage <- data_dosage %>%
  mutate(Jour = case_when(
    Condition == "BDC 6"  ~ 0,
    Condition == "HDT 6"  ~ 12,
    Condition == "HDT 57" ~ 63
  ))

# plot ferritine
p_ferritine <- ggplot(data_dosage, aes(x = Jour, y = `Ferritine\n (µg/l)`, group = Sujet, 
                        color = Sujet, shape = Sujet, linetype = Sujet)) +
  geom_vline(xintercept = 6, 
             color = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0, 12, 63), labels = c("BDC 6", "HDT 6", "HDT 57")) +
  scale_color_manual(values = c(
    "D" = "grey55", "F" = "grey20", "J" = "grey20",
    "T" = "grey35", "W" = "grey10"
  )) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  labs(title = "Evolution de la concentration en ferritine des sujets au cours du protocole",
       x = "Temps (condition)", y = "Ferritine (µg/L)",
       color = "Sujet", shape = "Sujet", linetype = "Sujet")

# plot CRP
p_CRP <- ggplot(data_dosage, aes(x = Jour, y = `CRP\n (mg/l)`, group = Sujet, 
                        color = Sujet, shape = Sujet, linetype = Sujet)) +
  geom_vline(xintercept = 6, 
             color = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0, 12, 63), labels = c("BDC 6", "HDT 6", "HDT 57")) +
  scale_color_manual(values = c(
    "D" = "grey55", "F" = "grey20", "J" = "grey20",
    "T" = "grey35", "W" = "grey10"
  )) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  labs(title = "Evolution de la concentration en CRP des sujets au cours du protocole",
       x = "Temps (condition)", y = "CRP (mg/L)",
       color = "Sujet", shape = "Sujet", linetype = "Sujet")

# plot Haptoglobine
p_haptoglobine <- ggplot(data_dosage, aes(x = Jour, y = `Haptoglobine \n(g/l)`, group = Sujet, 
                        color = Sujet, shape = Sujet, linetype = Sujet)) +
  geom_vline(xintercept = 6, 
             color = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0, 12, 63), labels = c("BDC 6", "HDT 6", "HDT 57")) +
  scale_color_manual(values = c(
    "D" = "grey55", "F" = "grey20", "J" = "grey20",
    "T" = "grey35", "W" = "grey10"
  )) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  labs(title = "Evolution de la concentration en haptoglobine des sujets au cours du protocole",
       x = "Temps (condition)", y = "Haptoglobine (g/L)",
       color = "Sujet", shape = "Sujet", linetype = "Sujet")

# Par sujet
data_dosage %>%
  filter(Condition %in% c("BDC 6", "HDT 6")) %>%
  group_by(Sujet, Condition) %>%
  summarise(Ferritine = mean(`Ferritine\n (µg/l)`), CRP = mean(`CRP
 (mg/l)`), Haptoglobine = mean(`Haptoglobine 
(g/l)`))

# Tous sujets confondus
data_dosage %>%
  filter(Condition %in% c("BDC 6", "HDT 6")) %>%
  group_by(Condition) %>%
  summarise(
    Ferritine_mean = mean(`Ferritine
 (µg/l)`), Ferritine_sd = sd(`Ferritine
 (µg/l)`),
    CRP_mean = mean(`CRP
 (mg/l)`), CRP_sd = sd(`CRP
 (mg/l)`),
    Haptoglobine_mean = mean(`Haptoglobine 
(g/l)`), Haptoglobine_sd = sd(`Haptoglobine 
(g/l)`)
  )


# Export pdf
combined_plot <- p_masse_totale / p_masse_maigre / p_masse_grasse
ggsave("masse_totale_maigre_grasse.pdf", 
       plot = combined_plot, 
       width = 8, height = 10)

combined_plot <- p_ferritine / p_CRP / p_haptoglobine
ggsave("masse_ferritine_CRP_haptoglobine.pdf", 
       plot = combined_plot, 
       width = 8, height = 10)