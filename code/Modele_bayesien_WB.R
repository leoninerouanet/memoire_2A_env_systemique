library(readxl)
library(ggplot2)
library(dplyr)
library(R2jags)
library(patchwork)
library(tidybayes)
library(ggdist)
library(tidyr)
library(stringr)
library(HDInterval)
library(bayesplot)
library(ggdag)
library(knitr)
library(kableExtra)

set.seed(123)

# Préparation des données

#setwd("")
data <- read_excel("data/western blot/Quantification_WB.xlsx")

data <- data %>%
  mutate(Signal_moyen = rowMeans(across(c(n1, n2)), na.rm = TRUE))

data_wide <- data %>%
  select(Sujet, Condition, Protéine, Signal_moyen) %>%
  pivot_wider(
    names_from = Protéine,
    values_from = Signal_moyen
  )

data_fc <- data_wide %>%
  mutate(
    across(
      -c(Sujet, Condition),
      ~ .x / mean(.x[Condition == "BDC-6"], na.rm = TRUE)
    )
  )

data_fc$condition_id <- ifelse(data_fc$Condition == "HDT6", 1, 0)
data_fc$sujet_id <- as.numeric(as.factor(data_fc$Sujet))

proteines <- setdiff(colnames(data_fc), c("Sujet", "Condition", "condition_id", "sujet_id"))

# Initialisation dataframe final
df_draws_final <- data.frame()
model_fits     <- list()   
ppc_plots      <- list()

# Modèle JAGS
model_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- alpha[subject[i]] + beta * condition[i]
    y_sim[i] ~ dlnorm(mu[i], tau) # pour posterior predtictive check

  }
  
  for (s in 1:N_sujets) {
    alpha[s] ~ dnorm(mu_alpha_log, tau_alpha)
  }
  
  # Priors
  beta ~ dnorm(0, 3.3)
  mu_alpha_log <- 0 
  sigma ~ dnorm(0,1) T(0,)
  tau <- 1 / pow(sigma, 2)
  sigma_alpha ~ dnorm(0,1) T(0,)
  tau_alpha <- 1 / pow(sigma_alpha, 2)
  
  # --- Quantités dérivées ---
  ratio_effet <- exp(beta)
  
  # On compare le FC prédit en post (1) vs pré (0)
  # pour un sujet 'moyen' (exp(mu_alpha_log))
  mu_alpha_FC <- exp(mu_alpha_log)
  mu_alpha_post_FC <- exp(mu_alpha_log + beta)
}"

# Boucle sur chaque protéine
for (prot in proteines) {
  
  cat("\n===========================================\n")
  cat("Modèle pour la protéine :", prot, "\n")
  cat("===========================================\n")
  
  # Vérification que la colonne existe bien dans data_fc
  if (!prot %in% colnames(data_fc)) {
    warning(paste("Protéine introuvable dans data_fc :", prot, "— ignorée."))
    next
  }
  
  # Préparation des données pour cette protéine
  data_list_prot <- list(
    y         = data_fc[[prot]],
    condition = data_fc$condition_id,
    subject   = data_fc$sujet_id,
    N         = nrow(data_fc),
    N_sujets  = max(data_fc$sujet_id)
  )
  
  # Ajustement du modèle JAGS
  model_fit <- jags(
    data               = data_list_prot,
    parameters.to.save = c("alpha", "beta", "sigma", "sigma_alpha",
                           "ratio_effet", "mu_alpha_FC",
                           "mu_alpha_log", "mu_alpha_post_FC", "y_sim"),
    model.file         = textConnection(model_string),
    n.chains           = 4,
    n.iter             = 10000,
    n.burnin           = 2000
  )
  
  print(model_fit)
  model_fits[[prot]] <- model_fit
  # Extraction des draws posterior
  draws_mat  <- model_fit$BUGSoutput$sims.matrix
  post_draws <- draws_mat %>% spread_draws(c(ratio_effet, sigma_alpha))
  
  # Empilement dans df_draws_final
  df_temp <- as.data.frame(post_draws) %>%
    mutate(proteine = prot)
  
  df_draws_final <- bind_rows(df_draws_final, df_temp)
  
  # Posterior Predictive Check
  y_rep <- model_fit$BUGSoutput$sims.list$y_sim
  y_obs <- data_list_prot$y
  
  ppc_plots[[prot]] <- ppc_dens_overlay(y_obs, y_rep[1:100, ]) +
    ggtitle(paste("PPC —", prot)) +
    labs(x = "Fold Change", y = "Densité") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

# save(df_draws_final, file = "resultats_model.RData")
# load("resultats_model.RData")

# plot posterior predictive check
wrap_plots(ppc_plots, ncol = 2)

# Résultats

# Plot

draws_wide <- df_draws_final %>%
  select(.draw, ratio_effet, proteine) %>%
  pivot_wider(names_from  = proteine,
              values_from = ratio_effet) %>%
  select(-.draw) %>%
  as.matrix()

bayesplot::color_scheme_set("gray")

draws_wide <- draws_wide[, c("Myosine", "P-Akt", "P-4EBP1", "P-p70s6k", "Ubiquitine", "P-AMPK", "P-p38", "PGC1-alpha")]

WB <- mcmc_areas(
  draws_wide,
  prob       = 0.95,
  prob_outer = 0.99,
  point_est  = "median"
) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 1.2, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 1, color = "grey50", linewidth = 0.4) +
  labs(
    title    = "Probabilité de l'effet de la condition HDT6 par\nprotéine d'intérêt",
    x        = "Effet de la condition HDT6 (ratio HDT6/BDC-6)",
    y        = "Protéine d'intérêt"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plot(WB)

# ggsave(file="WB.pdf", plot=WB, width=4.7, height=3.8)

# Tableau statistiques récapitulatives par protéine

ROPE_low  <- 0.8
ROPE_high <- 1.2

table_recap <- df_draws_final %>%
  group_by(proteine) %>%
  summarise(
    Mediane     = median(ratio_effet),
    IC_low_95   = quantile(ratio_effet, 0.025),
    IC_high_95  = quantile(ratio_effet, 0.975),
    P_direction = max(mean(ratio_effet > 1), mean(ratio_effet < 1)),
    P_ROPE      = mean(ratio_effet >= ROPE_low & ratio_effet <= ROPE_high),
    P_sup_ROPE  = mean(ratio_effet >= ROPE_high),
    P_inf_ROPE  = mean(ratio_effet <= ROPE_low ),
    Sigma_alpha = round(mean(sigma_alpha),2),
    .groups     = "drop"
  ) %>%
  mutate(
    IC_95          = paste0("[", round(IC_low_95, 2), " ; ", round(IC_high_95, 2), "]"),
    Mediane        = round(Mediane, 2),
    P_direction    = round(P_direction, 3),
    P_ROPE         = round(P_ROPE, 3),
    P_sup_ROPE     = round(P_sup_ROPE, 3),
    P_inf_ROPE     = round(P_inf_ROPE, 3)
  ) %>%
  select(Protéine = proteine, Mediane, IC_95, P_direction, P_ROPE, P_inf_ROPE, P_sup_ROPE, Sigma_alpha)

print(table_recap)

# Affichage avec kable
table_recap %>%
  kable(
    format  = "latex",
    caption = paste0("Résumé des distributions postérieures du ratio HDT6/BDC-6 ",
                     "(ROPE = [", ROPE_low, " ; ", ROPE_high, "])"),
    align   = c("l", "c", "c", "c", "c", "c", "c", "c"),
    col.names = c("Protéine", "Médiane", "IC 95%",
                  "P(direction)", "P(ROPE)", "P(ROPE < 0.8)", "P(ROPE > 1.2)", "Sigma Alpha")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "bordered"),
    full_width        = FALSE,
    position          = "center"
  ) %>%
  column_spec(1, bold = TRUE) 

