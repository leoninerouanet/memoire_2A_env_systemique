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

set.seed(123)

# Préparation des données

#setwd("")
data_n1 <- read_excel("data/immunocytochimie/mesures_myotubes_n1.xlsx") %>%
  select(Moyenne_fibre, nom_original) %>%
  drop_na(Moyenne_fibre) %>%
  mutate(
    condition_raw = str_extract(nom_original, "^[^_]+"),
    sujet = str_extract(nom_original, "(?<=_)[^_]+"),
    diam_um = Moyenne_fibre * 100 / 58130.5
  ) %>%
  filter(!condition_raw %in% c("CMdiff", "CMS")) %>%
  mutate(
    condition_id = ifelse(condition_raw == "HDT6", 1, 0),
    sujet_id = as.numeric(as.factor(sujet))
  )

data_n2 <- read_excel("/data/immunocytochimie/mesures_myotubes_n2.xlsx") %>%
  select(Moyenne_fibre, nom_original) %>%
  drop_na(Moyenne_fibre) %>%
  mutate(
    condition_raw = str_extract(nom_original, "^[^_]+"),
    sujet = str_extract(nom_original, "(?<=_)[^_]+"),
    diam_um = Moyenne_fibre * 100 / 219.83
  ) %>%
  filter(!condition_raw %in% c("CMdiff", "CMS")) %>%
  mutate(
    condition_id = ifelse(condition_raw == "HDT6", 1, 0),
    sujet_id = as.numeric(as.factor(sujet))
  )

data_clean <- bind_rows(data_n1, data_n2)

# prior check
n_sim <- 100000
# Simulation des paramètres à partir des distributions priors
p_beta <- rnorm(n_sim, 0, 0.275) 
p_mu_alpha_log <- rnorm(n_sim, 3.4, 0.5) # exp3,4 =30
p_sigma <- abs(rnorm(n_sim, 0, 1)) 
p_sigma_alpha <- abs(rnorm(n_sim, 0, 1))
# Simulation des effets sujets (alphas)
p_alpha <- rnorm(n_sim, p_mu_alpha_log, p_sigma_alpha)
# Simulation de données (diamètre en µm)
# On simule un diamètre pour la condition 0 et un pour la 1
p_diam_0 <- rlnorm(n_sim, p_alpha, p_sigma)
p_diam_1 <- rlnorm(n_sim, p_alpha + p_beta, p_sigma)
# Dataframe pour le plot
prior_data <- data.frame(
  Condition = rep(c("BDC-4", "HDT6"), each = n_sim),
  Diametre = c(p_diam_0, p_diam_1)
)

# Visualisation
ggplot(prior_data, aes(x = Diametre)) +
  geom_density(alpha = 0.5, fill = "grey55") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Distribution des diamètres simulés sous les priors",
       x = "Diamètre (µm)",
       y = "Densité") +
  xlim(0, 100)

data_list <- list(
  d = data_clean$diam_um,
  condition = data_clean$condition_id,
  subject = data_clean$sujet_id,
  N = nrow(data_clean),
  N_sujets = max(data_clean$sujet_id)
)

# Modèle JAGS
model_string <- "
model {
  for (i in 1:N) {
    d[i] ~ dlnorm(mu[i], tau)
    mu[i] <- alpha[subject[i]] + beta * condition[i]
    #d_sim[i] ~ dlnorm(mu[i], tau) # pour posterior predtictive check

  }
  
  for (s in 1:N_sujets) {
    alpha[s] ~ dnorm(mu_alpha_log, tau_alpha)
  }
  
  # Priors
  beta ~ dnorm(0, 13)
  mu_alpha_log ~ dnorm(3.4, 4) 
  sigma ~ dnorm(0,1) T(0,)
  tau <- 1 / pow(sigma, 2)
  sigma_alpha ~ dnorm(0,1) T(0,)
  tau_alpha <- 1 / pow(sigma_alpha, 2)
  
  # --- Quantités dérivées ---
  ratio_effet <- exp(beta)
  
  # Différence absolue moyenne en µm
  # On compare le diamètre prédit en post (1) vs pré (0)
  # pour un sujet 'moyen' (exp(mu_alpha_log))
  diff_abs_um <- exp(mu_alpha_log + beta) - exp(mu_alpha_log)
  mu_alpha_um <- exp(mu_alpha_log)
  mu_alpha_post_um <- exp(mu_alpha_log + beta)
}"

model_fit <- jags(data = data_list, 
                  parameters.to.save = c("alpha", "beta", "sigma", "sigma_alpha", "ratio_effet", "diff_abs_um", "mu_alpha_um", "mu_alpha_log", "mu_alpha_post_um"),#, "d_sim"),
                  model.file = textConnection(model_string),
                  n.chains = 4, n.iter = 10000, n.burnin = 2000)

print(model_fit)

draws_mat <- model_fit$BUGSoutput$sims.matrix
post_draws <- draws_mat %>% 
  spread_draws(beta, ratio_effet, diff_abs_um, mu_alpha_um, sigma_alpha, sigma, mu_alpha_log, mu_alpha_post_um)
draws <- as.data.frame(post_draws)

# Résultats

# VALEUR diamètre dans la population
median(draws$mu_alpha_um)
hdi(draws$mu_alpha_um, 0.95)

# VALEUR diamètre POST
median(draws$mu_alpha_post_um)
hdi(draws$mu_alpha_post_um, 0.95)

# VALEUR beta
hdi_95     <- hdi(draws$diff_abs_um, credMass = 0.95)
hdi_80     <- hdi(draws$diff_abs_um, credMass = 0.80)
median_val <- median(draws$diff_abs_um)

p_inf_rope <- sum(draws[,'diff_abs_um'] < -1)/length(draws[,"diff_abs_um"])

p_in_rope <- sum(draws[,'diff_abs_um'] < 1 & draws[,'diff_abs_um'] > -1)/length(draws[,"diff_abs_um"])

# variabilité inter-sujet
exp(median(draws$sigma_alpha))
hdi(exp(draws$sigma_alpha), 0.95)

# variabilité intra-sujet
exp(median(draws$sigma))
hdi(exp(draws$sigma))

# plot distribution pré post par sujet

compute_dens <- function(vals, sujet_i, condition_i, y_center, sign) {
  d <- density(vals, bw = "SJ", n = 1024)
  data.frame(x = d$x, y = y_center + sign * d$y / max(d$y) * 0.4,
             sujet = sujet_i, condition = condition_i)
}

mcmc_samples <- as.matrix(model_fit$BUGSoutput$sims.matrix)
alpha_samples <- mcmc_samples[, grep("^alpha\\[", colnames(mcmc_samples))]
alpha_samples <- exp(alpha_samples)

sujets <- c("alpha[1]" = "D", "alpha[2]" = "F", "alpha[3]" = "J",
            "alpha[4]" = "T", "alpha[5]" = "W")

hdi_data <- do.call(rbind, lapply(names(sujets), function(par) {
  vals <- alpha_samples[, par]  
  data.frame(
    sujet      = sujets[par],
    hdi_95_lo  = hdi(vals, credMass = 0.95)[1],
    hdi_95_hi  = hdi(vals, credMass = 0.95)[2],
    hdi_80_lo  = hdi(vals, credMass = 0.80)[1],
    hdi_80_hi  = hdi(vals, credMass = 0.80)[2],
    median_val = median(vals)
  )
}))

beta_samples <- mcmc_samples[, "beta"]

alpha_post_samples <- exp(log(alpha_samples) + beta_samples)  

hdi_data_post <- do.call(rbind, lapply(seq_along(sujets), function(i) {
  vals <- alpha_post_samples[, i]
  data.frame(
    sujet      = sujets[i],
    hdi_95_lo  = hdi(vals, credMass = 0.95)[1],
    hdi_95_hi  = hdi(vals, credMass = 0.95)[2],
    hdi_80_lo  = hdi(vals, credMass = 0.80)[1],
    hdi_80_hi  = hdi(vals, credMass = 0.80)[2],
    median_val = median(vals)
  )
}))

sujet_labels <- c("F", "J", "W", "D", "T")
y_positions  <- seq(length(sujet_labels), 1)

# index dans alpha_samples (ordre original : D=1, F=2, J=3, T=4, W=5)
sujet_idx <- c(2, 3, 5, 1, 4)  # F, J, W, D, T

dens_data <- do.call(rbind, lapply(seq_along(sujet_labels), function(i) {
  rbind(
    compute_dens(alpha_samples[, sujet_idx[i]],      sujet_labels[i], "BDC",  y_positions[i],  1),
    compute_dens(alpha_post_samples[, sujet_idx[i]], sujet_labels[i], "HDT6", y_positions[i], -1)
  )
}))

med_data <- do.call(rbind, lapply(seq_along(sujet_labels), function(i) {
  rbind(
    data.frame(sujet = sujet_labels[i], condition = "BDC",
               x = median(alpha_samples[, sujet_idx[i]]),
               y0 = y_positions[i], y1 = y_positions[i] + 0.4),
    data.frame(sujet = sujet_labels[i], condition = "HDT6",
               x = median(alpha_post_samples[, sujet_idx[i]]),
               y0 = y_positions[i], y1 = y_positions[i] - 0.4)
  )
}))

mirror_plot <- ggplot() +
  geom_polygon(data = subset(dens_data, condition == "BDC"),
               aes(x = x, y = y, group = sujet, fill = condition), alpha = 0.5) +
  geom_polygon(data = subset(dens_data, condition == "HDT6"),
               aes(x = x, y = y, group = sujet, fill = condition), alpha = 0.5) +
  geom_hline(yintercept = y_positions, color = "grey60", linewidth = 0.3) +
  geom_segment(data = med_data,
               aes(x = x, xend = x, y = y0, yend = y1, color = condition),
               linewidth = 0.4) +
  scale_fill_manual(values  = c("BDC" = "grey80", "HDT6" = "grey30"),
                    labels  = c("BDC" = "BDC-6", "HDT6" = "HDT6"),
                    name = "Condition") +
  scale_color_manual(values = c("BDC" = "grey40", "HDT6" = "grey10"),
                     labels  = c("BDC" = "BDC-6", "HDT6" = "HDT6"),
                     name = "Condition") +
  scale_y_continuous(breaks = y_positions, labels = sujet_labels) +
  theme_bw() +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y       = element_blank()) +
  labs(x = "Diamètre moyen des myotubes (µm)", y = "Sujet",
       title = "Probabilité du diamètre moyen des\nmyotubes des sujets par condition")

plot(mirror_plot)

# plot effet avec rope

bayesplot::color_scheme_set("gray")

x_min <- min(draws$diff_abs_um)
x_max <- max(draws$diff_abs_um)

rope_lo <- -1
rope_hi <-  1

p1 <- ggplot(draws, aes(x = diff_abs_um)) +
  annotate("segment", x = rope_lo, xend = rope_lo, y = -Inf, yend = Inf,
           linetype = "dashed", color = "grey40", linewidth = 0.4) +
  annotate("segment", x = rope_hi, xend = rope_hi, y = -Inf, yend = Inf,
           linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_density(alpha = 0.5, fill = "grey55", color = NA) +
  xlim(x_min, x_max) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.title.x     = element_blank()) +
  labs(y = "Densité")

p2 <- ggplot() +
  annotate("segment", x = rope_lo, xend = rope_lo, y = -Inf, yend = Inf,
           linetype = "dashed", color = "grey40", linewidth = 0.4) +
  annotate("segment", x = rope_hi, xend = rope_hi, y = -Inf, yend = Inf,
           linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_segment(aes(x = hdi_95[1], xend = hdi_95[2], y = 0, yend = 0),
               color = "grey40", linewidth = 0.5) +
  geom_segment(aes(x = hdi_80[1], xend = hdi_80[2], y = 0, yend = 0),
               color = "grey30", linewidth = 2) +
  geom_point(aes(x = median_val, y = 0),
             color = "grey40", fill = "grey90", size = 4, shape = 21, stroke = 0.5) +
  xlim(x_min, x_max) +
  ylim(-1, 1) +
  theme_bw() +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank()) +
  labs(x = "Effet de la condition HDT6 (µm)", y = "")

effet <- (p1 / p2 +
             plot_layout(heights = c(4, 1), axes = "collect_x") +
             plot_annotation(title = "Probabilité de l'effet de la condition HDT6\nsur le diamètre moyen des myotubes"))

plot(effet)

# DAG
dag_model <- dagify(
  di ~ mu_i + sigma,
  mu_i ~ alpha_si + b + cond_i,
  alpha_si ~ mu_alpha + sigma_alpha,
  coords = list(
    x = c(mu_alpha = 1, sigma_alpha = 2, sigma = 6,
          alpha_si = 1.5, b = 4.5, cond_i = 3,
          mu_i = 3, di = 4.5),
    y = c(mu_alpha = 3, sigma_alpha = 3, sigma = 2,
          alpha_si = 2, b = 2, cond_i = 2,
          mu_i = 1, di = 0)
  ),
  labels = c(
    di = "d[i]",
    mu_i = "mu[i]",
    sigma = "sigma",
    alpha_si = "alpha[s[i]]",
    b = "beta",
    cond_i = "cond[i]",
    mu_alpha = "mu[alpha]",
    sigma_alpha = "sigma[alpha]"
  )
)

dag <- ggdag(dag_model) +
  geom_dag_edges(edge_color = "grey40", edge_width = 0.7) +
  geom_dag_node(color = "grey20", size = 16) +
  geom_dag_text(aes(label = label), parse = TRUE, color = "white", size = 3.5) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = "Modèle semi-hiérarchique")

# courbe posterior predictive check
# décommenter d_sim

y_rep <- model_fit$BUGSoutput$sims.list$d_sim
y <- data_list$d 

ppc <- ppc_dens_overlay(y, y_rep[1:100, ]) +
   ggtitle("Posterior Predictive Check") +
   labs(x = "Diamètre des myotubes (µm)",
        y = "Densité") +
   theme_bw() +
   theme(
     panel.grid.minor = element_blank()
   )
 
# ggsave(file="sujet.pdf", plot=mirror_plot, width=4.3, height=3)
# ggsave(file="effet.pdf", plot=effet, width=3.7, height=3.2)
# ggsave(file="dag.pdf", plot=dag, width=6, height=4.26)
# ggsave(file="ppc.pdf", plot=ppc, width=6, height=4.26)

