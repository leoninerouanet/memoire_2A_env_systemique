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

setwd("/home/romain/Documents/lele wambs")
output_dir <- file.path(getwd(), "wambs_outputs_myotubes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

# ================================
# DATA
# ================================
prepare_data <- function() {
  data_n1 <- read_excel("/home/romain/Documents/lele wambs/mesures_myotubes_n1_final.xlsx") %>%
    select(Moyenne_fibre, nom_original) %>%
    drop_na() %>%
    mutate(
      condition_raw = str_extract(nom_original, "^[^_]+"),
      sujet = str_extract(nom_original, "(?<=_)[^_]+"),
      diam_um = Moyenne_fibre * 100 / 58130.5
    )
  
  data_n2 <- read_excel("/home/romain/Documents/lele wambs/mesures_myotubes_n2_final.xlsx") %>%
    select(Moyenne_fibre, nom_original) %>%
    drop_na() %>%
    mutate(
      condition_raw = str_extract(nom_original, "^[^_]+"),
      sujet = str_extract(nom_original, "(?<=_)[^_]+"),
      diam_um = Moyenne_fibre * 100 / 219.83
    )
  
  data_clean <- bind_rows(data_n1, data_n2) %>%
    filter(!condition_raw %in% c("CMdiff","CMS")) %>%
    mutate(
      condition_id = ifelse(condition_raw=="HDT6",1,0),
      sujet_id = as.numeric(as.factor(sujet))
    )
  
  list(
    d = data_clean$diam_um,
    condition = data_clean$condition_id,
    subject = data_clean$sujet_id,
    N = nrow(data_clean),
    N_sujets = max(data_clean$sujet_id)
  )
}

# ================================
# MODEL
# ================================
run_model <- function(data_list, priors, name, n_iter=10000) {
  
  cat("\nRunning:", name, "\n")
  
  sigma_prior <- if (priors$sigma_dist == "halfcauchy") {
    paste0("sigma ~ dt(0, ", 1/(priors$sigma_scale^2), ", 1) T(0,)")
  } else {
    paste0("sigma ~ dnorm(0, ", priors$sigma_prec, ") T(0,)")
  }
  
  sigma_alpha_prior <- if (priors$sigma_alpha_dist == "halfcauchy") {
    paste0("sigma_alpha ~ dt(0, ", 1/(priors$sigma_alpha_scale^2), ", 1) T(0,)")
  } else {
    paste0("sigma_alpha ~ dnorm(0, ", priors$sigma_alpha_prec, ") T(0,)")
  }
  
  model_string <- paste0("
  model {
    for (i in 1:N) {
      d[i] ~ dlnorm(mu[i], tau)
      mu[i] <- alpha[subject[i]] + beta * condition[i]
    }
    
    for (s in 1:N_sujets) {
      alpha[s] ~ dnorm(mu_alpha_log, tau_alpha)
    }
    
    beta ~ dnorm(", priors$beta_mean, ",", priors$beta_prec, ")
    mu_alpha_log ~ dnorm(", priors$mu_mean, ",", priors$mu_prec, ")
    
    sigma ~ dnorm(0,", priors$sigma_prec, ") T(0,)
    tau <- 1/pow(sigma,2)
    
    sigma_alpha ~ dnorm(0,", priors$sigma_alpha_prec, ") T(0,)
    tau_alpha <- 1/pow(sigma_alpha,2)
  }")
  
  fit <- jags(
    data=data_list,
    parameters.to.save=c("beta","sigma","sigma_alpha","mu_alpha_log"),
    model.file=textConnection(model_string),
    n.chains=4,
    n.iter=n_iter,
    n.burnin=2000
  )
  
  return(fit)
}

# ================================
# PDF
# ================================
save_plots_pdf <- function(fit, name) {
  
  pdf(file.path(output_dir, paste0(name, ".pdf")), width=11, height=8.5)
  
  params <- c("beta","sigma","sigma_alpha","mu_alpha_log")
  
  for (p in params) {
    samples <- fit$BUGSoutput$sims.matrix[,p]
    
    traceplot(fit, varname = p, axes=FALSE, ask = FALSE)
    hist(samples, main=p, xlab="", ylab="")
    acf(samples, main=p)
    plot(density(samples), main=p, xlab="", ylab="")
  }
  
  dev.off()
}

# ================================
# EXTRACTION COMPLETE
# ================================
extract_full <- function(fit, model_name) {
  
  draws <- fit$BUGSoutput$sims.matrix
  
  res <- data.frame(
    model = model_name,
    param = colnames(draws),
    mean = apply(draws,2,mean),
    sd = apply(draws,2,sd),
    q2.5 = apply(draws,2,quantile,0.025),
    q97.5 = apply(draws,2,quantile,0.975)
  )
  
  return(res)
}

# ================================
# DEVIATION CALCULATION
# ================================
compute_deviation <- function(base, other, label) {
  
  params <- c("beta","sigma","sigma_alpha","mu_alpha_log")
  
  res <- data.frame()
  
  for (p in params) {
    b <- mean(base$BUGSoutput$sims.matrix[,p])
    o <- mean(other$BUGSoutput$sims.matrix[,p])
    
    dev <- ((b - o)/b)*100 # En pourcentage !! Pas besoin de faire x100, déjà fait
    
    res <- rbind(res, data.frame(
      param=p,
      base=b,
      other=o,
      deviation=dev,
      comparison=label
    ))
  }
  
  return(res)
}

# ================================
# DOUBLER ITERATIONS
# ================================

run_double_iteration <- function(data_list, priors, base_fit, name) {
  
  fit_double <- run_model(data_list, priors, name, n_iter = 20000)
  
  dev <- compute_deviation(base_fit, fit_double, paste0(name, "_double_iter"))
  
  return(list(fit = fit_double, deviation = dev))
}

# ================================
# PRIORS 
# ================================
define_priors <- function() {
  
  # ====================
  # BASE MODEL
  # ====================
  base <- list(
    beta_mean = log(1.1),
    beta_prec = 1/(0.275^2),
    
    mu_mean = 3.4,
    mu_prec = 4,
    
    sigma_dist = "halfnorm",
    sigma_prec = 1,
    
    sigma_alpha_dist = "halfnorm",
    sigma_alpha_prec = 1
  )
  
  # ====================
  # NON INFORMATIF
  # ====================
  noninf <- list(
    
    beta = modifyList(base, list(
      beta_mean = 0,
      beta_prec = 1/(5^2)
    )),
    
    mu = modifyList(base, list(
      mu_mean = 0,
      mu_prec = 1/(5^2)
    )),
    
    sigma = modifyList(base, list(
      sigma_dist = "halfcauchy",
      sigma_scale = 2.5
    )),
    
    sigma_alpha = modifyList(base, list(
      sigma_alpha_dist = "halfcauchy",
      sigma_alpha_scale = 2.5
    ))
  )
  
  # ====================
  # SENSITIVITY PAR PARAM
  # ====================
  sensitivity <- list()
  
  # ---- BETA ----
  beta_vals <- list(
    log(1.1),
    log(1.2),
    log(0.9),
    log(0.8)
  )
  
  for (i in seq_along(beta_vals)) {
    sensitivity[[paste0("beta_", i)]] <- modifyList(base, list(
      beta_mean = beta_vals[[i]]
    ))
  }
  
  # ---- MU_ALPHA_LOG ----
  mu_vals <- list(
    log(40),
    log(50),
    log(20),
    log(10)
  )
  
  for (i in seq_along(mu_vals)) {
    sensitivity[[paste0("mu_", i)]] <- modifyList(base, list(
      mu_mean = mu_vals[[i]],
      mu_prec = 1/(0.5^2)
    ))
  }
  
  # ---- SIGMA ----
  sigma_vals <- list(
    log(2),
    log(3)
  )
  
  for (i in seq_along(sigma_vals)) {
    sensitivity[[paste0("sigma_", i)]] <- modifyList(base, list(
      sigma_dist = "halfnorm",
      sigma_meanlog = sigma_vals[[i]],
      sigma_prec = 1
    ))
  }
  
  # ---- SIGMA_ALPHA ----
  sigma_alpha_vals <- list(
    log(1.5),
    log(2)
  )
  
  for (i in seq_along(sigma_alpha_vals)) {
    sensitivity[[paste0("sigma_alpha_", i)]] <- modifyList(base, list(
      sigma_alpha_dist = "halfnorm",
      sigma_alpha_meanlog = sigma_alpha_vals[[i]],
      sigma_alpha_prec = 1
    ))
  }
  
  return(list(
    base = base,
    noninf = noninf,
    sensitivity = sensitivity
  ))
}
# ================================
# MAIN 
# ================================
data_list <- prepare_data()
priors <- define_priors()

results_all <- data.frame()
deviation_all <- data.frame()

# BASE MODEL
base_model <- run_model(data_list, priors$base, "base")
save_plots_pdf(base_model, "base_model")

results_all <- rbind(results_all, extract_full(base_model,"base"))

# NON INFORMATIVE (1x4)
for (n in names(priors$noninf)) {
  
  fit <- run_model(data_list, priors$noninf[[n]], paste0("noninf_",n))
  save_plots_pdf(fit, paste0("noninf_",n))
  
  results_all <- rbind(results_all, extract_full(fit,n))
  deviation_all <- rbind(deviation_all, compute_deviation(base_model, fit, paste0("noninf_",n)))
}

# SENSITIVITY
for (s in names(priors$sensitivity)) {
  
  fit <- run_model(data_list, priors$sensitivity[[s]], s)
  save_plots_pdf(fit, s)
  
  results_all <- rbind(results_all, extract_full(fit,s))
  deviation_all <- rbind(deviation_all, compute_deviation(base_model, fit, s))
}

double_model <- run_model(
  data_list,
  priors$base,
  "base_double_iter",
  n_iter = 20000
)

deviation_all <- rbind(
  deviation_all,
  compute_deviation(base_model, double_model, "base_double_iter")
)
deviation_all$latex_formula <- sprintf(
  "$[(%.3f - %.3f)/%.3f] \\times 100 = %.0f\\%%$",
  deviation_all$base,
  deviation_all$other,
  deviation_all$base,
  deviation_all$deviation
)
# SAVE CSV
write.csv(results_all, file.path(output_dir,"full_results.csv"), row.names=FALSE)
write.csv(deviation_all, file.path(output_dir,"full_deviation.csv"), row.names=FALSE)

