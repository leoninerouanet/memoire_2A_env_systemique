library(readxl)
library(ggplot2)
library(dplyr)
library(R2jags)
library(tidyr)
library(stringr)
library(HDInterval)
library(bayesplot)

setwd("/home/romain/Documents/lele wambs")
output_dir <- file.path(getwd(), "wambs_outputs_proteines")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

# ================================
# DATA (à partir de ton pipeline)
# ================================
data <- read_excel("Quantification_WB.xlsx")

data <- data %>%
  mutate(Signal_moyen = rowMeans(across(c(n1, n2)), na.rm = TRUE))

data_wide <- data %>%
  select(Sujet, Condition, Protéine, Signal_moyen) %>%
  pivot_wider(names_from = Protéine, values_from = Signal_moyen)

data_fc <- data_wide %>%
  mutate(across(
    -c(Sujet, Condition),
    ~ .x / mean(.x[Condition == "BDC-6"], na.rm = TRUE)
  ))

data_fc$condition_id <- ifelse(data_fc$Condition == "HDT6", 1, 0)
data_fc$sujet_id <- as.numeric(as.factor(data_fc$Sujet))

proteines <- setdiff(colnames(data_fc),
                     c("Sujet","Condition","condition_id","sujet_id"))

# ================================
# PREPARE DATA PAR PROT
# ================================
prepare_data_proteines <- function(data_fc, prot) {
  list(
    y = data_fc[[prot]],
    condition = data_fc$condition_id,
    subject = data_fc$sujet_id,
    N = nrow(data_fc),
    N_sujets = max(data_fc$sujet_id)
  )
}

# ================================
# MODEL
# ================================
run_model <- function(data_list, priors, name, n_iter=10000) {
  
  cat("\nRunning:", name, "\n")
  
  model_string <- paste0("
  model {
    for (i in 1:N) {
      y[i] ~ dlnorm(mu[i], tau)
      mu[i] <- alpha[subject[i]] + beta * condition[i]
    }
    
    for (s in 1:N_sujets) {
      alpha[s] ~ dnorm(mu_alpha_log, tau_alpha)
    }
    
    beta ~ dnorm(", priors$beta_mean, ",", priors$beta_prec, ")
    mu_alpha_log <- 0
    
    sigma ~ dnorm(0,", priors$sigma_prec, ") T(0,)
    tau <- 1/pow(sigma,2)
    
    sigma_alpha ~ dnorm(0,", priors$sigma_alpha_prec, ") T(0,)
    tau_alpha <- 1/pow(sigma_alpha,2)
  }")
  
  jags(
    data=data_list,
    parameters.to.save=c("beta","sigma","sigma_alpha","mu_alpha_log"),
    model.file=textConnection(model_string),
    n.chains=4,
    n.iter=n_iter,
    n.burnin=2000
  )
}

# ================================
# PRIORS
# ================================
define_priors <- function() {
  
  base <- list(
    beta_mean = 0,
    beta_prec = 1/(0.55^2),
    sigma_prec = 1,
    sigma_alpha_prec = 1
  )
  
  noninf <- list(
    beta = modifyList(base, list(beta_mean = 0, beta_prec = 1/(2^2)))
  )
  
  sensitivity <- list(
    beta_up_1   = modifyList(base, list(beta_mean = log(2))),
    beta_up_2   = modifyList(base, list(beta_mean = log(3))),
    beta_down_1 = modifyList(base, list(beta_mean = log(0.5))),
    beta_down_2 = modifyList(base, list(beta_mean = log(1/3)))
  )
  
  list(base=base, noninf=noninf, sensitivity=sensitivity)
}

# ================================
# EXTRACTION BETA
# ================================
extract_beta <- function(fit, model_name, prot) {
  
  draws <- fit$BUGSoutput$sims.matrix[,"beta"]
  
  data.frame(
    proteine = prot,
    model = model_name,
    mean = mean(draws),
    sd = sd(draws),
    q2.5 = quantile(draws,0.025),
    q97.5 = quantile(draws,0.975)
  )
}

# ================================
# DEVIATION GENERIQUE
# ================================
compute_deviation_beta <- function(base, other, label, prot, param) {
  
  b <- mean(base$BUGSoutput$sims.matrix[, param])
  o <- mean(other$BUGSoutput$sims.matrix[, param])
  
  dev <- ((b - o) / b) * 100
  
  data.frame(
    proteine   = prot,
    param      = param,
    base       = b,
    other      = o,
    deviation  = dev,
    comparison = label
  )
}

# ================================
# PDF BETA
# ================================
save_plots_pdf_beta <- function(fit, name, prot) {
  
  pdf(file.path(output_dir, paste0(name, "_", prot, ".pdf")), width=11, height=8.5)
  
  params <- c("beta","sigma","sigma_alpha")
  
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
# MAIN LOOP
# ================================
priors <- define_priors()

results_all <- data.frame()
deviation_all <- data.frame()

for (prot in proteines) {
  
  cat("\n====================\n")
  cat("Protéine :", prot, "\n")
  cat("====================\n")
  
  data_list <- prepare_data_proteines(data_fc, prot)
  
  # BASE
  base_model <- run_model(data_list, priors$base, paste0("base_",prot))
  save_plots_pdf_beta(base_model,"base",prot)
  
  results_all <- rbind(results_all,
                       extract_beta(base_model,"base",prot))
  
  # NON INFORMATIVE
  for (n in names(priors$noninf)) {
    
    fit <- run_model(data_list, priors$noninf[[n]],
                     paste0("noninf_",n,"_",prot))
    
    save_plots_pdf_beta(fit,paste0("noninf_",n),prot)
    
    results_all <- rbind(results_all,
                         extract_beta(fit,n,prot))
    
    for (param in c("beta", "sigma", "sigma_alpha")) {
      deviation_all <- rbind(deviation_all,
                             compute_deviation_beta(base_model, fit,
                                                     paste0("noninf_", n),
                                                     prot, param))
    }
  }
  
  # SENSITIVITY
  for (s in names(priors$sensitivity)) {
    
    fit <- run_model(data_list, priors$sensitivity[[s]],
                     paste0(s,"_",prot))
    
    save_plots_pdf_beta(fit,s,prot)
    
    results_all <- rbind(results_all,
                         extract_beta(fit,s,prot))
    
    deviation_all <- rbind(deviation_all,
                           compute_deviation_beta(base_model,fit,s,prot, "beta"))
  }
  
  # DOUBLE ITER
  double_model <- run_model(data_list, priors$base,
                            paste0("double_", prot),
                            n_iter = 20000)
  
  for (param in c("beta", "sigma", "sigma_alpha")) {
    deviation_all <- rbind(deviation_all,
                           compute_deviation_beta(base_model,
                                                   double_model,
                                                   "double_iter",
                                                   prot,
                                                   param))
  }
}

deviation_all$latex_formula <- sprintf(
  "$[(%.3f - %.3f)/%.3f] \\times 100 = %.0f\\%%$",
  deviation_all$base,
  deviation_all$other,
  deviation_all$base,
  deviation_all$deviation
)
# ================================
# EXPORT CSV
# ================================
write.csv(results_all,
          file.path(output_dir,"beta_results_all_proteins.csv"),
          row.names=FALSE)

write.csv(deviation_all,
          file.path(output_dir,"beta_deviation_all_proteins.csv"),
          row.names=FALSE)
