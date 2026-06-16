
set.seed(2025)

pacotes <- c('readxl','numDeriv','rstan','bayesplot', 'Metrics','ggplot2','rstudioapi')

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = TRUE)
  }
  sapply(pacotes, require, character = TRUE) 
} else {
  sapply(pacotes, require, character = TRUE) 
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read_excel("kinetic.xlsx")
data$T <- data$T + 273.15

n <- nrow(data)
Tref <- 1 / (sum(1 / unique(data$T)) / length(unique(data$T)))

stan_data <- list(
  n = n,
  yCO2 = data$yCO2,
  yH2 = data$yH2,
  T = data$T,
  rCH4_obs = data$rCH4,
  Tref = Tref
)

# Stan Model
stan_model2_code <- "
data {
  int<lower=0> n;
  vector[n] yCO2;
  vector[n] yH2;
  vector[n] T;
  vector[n] rCH4_obs;
  real Tref;
}

parameters {
  real lnk5;
  real Ea5_RT;
  real lnk6;
  real Ea6_RT;
  real lnk10;
  real Ea10_RT;
  real lnK1;
  real dH1_RT;
  real lnK2;
  real dH2_RT;
  real lnK3;
  real dH3_RT;
  real lnK4;
  real dH4_RT;
}

model {
  lnk5 ~ normal(-6.70, 3);
  Ea5_RT ~ normal(40.63, 20);
  lnk6 ~ normal(22.92, 10);
  Ea6_RT ~ normal(58.78, 25);
  lnk10 ~ normal(68.31, 30);
  Ea10_RT ~ normal(56.33, 25);
  lnK1 ~ normal(91.60, 20);
  dH1_RT ~ normal(-7.16, 3);
  lnK2 ~ normal(-81.71, 40);
  dH2_RT ~ normal(-24.08, 10);
  lnK3 ~ normal(42.35, 20);
  dH3_RT ~ normal(-45.96, 25);
  lnK4 ~ normal(42.35, 20);
  dH4_RT ~ normal(-45.96, 25);
  
  for (i in 1:n) {
    real k5 = exp(lnk5 + Ea5_RT * (T[i] - Tref) / T[i]);
    real k6 = exp(lnk6 + Ea6_RT * (T[i] - Tref) / T[i]);
    real k10 = exp(lnk10 + Ea10_RT * (T[i] - Tref) / T[i]);
    real K1 = exp(lnK1 + dH1_RT * (T[i] - Tref) / T[i]);
    real K2 = exp(lnK2 + dH2_RT * (T[i] - Tref) / T[i]);
    real K3 = exp(lnK3 + dH3_RT * (T[i] - Tref) / T[i]);
    real K4 = exp(lnK4 + dH4_RT * (T[i] - Tref) / T[i]);
    
    real den = (1 + sqrt(K1)*sqrt(yH2[i]) + K2*sqrt(K1)*yCO2[i]*sqrt(yH2[i]) +
    sqrt(k10*k5*K4*K2*sqrt(K1))*((1/k10)+(1/k6)+(1/k5)+(K3/(k5*K4)))*sqrt(yCO2[i])*(yH2[i]^0.25))^2;
    
    real rCH4 = sqrt(k10*k5*K4*K2*(K1^1.5)) * sqrt(yCO2[i]) * (yH2[i]^0.75) / den;
    
    rCH4_obs[i] ~ normal(rCH4, 0.000001);  
  }
}
"

writeLines(stan_model2_code, "kinetic_model2.stan")

stan_model <- stan_model(file = "kinetic_model2.stan")

fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4, control = list(max_treedepth = 15, adapt_delta = 0.95))

print(fit)

traceplot(fit, pars = c("lnk5", "Ea5_RT", "lnk6", "Ea6_RT", "lnk10", "Ea10_RT", "lnK1", "dH1_RT", "lnK2", "dH2_RT", "lnK3", "dH3_RT", "lnK4", "dH4_RT"))

pairs(fit, pars = c("lnk5", "Ea5_RT", "lnk6", "Ea6_RT", "lnk10", "Ea10_RT", "lnK1", "dH1_RT", "lnK2", "dH2_RT", "lnK3", "dH3_RT", "lnK4", "dH4_RT"))

mcmc_dens(fit, facet_args = list(ncol = 2), pars = c("lp__")) # need bayesplot

posterior_samples <- extract(fit)
y_pred <- matrix(NA, nrow = length(posterior_samples$lnk5), ncol = n)
for (i in 1:length(posterior_samples$lnk5)) {
  k5 <- exp(posterior_samples$lnk5[i] + posterior_samples$Ea5_RT[i] * (data$T - Tref) / data$T)
  k6 <- exp(posterior_samples$lnk6[i] + posterior_samples$Ea6_RT[i] * (data$T - Tref) / data$T)
  k10 <- exp(posterior_samples$lnk10[i] + posterior_samples$Ea10_RT[i] * (data$T - Tref) / data$T)
  K1 <- exp(posterior_samples$lnK1[i] + posterior_samples$dH1_RT[i] * (data$T - Tref) / data$T)
  K2 <- exp(posterior_samples$lnK2[i] + posterior_samples$dH2_RT[i] * (data$T - Tref) / data$T)
  K3 <- exp(posterior_samples$lnK3[i] + posterior_samples$dH3_RT[i] * (data$T - Tref) / data$T)
  K4 <- exp(posterior_samples$lnK4[i] + posterior_samples$dH4_RT[i] * (data$T - Tref) / data$T)
  den <- (1 + sqrt(K1)*sqrt(data$yH2) + K2*sqrt(K1)*data$yCO2*sqrt(data$yH2) +
            sqrt(k10*k5*K4*K2*sqrt(K1))*((1/k10)+(1/k6)+(1/k5)+(K3/(k5*K4)))*sqrt(data$yCO2)*(data$yH2^0.25))^2
  y_pred[i, ] <- sqrt(k10*k5*K4*K2*(K1^1.5))*sqrt(data$yCO2)*(data$yH2^0.75)/den
}
mean_predictions <- apply(y_pred, 2, mean)
plot(data$rCH4, mean_predictions, xlab = "Observed rCH4", ylab = "Predicted rCH4")
abline(0, 1, col = "red")
upper_bound <- data$rCH4 * 1.05
lower_bound <- data$rCH4 * 0.95
lines(data$rCH4, upper_bound, col = "gray", lty = 2)
lines(data$rCH4, lower_bound, col = "gray", lty = 2)

# statistics

observed <- data$rCH4
predicted <- mean_predictions

RSS <- sum((observed - predicted)^2)
TSS <- sum((observed - mean(observed))^2)
R2 <- 1 - (RSS / TSS)

p <- 14
R2_adj <- 1 - ((1 - R2) * (n - 1) / (n - p - 1))

rmse_value <- rmse(observed, predicted)

cat("R²:", R2, "\n")
cat("Adjusted R²:", R2_adj, "\n")
cat("RMSE:", rmse_value, "\n")
