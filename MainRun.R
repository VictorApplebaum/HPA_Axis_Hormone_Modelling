library(rstan)

setwd("~/HormoneProject/Final")
load('df_full.RData')



#Compiling the model sometimes takes a while, so if you have already done this
#you can set LOAD_MODEL to true which will load the last compile
LOAD_MODEL <- TRUE
#Set FIRST_RUN to TRUE to start a new MCMC chain
#If you have already run the MCMC chain and want to continue from where
#you left off, set to FALSE
FIRST_RUN <- FALSE




if (LOAD_MODEL == TRUE){
  load('Model.RData')
} else {
  model <- stan_model('Model.stan')
  save(model,file='Model.RData')
}


if (FIRST_RUN == TRUE){
  Inits <- list()
  # Hormone-level parameters
  Inits$decay_mean  <- c(2.508046, 1.176546, 1.836247, 2.010751, 5.360323, 3.695106, 3.894436)
  Inits$decay_sd    <- c(0.03922010, 0.61520881, 0.41256277, 0.30857612, 0.22413928, 0.18277102, 0.03955081)
  Inits$init_mean   <- c(-3.657463, -3.445528, -4.144719, -3.929321, -2.773490, -4.299015, -3.222718)
  Inits$init_sd     <- c(1.7328474, 1.7673241, 1.6099945, 1.1651299, 0.7949594, 1.5117977, 0.7657641)
  Inits$sigma_mean  <- c(-2.074480, -2.279375, -2.415676, -2.484447, -2.063951, -1.947232, -2.037080)
  Inits$sigma_sd    <- c(0.2663305, 0.1940560, 0.2314435, 0.2635719, 0.2187760, 0.3258356, 0.2605311)
  
  Inits$decay_v     <- array(0, c(df$N_hormones, df$N_vars))
  Inits$init_v      <- array(0, c(df$N_hormones, df$N_vars))
  Inits$sigma_v     <- array(0, c(df$N_hormones, df$N_vars))
  
  Inits$decay_unscaled <- array(0, c(df$N_i, df$N_hormones))
  Inits$init_unscaled  <- array(0, c(df$N_i, df$N_hormones))
  Inits$sigma_unscaled <- array(0, c(df$N_i, df$N_hormones))
  
  # Hormone enzyme parameters
  enzymes <- c("IIBetaHSD1", "IIBetaHSD2", "VBetaReductase",
               "VAlphaReductase", "XVIIIHydroxylase", "XIBetaHydroxylase")
  
  for (enz in enzymes) {
    Inits[[paste0(enz, "_v")]]        <- rep(0, df$N_vars)
    Inits[[paste0(enz, "_unscaled")]] <- rep(0, df$N_i)
  }
  
  Inits$IIBetaHSD1_mean <- -4.236933
  Inits$IIBetaHSD1_sd <- 0.1411491
  Inits$IIBetaHSD2_mean <- -0.9739528
  Inits$IIBetaHSD2_sd <- 0.5012365
  Inits$VBetaReductase_mean <- -1.608593
  Inits$VBetaReductase_sd <- 0.5931526
  Inits$VAlphaReductase_mean <- -2.000691
  Inits$VAlphaReductase_sd <- 0.787518
  Inits$XVIIIHydroxylase_mean <- -4.827795
  Inits$XVIIIHydroxylase_sd <- 0.5456815
  Inits$XIBetaHydroxylase_mean <- -4.236933
  Inits$XIBetaHydroxylase_sd <- 0.1411491
  # Pulse 1 parameters
  Inits$pulse_1_length_mean              <- 15
  Inits$pulse_1_length_sd                <- 4
  Inits$pulse_1_length_v                 <- rep(0, df$N_vars)
  Inits$pulse_1_length          <- rep(15, df$N_i)
  
  Inits$pulse_1_Cortisol_effect_mean     <- -4.969329
  Inits$pulse_1_Cortisol_effect_sd       <- 1.122794
  Inits$pulse_1_Cortisol_effect_v        <- rep(0, df$N_vars)
  Inits$pulse_1_Cortisol_effect_unscaled <- rep(-4.969329, df$N_i)
  
  Inits$pulse_1_CCS_effect_mean          <- -6.105902
  Inits$pulse_1_CCS_effect_sd            <- 0.9885756
  Inits$pulse_1_CCS_effect_v             <- rep(0, df$N_vars)
  Inits$pulse_1_CCS_effect_unscaled      <- rep(-6.105902, df$N_i)
  
  # Pulse 2 parameters
  Inits$pulse_2_length_mean              <- 15
  Inits$pulse_2_length_sd                <- 10
  Inits$pulse_2_length_v                 <- rep(0, df$N_vars)
  Inits$pulse_2_length          <- rep(15, df$N_i)
  
  Inits$pulse_2_init_mean                <- 50
  Inits$pulse_2_init_sd                  <- 10
  Inits$pulse_2_init_v                   <- rep(0, df$N_vars)
  Inits$pulse_2_init            <- rep(50, df$N_i)
  
  Inits$pulse_2_Cortisol_effect_mean     <- -3.485445
  Inits$pulse_2_Cortisol_effect_sd       <- 0.5805732
  Inits$pulse_2_Cortisol_effect_v        <- rep(0, df$N_vars)
  Inits$pulse_2_Cortisol_effect_unscaled <- rep(-3.485445, df$N_i)
  
  Inits$pulse_2_CCS_effect_mean          <- -4.389441
  Inits$pulse_2_CCS_effect_sd            <- 0.6711251
  Inits$pulse_2_CCS_effect_v             <- rep(0, df$N_vars)
  Inits$pulse_2_CCS_effect_unscaled      <- rep(-4.389441, df$N_i)
  
  # Residuals
  Inits$eps_raw <- array(0, c(df$T_max, df$N_i, df$N_hormones))
  
  # Theta (vector<lower=0>[N_hormones])
  Inits$theta <- c(-0.09831517, -0.14980614, -0.27792770, -0.33416476, -0.15964856, -0.13836513, -0.23130715)
  Inits$eps_raw = array(0,dim=c(dim(df$measurements)[2],dim(df$measurements)[1],dim(df$measurements)[3]))
} else {
  Inits <- list()
  load("Fit.RData")
  ex_init <- extract(samps)
  for (parameter in ls(ex_init)){
    #print(length(dim(ex_init[[parameter]])))
    if (length(dim(ex_init[[parameter]])) == 1){
      Inits[[parameter]] = ex_init[[parameter]][60]
    } else if (length(dim(ex_init[[parameter]])) == 2){
      Inits[[parameter]] = ex_init[[parameter]][60,]
    } else  if (length(dim(ex_init[[parameter]])) == 3){
      Inits[[parameter]] = ex_init[[parameter]][60,,]
    } else  if (length(dim(ex_init[[parameter]])) == 4){
      Inits[[parameter]] = ex_init[[parameter]][60,,,]
    }
  }
}
init_f <- function(){
  Inits
}

#opt <- optimizing(model,df,init=init_f,verbose=TRUE,iter = 100000, as_vector=FALSE)
#init_2 <- function(){as.list(opt$par)}
samps <- sampling(model,
                  df,
                  chains=6,
                  cores=6,
                  iter=400,
                  warmup = 300,
                  thin = 1,
                  init = init_f
)
save(samps,file='Fit.RData')

traceplot(samps,'decay[1, 1]')
traceplot(samps,'init[1, 1]')
traceplot(samps,'sigma[1, 1]')


traceplot(samps,'IIBetaHSD2[3]')
traceplot(samps,'IIBetaHSD2_mean')
traceplot(samps,'IIBetaHSD2_sd')
traceplot(samps,'VBetaReductase_mean')
traceplot(samps,'VBetaReductase_sd')
traceplot(samps,'IIBetaHSD1_mean')
traceplot(samps,'IIBetaHSD1_sd')
traceplot(samps,'pulse_1_length[3]')
traceplot(samps,'pulse_2_length[1]')
traceplot(samps,'pulse_2_init_mean')
traceplot(samps,'pulse_2_init_sd')
traceplot(samps,'pulse_2_init[2]')
traceplot(samps,'pulse_2_init[4]')
traceplot(samps,'pulse_1_CCS_effect[1]')
traceplot(samps,'pulse_1_Cortisol_effect[2]')
traceplot(samps,'pulse_1_length_mean')
traceplot(samps,'pulse_2_CCS_effect_sd')
traceplot(samps,'XVIIIHydroxylase_mean')
traceplot(samps,'XVIIIHydroxylase_sd')
traceplot(samps,'XVIIIHydroxylase[4]')

traceplot(samps,'init_mean')
traceplot(samps,'init_sd')
traceplot(samps,'init[3,3]')

traceplot(samps,'decay_mean')
traceplot(samps,'decay_sd')
traceplot(samps,'decay[3,3]')
traceplot(samps,'decay[1,1]')

traceplot(samps,'pulse_2_init[2]')
traceplot(samps,'pulse_2_init_mean')
traceplot(samps,'pulse_2_init_sd')

traceplot(samps,'sigma[2,2]')
traceplot(samps,'sigma_mean')
traceplot(samps,'sigma_sd')
traceplot(samps,'pulse_1_length_v')
traceplot(samps,'pulse_1_CCS_effect_v')
traceplot(samps,'pulse_1_Cortisol_effect_v')
traceplot(samps,'pulse_2_CCS_effect_v')
traceplot(samps,'pulse_2_Cortisol_effect_v')
traceplot(samps,'IIBetaHSD2_v')
traceplot(samps,'VBetaReductase_v')
traceplot(samps,'VAlphaReductase_v')
traceplot(samps,'XVIIIHydroxylase_v')
traceplot(samps,'XIBetaHydroxylase_v')

traceplot(samps,'decay_v[1, 3]')
traceplot(samps,'decay_v[2, 3]')
traceplot(samps,'decay_v[3, 3]')
traceplot(samps,'decay_v[4, 3]')
traceplot(samps,'decay_v[5, 3]')
traceplot(samps,'decay_v[4, 3]')
traceplot(samps,'decay_v[5, 3]')
traceplot(samps,'decay_v[1, 1]')

traceplot(samps,'pulse_2_init_v')
traceplot(samps,'pulse_2_length_v')
traceplot(samps,'sigma_v[1, 1]')

traceplot(samps,'theta')



traceplot(samps,c('IIBetaHSD2_mean','VBetaReductase_mean','VAlphaReductase_mean','XVIIIHydroxylase_mean','XIBetaHydroxylase_mean','IIBetaHSD1_mean'))
traceplot(samps,c('IIBetaHSD2_sd','VBetaReductase_sd','VAlphaReductase_sd','XVIIIHydroxylase_sd','XIBetaHydroxylase_sd','IIBetaHSD1_sd'))

exfit <- extract(samps)

########################Posterior Odds Ratios################################################################
POR <- array(NA,c(13+3*7,df$N_vars))
for (v in 1:df$N_vars){
  POR[,v] <- c(
    sum(exfit$pulse_1_length_v[,v] > 0)/sum(exfit$pulse_1_length_v[,v] < 0),
    sum(exfit$pulse_1_Cortisol_effect_v[,v] > 0)/sum(exfit$pulse_1_Cortisol_effect_v[,v] < 0),
    sum(exfit$pulse_1_CCS_effect_v[,v] > 0)/sum(exfit$pulse_1_CCS_effect_v[,v] < 0),
    sum(exfit$pulse_2_length_v[,v] > 0)/sum(exfit$pulse_2_length_v[,v] < 0),
    sum(exfit$pulse_2_init_v[,v] > 0)/sum(exfit$pulse_2_init_v[,v] < 0),
    sum(exfit$pulse_2_Cortisol_effect_v[,v] > 0)/sum(exfit$pulse_2_Cortisol_effect_v[,v] < 0),
    sum(exfit$pulse_2_CCS_effect_v[,v] > 0)/sum(exfit$pulse_2_CCS_effect_v[,v] < 0),
    sum(exfit$IIBetaHSD1_v[,v] > 0)/sum(exfit$IIBetaHSD1_v[,v] < 0),
    sum(exfit$IIBetaHSD2_v[,v] > 0)/sum(exfit$IIBetaHSD2_v[,v] < 0),
    sum(exfit$VBetaReductase_v[,v] > 0)/sum(exfit$VBetaReductase_v[,v] < 0),
    sum(exfit$VAlphaReductase_v[,v] > 0)/sum(exfit$VAlphaReductase_v[,v] < 0),
    sum(exfit$XVIIIHydroxylase_v[,v] > 0)/sum(exfit$XVIIIHydroxylase_v[,v] < 0),
    sum(exfit$XIBetaHydroxylase_v[,v] > 0)/sum(exfit$XIBetaHydroxylase_v[,v] < 0),
    sum(exfit$decay_v[,1,v] > 0)/sum(exfit$decay_v[,1,v] < 0),
    sum(exfit$decay_v[,2,v] > 0)/sum(exfit$decay_v[,2,v] < 0),
    sum(exfit$decay_v[,3,v] > 0)/sum(exfit$decay_v[,3,v] < 0),
    sum(exfit$decay_v[,4,v] > 0)/sum(exfit$decay_v[,4,v] < 0),
    sum(exfit$decay_v[,5,v] > 0)/sum(exfit$decay_v[,5,v] < 0),
    sum(exfit$decay_v[,6,v] > 0)/sum(exfit$decay_v[,6,v] < 0),
    sum(exfit$decay_v[,7,v] > 0)/sum(exfit$decay_v[,7,v] < 0),
    sum(exfit$init_v[,1,v] > 0)/sum(exfit$init_v[,1,v] < 0),
    sum(exfit$init_v[,2,v] > 0)/sum(exfit$init_v[,2,v] < 0),
    sum(exfit$init_v[,3,v] > 0)/sum(exfit$init_v[,3,v] < 0),
    sum(exfit$init_v[,4,v] > 0)/sum(exfit$init_v[,4,v] < 0),
    sum(exfit$init_v[,5,v] > 0)/sum(exfit$init_v[,5,v] < 0),
    sum(exfit$init_v[,6,v] > 0)/sum(exfit$init_v[,6,v] < 0),
    sum(exfit$init_v[,7,v] > 0)/sum(exfit$init_v[,7,v] < 0),
    sum(exfit$sigma_v[,1,v] > 0)/sum(exfit$sigma_v[,1,v] < 0),
    sum(exfit$sigma_v[,2,v] > 0)/sum(exfit$sigma_v[,2,v] < 0),
    sum(exfit$sigma_v[,3,v] > 0)/sum(exfit$sigma_v[,3,v] < 0),
    sum(exfit$sigma_v[,4,v] > 0)/sum(exfit$sigma_v[,4,v] < 0),
    sum(exfit$sigma_v[,5,v] > 0)/sum(exfit$sigma_v[,5,v] < 0),
    sum(exfit$sigma_v[,6,v] > 0)/sum(exfit$sigma_v[,6,v] < 0),
    sum(exfit$sigma_v[,7,v] > 0)/sum(exfit$sigma_v[,7,v] < 0)
  )
}
POR

library(reshape2)  # For reshaping the data
library(RColorBrewer)  # For colorblind-friendly palettes
library(latex2exp)
library(ggplot2)
library(dplyr)

# Example data (replace with actual data)
set.seed(42)
posterior_odds <- POR

# Define limits
cap_high <- 100
cap_low <- 1 / 100

# Cap values
posterior_odds <- pmin(pmax(posterior_odds, cap_low), cap_high)

# Convert matrix to a data frame for ggplot
df_por <- melt(posterior_odds)
colnames(df_por) <- c("Row", "Column", "Value")

# Apply a log transformation centered at 1
df_por$LogValue <- log(df_por$Value)

# Assign Jeffreys' scale categories
df_por <- df_por %>%
  mutate(Category = case_when(
    LogValue < log(1/100) ~ "Decisive -",
    LogValue < log(1/30)  ~ "V. Strong -",
    LogValue < log(1/10)  ~ "Strong -",
    LogValue < log(1/3)   ~ "Substantial -",
    LogValue < log(3)     ~ "Weak",
    LogValue < log(10)    ~ "Substantial +",
    LogValue < log(30)    ~ "Strong +",
    LogValue < log(100)   ~ "V. Strong +",
    TRUE                  ~ "Decisive +"
  ))

# Choose a colorblind-friendly colormap
color_palette <- rev(brewer.pal(11, "RdBu"))

# Set a threshold to switch text color
df_por <- df_por %>%
  mutate(TextColor = ifelse(abs(LogValue) > log(10), "white", "black"))

# Custom y-axis labels formatted correctly
y_labels <- rev(c(
  "l[1]",
  "alpha^{F1}",
  "alpha^{CCS1}",
  "l[2]",
  "xi",
  "alpha^{F2}",
  "alpha^{CCS2}",
  "alpha^{F3}",
  "alpha^{E}",
  "alpha^{THF}",
  "alpha^{aTHF}",
  "alpha^{\"18OHF\"}",
  "alpha^{Aldo}",
  "beta^{F}",
  "beta^{E}",
  "beta^{THF}",
  "beta^{aTHF}",
  "beta^{\"18OHF\"}",
  "beta^{CCS}",
  "beta^{Aldo}",
  "Y[0]^{F}",
  "Y[0]^{E}",
  "Y[0]^{THF}",
  "Y[0]^{aTHF}",
  "Y[0]^{\"18OHF\"}",
  "Y[0]^{CCS}",
  "Y[0]^{Aldo}",
  "sigma^{F}",
  "sigma^{E}",
  "sigma^{THF}",
  "sigma^{aTHF}",
  "sigma^{\"18OHF\"}",
  "sigma^{CCS}",
  "sigma^{Aldo}"
))


# Custom x-axis labels
x_labels <- c("Gender", "Year of Birth", "Smoking", "Systolic BP", "Diastolic BP", "BMI")

# Convert Row to a factor with explicit levels
df_por$Row <- factor(df_por$Row, levels = rev(unique(df_por$Row)), labels = y_labels)

######################## BOTH HEATMAPS #################################

cap_high <- 100
cap_low  <- 1 / 100

posterior_odds <- pmin(pmax(POR, cap_low), cap_high)

df_por <- reshape2::melt(posterior_odds)
colnames(df_por) <- c("Row", "Column", "Value")

df_por$LogValue <- log(df_por$Value)

# ------------------- Jeffreys categories -------------------
df_por$Category <- dplyr::case_when(
  df_por$LogValue < log(1/100) ~ "Decisive -",
  df_por$LogValue < log(1/30)  ~ "V. Strong -",
  df_por$LogValue < log(1/10)  ~ "Strong -",
  df_por$LogValue < log(1/3)   ~ "Substantial -",
  df_por$LogValue < log(3)     ~ "Weak",
  df_por$LogValue < log(10)    ~ "Substantial +",
  df_por$LogValue < log(30)    ~ "Strong +",
  df_por$LogValue < log(100)   ~ "V. Strong +",
  TRUE                         ~ "Decisive +"
)

# ------------------- Numeric POR labels -------------------
df_por$POR_label <- sprintf("%.2f", df_por$LogValue)

# Text colour
df_por$TextColor <- ifelse(abs(df_por$LogValue) > log(10), "white", "black")

# Axis formatting (reuse yours)
df_por$Row <- factor(df_por$Row,
                     levels = rev(unique(df_por$Row)),
                     labels = y_labels)

color_palette <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

# ------------------- Plot 1: Jeffreys -------------------
p1 <- ggplot(df_por, aes(x = factor(Column, labels = x_labels),
                         y = Row, fill = LogValue)) +
  geom_tile() +
  geom_text(aes(label = Category),
            color = df_por$TextColor, size = 3, fontface = "bold") +
  scale_fill_gradientn(colors = color_palette,
                       limits = c(-log(cap_high), log(cap_high))) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  labs(x = "Predictors", y = "Parameters", fill = "Log Odds") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------- Plot 2: Numeric POR -------------------
p2 <- ggplot(df_por, aes(x = factor(Column, labels = x_labels),
                         y = Row, fill = LogValue)) +
  geom_tile() +
  geom_text(aes(label = POR_label),
            color = df_por$TextColor, size = 3, fontface = "bold") +
  scale_fill_gradientn(colors = color_palette,
                       limits = c(-log(cap_high), log(cap_high))) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  labs(x = "Predictors", y = "Parameters", fill = "Log Odds") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save both
ggsave("HeatPlot_Jeffreys.png", p1, width = 10, height = 8)
ggsave("HeatPlot_NumericPOR.png", p2, width = 10, height = 8)

#############################Credible Intervals######################################
N_quantiles <- 199
pulse_1_length_v <- array(NA,c(N_quantiles,df$N_vars))
pulse_1_Cortisol_effect_v <- array(NA,c(N_quantiles,df$N_vars))
pulse_1_CCS_effect_v <- array(NA,c(N_quantiles,df$N_vars))
pulse_2_length_v <- array(NA,c(N_quantiles,df$N_vars))
pulse_2_init_v <- array(NA,c(N_quantiles,df$N_vars))
pulse_2_Cortisol_effect_v <- array(NA,c(N_quantiles,df$N_vars))
pulse_2_CCS_effect_v <- array(NA,c(N_quantiles,df$N_vars))
IIBetaHSD1_v <- array(NA,c(N_quantiles,df$N_vars))
IIBetaHSD2_v <- array(NA,c(N_quantiles,df$N_vars))
VBetaReductase_v <- array(NA,c(N_quantiles,df$N_vars))
VAlphaReductase_v <- array(NA,c(N_quantiles,df$N_vars))
XVIIIHydroxylase_v <- array(NA,c(N_quantiles,df$N_vars))
XIBetaHydroxylase_v <- array(NA,c(N_quantiles,df$N_vars))
decay_v <- array(NA,c(N_quantiles,df$N_hormones,df$N_vars))
init_v <- array(NA,c(N_quantiles,df$N_hormones,df$N_vars))
sigma_v <- array(NA,c(N_quantiles,df$N_hormones,df$N_vars))
for (v in 1:df$N_vars){
  pulse_1_length_v[,v] <- quantile(exfit$pulse_1_length_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  pulse_1_Cortisol_effect_v[,v] <- quantile(exfit$pulse_1_Cortisol_effect_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  pulse_1_CCS_effect_v[,v] <- quantile(exfit$pulse_1_CCS_effect_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  pulse_2_length_v[,v] <- quantile(exfit$pulse_2_length_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  pulse_2_init_v[,v] <- quantile(exfit$pulse_2_init_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  pulse_2_Cortisol_effect_v[,v] <- quantile(exfit$pulse_2_Cortisol_effect_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  pulse_2_CCS_effect_v[,v] <- quantile(exfit$pulse_2_CCS_effect_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  IIBetaHSD1_v[,v] <- quantile(exfit$IIBetaHSD1_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  IIBetaHSD2_v[,v] <- quantile(exfit$IIBetaHSD2_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  VBetaReductase_v[,v] <- quantile(exfit$VBetaReductase_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  VAlphaReductase_v[,v]  <- quantile(exfit$VAlphaReductase_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  XVIIIHydroxylase_v[,v] <- quantile(exfit$XVIIIHydroxylase_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  XIBetaHydroxylase_v[,v] <- quantile(exfit$XIBetaHydroxylase_v[,v],probs=(1:N_quantiles)/(N_quantiles+1))
  for (hormone in 1:df$N_hormones){
    decay_v[,hormone,v] <- quantile(exfit$decay_v[,hormone,v],probs=(1:N_quantiles)/(N_quantiles+1))
    init_v[,hormone,v] <- quantile(exfit$init_v[,hormone,v],probs=(1:N_quantiles)/(N_quantiles+1))
    sigma_v[,hormone,v] <- quantile(exfit$sigma_v[,hormone,v],probs=(1:N_quantiles)/(N_quantiles+1))
  }
}

v_effects <- array(NA,c(13+3*7,2*df$N_vars))
for (v in 1:df$N_vars){
  upper_outputs = c(pulse_1_length_v[195,v],
                    pulse_1_Cortisol_effect_v[195,v],
                    pulse_1_CCS_effect_v[195,v],
                    pulse_2_length_v[195,v],
                    pulse_2_init_v[195,v],
                    pulse_2_Cortisol_effect_v[195,v],
                    pulse_2_CCS_effect_v[195,v],
                    IIBetaHSD1_v[195,v],
                    IIBetaHSD2_v[195,v],
                    VBetaReductase_v[195,v],
                    VAlphaReductase_v[195,v],
                    XVIIIHydroxylase_v[195,v],
                    XIBetaHydroxylase_v[195,v],
                    decay_v[195,1:df$N_hormones,v],
                    init_v[195,1:df$N_hormones,v],
                    sigma_v[195,1:df$N_hormones,v]
  )
  lower_outputs = c(pulse_1_length_v[5,v],
                    pulse_1_Cortisol_effect_v[5,v],
                    pulse_1_CCS_effect_v[5,v],
                    pulse_2_length_v[5,v],
                    pulse_2_init_v[5,v],
                    pulse_2_Cortisol_effect_v[5,v],
                    pulse_2_CCS_effect_v[5,v],
                    IIBetaHSD1_v[5,v],
                    IIBetaHSD2_v[5,v],
                    VBetaReductase_v[5,v],
                    VAlphaReductase_v[5,v],
                    XVIIIHydroxylase_v[5,v],
                    XIBetaHydroxylase_v[5,v],
                    decay_v[5,1:df$N_hormones,v],
                    init_v[5,1:df$N_hormones,v],
                    sigma_v[5,1:df$N_hormones,v]
  )
  v_effects[,2*v] = upper_outputs
  v_effects[,2*v-1] = lower_outputs
}

# Convert to dataframe
df_ci <- reshape2::melt(v_effects)
colnames(df_ci) <- c("Row", "Column", "Value")

# Separate lower/upper
df_ci$Type <- ifelse(df_ci$Column %% 2 == 1, "lower", "upper")
df_ci$Var  <- ceiling(df_ci$Column / 2)

# Reshape
df_ci <- reshape2::melt(v_effects)
colnames(df_ci) <- c("Row", "Column", "Value")

# Identify lower/upper and variable index
df_ci$Type <- ifelse(df_ci$Column %% 2 == 1, "lower", "upper")
df_ci$Var  <- ceiling(df_ci$Column / 2)

# Wide format: one row per (Row, Var)
df_ci <- reshape2::dcast(df_ci, Row + Var ~ Type, value.var = "Value")

# Format labels: (lower, upper)
df_ci$CI_label <- sprintf("(%.2f, %.2f)", df_ci$lower, df_ci$upper)

# Axis labels
df_ci$Row <- factor(df_ci$Row,
                    levels = rev(unique(df_ci$Row)),
                    labels = y_labels)
df_ci$Var <- factor(df_ci$Var, labels = x_labels)

p_ci_text <- ggplot(df_ci, aes(x = Var, y = Row)) +
  geom_tile(fill = "white", color = "grey80") +
  geom_text(aes(label = CI_label), size = 3) +
  scale_y_discrete(labels = function(x) parse(text = x)) +  # <-- THIS FIX
  theme_minimal() +
  labs(
    x = "Predictors",
    y = "Parameters"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("CredibleIntervals.png", p_ci_text, width = 10, height = 8)
