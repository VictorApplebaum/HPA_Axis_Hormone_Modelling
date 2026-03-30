args <- commandArgs(trailingOnly = TRUE)
FoldNo <- as.numeric(sub("--FoldNo=", "", args[1]))

print(FoldNo)



library(rstan)

LOAD_MODEL <- TRUE
FIRST_RUN <- FALSE


setwd("~/HormoneProject/Final")
load('df_full.RData')

load('test_individuals.RData')
test_inds <- test_individuals[FoldNo,]#select test individuals for this fold
#test_inds <- 6:213#for testing code quickly
train_inds <- (1:213)[!(1:213) %in% test_inds]

#Restrict training data to those individuals not in the test
df$measurements <- df$measurements[train_inds,,]
df$measurements_indicator <- df$measurements_indicator[train_inds,,]
df$mdata <- df$mdata[train_inds,]
df$N_i <- length(train_inds)

  


if (LOAD_MODEL == TRUE){
  load('Model.RData')
} else {
  model <- stan_model('Model.stan')
  save(model,file='Model.RData')
}

if (FIRST_RUN == TRUE){
  load('Fit.RData')
  Inits <- list()
  ex_init <- extract(samps)
  # Hormone-level parameters
  Inits$decay_mean  <- ex_init$decay_mean[1000,]
  Inits$decay_sd    <- ex_init$decay_sd[1000,]
  Inits$init_mean   <- ex_init$init_mean[1000,]
  Inits$init_sd     <- ex_init$init_sd[1000,]
  Inits$sigma_mean  <- ex_init$sigma_mean[1000,]
  Inits$sigma_sd    <- ex_init$sigma_sd[1000,]
  
  Inits$decay_v     <- ex_init$decay_v[1000,,]
  Inits$init_v      <- ex_init$init_v[1000,,]
  Inits$sigma_v     <- ex_init$sigma_v[1000,,]
  
  Inits$decay_unscaled <- ex_init$decay_unscaled[1000,train_inds,]
  Inits$init_unscaled  <- ex_init$init_unscaled[1000,train_inds,]
  Inits$sigma_unscaled <- ex_init$sigma_unscaled[1000,train_inds,]
  
  # Hormone enzyme parameters
  enzymes <- c("IIBetaHSD1", "IIBetaHSD2", "VBetaReductase",
               "VAlphaReductase", "XVIIIHydroxylase", "XIBetaHydroxylase")
  
  for (enz in enzymes) {
    Inits[[paste0(enz, "_v")]]        <- ex_init[[paste0(enz, "_v")]][1000,]
    Inits[[paste0(enz, "_unscaled")]] <- ex_init[[paste0(enz, "_unscaled")]][1000,train_inds]
  }
  
  Inits$IIBetaHSD1_mean <- ex_init$IIBetaHSD1_mean[1000]
  Inits$IIBetaHSD1_sd <- ex_init$IIBetaHSD1_sd[1000]
  Inits$IIBetaHSD2_mean <- ex_init$IIBetaHSD2_mean[1000]
  Inits$IIBetaHSD2_sd <-ex_init$IIBetaHSD2_sd[1000]
  Inits$VBetaReductase_mean <- ex_init$VBetaReductase_mean[1000]
  Inits$VBetaReductase_sd <- ex_init$VBetaReductase_sd[1000]
  Inits$VAlphaReductase_mean <- ex_init$VAlphaReductase_mean[1000]
  Inits$VAlphaReductase_sd <- ex_init$VAlphaReductase_sd[1000]
  Inits$XVIIIHydroxylase_mean <- ex_init$XVIIIHydroxylase_mean[1000]
  Inits$XVIIIHydroxylase_sd <- ex_init$XVIIIHydroxylase_sd[1000]
  Inits$XIBetaHydroxylase_mean <- ex_init$XIBetaHydroxylase_mean[1000]
  Inits$XIBetaHydroxylase_sd <- ex_init$XIBetaHydroxylase_sd[1000]
  
  # Pulse 1 parameters
  Inits$pulse_1length_mean              <- ex_init$pulse_1_length_mean[1000]
  Inits$pulse_1_length_sd                <- ex_init$pulse_1_length_sd[1000]
  Inits$pulse_1_length_v                 <- ex_init$pulse_1_length_v[1000,]
  Inits$pulse_1_length          <- ex_init$pulse_1_length[1000,train_inds]
  
  Inits$pulse_1_Cortisol_effect_mean     <- ex_init$pulse_1_Cortisol_effect_mean[1000]
  Inits$pulse_1_Cortisol_effect_sd       <- ex_init$pulse_1_Cortisol_effect_sd[1000]
  Inits$pulse_1_Cortisol_effect_v        <- ex_init$pulse_1_Cortisol_effect_v[1000,]
  Inits$pulse_1_Cortisol_effect_unscaled <- ex_init$pulse_1_Cortisol_effect_unscaled[1000,train_inds]
  
  Inits$pulse_1_CCS_effect_mean          <- ex_init$pulse_1_CCS_effect_mean[1000]
  Inits$pulse_1_CCS_effect_sd            <- ex_init$pulse_1_CCS_effect_sd[1000]
  Inits$pulse_1_CCS_effect_v             <- ex_init$pulse_1_CCS_effect_v[1000,]
  Inits$pulse_1_CCS_effect_unscaled      <- ex_init$pulse_1_CCS_effect_unscaled[1000,train_inds]
  
  # Pulse 2 parameters
  Inits$pulse_2_length_mean              <- ex_init$pulse_2_length_mean[1000]
  Inits$pulse_2_length_sd                <- ex_init$pulse_2_length_sd[1000]
  Inits$pulse_2_length_v                 <- ex_init$pulse_2_length_v[1000,]
  Inits$pulse_2_length          <- ex_init$pulse_2_length[1000,train_inds]
  
  Inits$pulse_2_init_mean                <- ex_init$pulse_2_init_mean[1000]
  Inits$pulse_2_init_sd                  <- ex_init$pulse_2_init_sd[1000]
  Inits$pulse_2_init_v                   <- ex_init$pulse_2_init_v[1000,]
  Inits$pulse_2_init            <- ex_init$pulse_2_init[1000,train_inds]
  
  Inits$pulse_2_Cortisol_effect_mean     <- ex_init$pulse_2_Cortisol_effect_mean[1000]
  Inits$pulse_2_Cortisol_effect_sd       <- ex_init$pulse_2_Cortisol_effect_sd[1000]
  Inits$pulse_2_Cortisol_effect_v        <- ex_init$pulse_2_Cortisol_effect_v[1000,]
  Inits$pulse_2_Cortisol_effect_unscaled <- ex_init$pulse_2_Cortisol_effect_unscaled[1000,train_inds]
  
  Inits$pulse_2_CCS_effect_mean          <- ex_init$pulse_2_CCS_effect_mean[1000]
  Inits$pulse_2_CCS_effect_sd            <- ex_init$pulse_2_CCS_effect_sd[1000]
  Inits$pulse_2_CCS_effect_v             <- ex_init$pulse_2_CCS_effect_v[1000,]
  Inits$pulse_2_CCS_effect_unscaled      <- ex_init$pulse_2_CCS_effect_unscaled[1000,train_inds]
  
  # Residuals
  Inits$eps_raw <- ex_init$eps_raw[1000,,train_inds,]
  
  # Theta (vector<lower=0>[N_hormones])
  Inits$theta <- ex_init$theta[1000,]
} else {
  load(paste0('FoldSaves/FitFold',FoldNo,'.RData'))
  endpoints <- c(300,100,300,300,100,600,150,50,400,100,50,500,600,600,250,300,600,300,110,1)
  endpoint <- endpoints[FoldNo]
  #load(paste0('FoldSaves/FitFold',5,'.RData'))
  #endpoint <- 300
  Inits <- list()
  ex_init <- extract(samps)
  # Hormone-level parameters
  Inits$decay_mean  <- ex_init$decay_mean[endpoint,]
  Inits$decay_sd    <- ex_init$decay_sd[endpoint,]
  Inits$init_mean   <- ex_init$init_mean[endpoint,]
  Inits$init_sd     <- ex_init$init_sd[endpoint,]
  Inits$sigma_mean  <- ex_init$sigma_mean[endpoint,]
  Inits$sigma_sd    <- ex_init$sigma_sd[endpoint,]
  
  Inits$decay_v     <- ex_init$decay_v[endpoint,,]
  Inits$init_v      <- ex_init$init_v[endpoint,,]
  Inits$sigma_v     <- ex_init$sigma_v[endpoint,,]
  
  Inits$decay_unscaled <- ex_init$decay_unscaled[endpoint,,]
  Inits$init_unscaled  <- ex_init$init_unscaled[endpoint,,]
  Inits$sigma_unscaled <- ex_init$sigma_unscaled[endpoint,,]#+log(3)
  
  # Hormone enzyme parameters
  enzymes <- c("IIBetaHSD1", "IIBetaHSD2", "VBetaReductase",
               "VAlphaReductase", "XVIIIHydroxylase", "XIBetaHydroxylase")
  
  for (enz in enzymes) {
    Inits[[paste0(enz, "_v")]]        <- ex_init[[paste0(enz, "_v")]][endpoint,]
    Inits[[paste0(enz, "_unscaled")]] <- ex_init[[paste0(enz, "_unscaled")]][endpoint,]
  }
  
  Inits$IIBetaHSD1_mean <- ex_init$IIBetaHSD1_mean[endpoint]
  Inits$IIBetaHSD1_sd <- ex_init$IIBetaHSD1_sd[endpoint]
  Inits$IIBetaHSD2_mean <- ex_init$IIBetaHSD2_mean[endpoint]
  Inits$IIBetaHSD2_sd <-ex_init$IIBetaHSD2_sd[endpoint]
  Inits$VBetaReductase_mean <- ex_init$VBetaReductase_mean[endpoint]
  Inits$VBetaReductase_sd <- ex_init$VBetaReductase_sd[endpoint]
  Inits$VAlphaReductase_mean <- ex_init$VAlphaReductase_mean[endpoint]
  Inits$VAlphaReductase_sd <- ex_init$VAlphaReductase_sd[endpoint]
  Inits$XVIIIHydroxylase_mean <- ex_init$XVIIIHydroxylase_mean[endpoint]
  Inits$XVIIIHydroxylase_sd <- ex_init$XVIIIHydroxylase_sd[endpoint]
  Inits$XIBetaHydroxylase_mean <- ex_init$XIBetaHydroxylase_mean[endpoint]
  Inits$XIBetaHydroxylase_sd <- ex_init$XIBetaHydroxylase_sd[endpoint]
  
  # Pulse 1 parameters
  Inits$pulse_1_length_mean              <- ex_init$pulse_1_length_mean[endpoint]
  Inits$pulse_1_length_sd                <- ex_init$pulse_1_length_sd[endpoint]
  Inits$pulse_1_length_v                 <- ex_init$pulse_1_length_v[endpoint,]
  Inits$pulse_1_length          <- ex_init$pulse_1_length[endpoint,]
  
  Inits$pulse_1_Cortisol_effect_mean     <- ex_init$pulse_1_Cortisol_effect_mean[endpoint]
  Inits$pulse_1_Cortisol_effect_sd       <- ex_init$pulse_1_Cortisol_effect_sd[endpoint]
  Inits$pulse_1_Cortisol_effect_v        <- ex_init$pulse_1_Cortisol_effect_v[endpoint,]
  Inits$pulse_1_Cortisol_effect_unscaled <- ex_init$pulse_1_Cortisol_effect_unscaled[endpoint,]
  
  Inits$pulse_1_CCS_effect_mean          <- ex_init$pulse_1_CCS_effect_mean[endpoint] + 1
  Inits$pulse_1_CCS_effect_sd            <- ex_init$pulse_1_CCS_effect_sd[endpoint]
  Inits$pulse_1_CCS_effect_v             <- ex_init$pulse_1_CCS_effect_v[endpoint,]
  Inits$pulse_1_CCS_effect_unscaled      <- ex_init$pulse_1_CCS_effect_unscaled[endpoint,]
  
  # Pulse 2 parameters
  Inits$pulse_2_length_mean              <- ex_init$pulse_2_length_mean[endpoint]
  Inits$pulse_2_length_sd                <- ex_init$pulse_2_length_sd[endpoint]
  Inits$pulse_2_length_v                 <- ex_init$pulse_2_length_v[endpoint,]
  Inits$pulse_2_length          <- ex_init$pulse_2_length[endpoint,]
  
  Inits$pulse_2_init_mean                <- ex_init$pulse_2_init_mean[endpoint]
  Inits$pulse_2_init_sd                  <- ex_init$pulse_2_init_sd[endpoint]
  Inits$pulse_2_init_v                   <- ex_init$pulse_2_init_v[endpoint,]
  Inits$pulse_2_init            <- ex_init$pulse_2_init[endpoint,]
  
  Inits$pulse_2_Cortisol_effect_mean     <- ex_init$pulse_2_Cortisol_effect_mean[endpoint]
  Inits$pulse_2_Cortisol_effect_sd       <- ex_init$pulse_2_Cortisol_effect_sd[endpoint]
  Inits$pulse_2_Cortisol_effect_v        <- ex_init$pulse_2_Cortisol_effect_v[endpoint,]
  Inits$pulse_2_Cortisol_effect_unscaled <- ex_init$pulse_2_Cortisol_effect_unscaled[endpoint,]
  
  Inits$pulse_2_CCS_effect_mean          <- ex_init$pulse_2_CCS_effect_mean[endpoint]
  Inits$pulse_2_CCS_effect_sd            <- ex_init$pulse_2_CCS_effect_sd[endpoint]
  Inits$pulse_2_CCS_effect_v             <- ex_init$pulse_2_CCS_effect_v[endpoint,]
  Inits$pulse_2_CCS_effect_unscaled      <- ex_init$pulse_2_CCS_effect_unscaled[endpoint,]
  
  # Residuals
  Inits$eps_raw <- ex_init$eps_raw[endpoint,,,]
  #Inits$eps_raw = array(0,dim(ex_init$eps_raw[endpoint,,,]))
  
  # Theta (vector<lower=0>[N_hormones])
  Inits$theta <- ex_init$theta[endpoint,]
}
init_f <- function(){
  Inits
}

#opt <- optimizing(model,df,init=init_f,verbose=TRUE,iter = 100000, as_vector=FALSE)
#init_2 <- function(){as.list(opt$par)}
samps <- sampling(model,
                  df,
                  chains=12,
                  cores=12,
                  iter= 200,
                  warmup = 150,
                  thin = 1,
                  init = init_f
)
save(samps,file=paste0('FoldSaves/FitFold',FoldNo,'.RData'))
traceplot(samps)

exfit <- extract(samps)

###############Approximate posterior distribution################################
#extract population parameters from fit
#and take logs of positive only parameters
population_samples <- list()

population_samples$pulse_1_length_mean <- exfit$pulse_1_length_mean
population_samples$pulse_1_length_sd <- log(exfit$pulse_1_length_sd)
population_samples$pulse_1_Cortisol_effect_mean <- exfit$pulse_1_Cortisol_effect_mean
population_samples$pulse_1_Cortisol_effect_sd <- log(exfit$pulse_1_Cortisol_effect_sd)
population_samples$pulse_1_CCS_effect_mean <- exfit$pulse_1_CCS_effect_mean
population_samples$pulse_1_CCS_effect_sd <- log(exfit$pulse_1_CCS_effect_sd)
population_samples$pulse_2_length_mean <- exfit$pulse_2_length_mean
population_samples$pulse_2_length_sd <- log(exfit$pulse_2_length_sd)
population_samples$pulse_2_init_mean <- exfit$pulse_2_init_mean
population_samples$pulse_2_init_sd <- log(exfit$pulse_2_init_sd)
population_samples$pulse_2_Cortisol_effect_mean <- exfit$pulse_2_Cortisol_effect_mean
population_samples$pulse_2_Cortisol_effect_sd <- log(exfit$pulse_2_Cortisol_effect_sd)
population_samples$pulse_2_CCS_effect_mean <- exfit$pulse_2_CCS_effect_mean
population_samples$pulse_2_CCS_effect_sd <- log(exfit$pulse_2_CCS_effect_sd)
population_samples$IIBetaHSD2_mean <- exfit$IIBetaHSD2_mean
population_samples$IIBetaHSD2_sd <- log(exfit$IIBetaHSD2_sd)
population_samples$VBetaReductase_mean <- exfit$VBetaReductase_mean
population_samples$VBetaReductase_sd <- log(exfit$VBetaReductase_sd)
population_samples$VAlphaReductase_mean <- exfit$VAlphaReductase_mean
population_samples$VAlphaReductase_sd <- log(exfit$VAlphaReductase_sd)
population_samples$XVIIIHydroxylase_mean <- exfit$XVIIIHydroxylase_mean
population_samples$XVIIIHydroxylase_sd <- log(exfit$XVIIIHydroxylase_sd)
population_samples$XIBetaHydroxylase_mean <- exfit$XIBetaHydroxylase_mean
population_samples$XIBetaHydroxylase_sd <- log(exfit$XIBetaHydroxylase_sd)

for (hormone in 1:df$N_hormones){
  population_samples[[paste0('decay_mean.',hormone)]] <- exfit$decay_mean[,hormone]
  population_samples[[paste0('decay_sd.',hormone)]] <- log(exfit$decay_sd[,hormone])
  population_samples[[paste0('init_mean.',hormone)]] <- exfit$init_mean[,hormone]
  population_samples[[paste0('init_sd.',hormone)]] <- log(exfit$init_sd[,hormone])
  population_samples[[paste0('sigma_mean.',hormone)]] <- exfit$sigma_mean[,hormone]
  population_samples[[paste0('sigma_sd.',hormone)]] <- log(exfit$sigma_sd[,hormone])
}
for (v in 1:df$N_vars){
  population_samples[[paste0('pulse_1_length_v.',v)]] <- exfit$pulse_1_length_v[,v]
  population_samples[[paste0('pulse_1_Cortisol_effect_v.',v)]] <- exfit$pulse_1_Cortisol_effect_v[,v]
  population_samples[[paste0('pulse_1_CCS_effect_v.',v)]] <- exfit$pulse_1_CCS_effect_v[,v]
  population_samples[[paste0('pulse_2_length_v.',v)]] <- exfit$pulse_2_length_v[,v]
  population_samples[[paste0('pulse_2_init_v.',v)]] <- exfit$pulse_2_init_v[,v]
  population_samples[[paste0('pulse_2_Cortisol_effect_v.',v)]] <- exfit$pulse_2_Cortisol_effect_v[,v]
  population_samples[[paste0('pulse_2_CCS_effect_v.',v)]] <- exfit$pulse_2_CCS_effect_v[,v]
  population_samples[[paste0('IIBetaHSD2_v.',v)]] <- exfit$IIBetaHSD2_v[,v]
  population_samples[[paste0('VBetaReductase_v.',v)]] <- exfit$VBetaReductase_v[,v]
  population_samples[[paste0('VAlphaReductase_v.',v)]] <- exfit$VAlphaReductase_v[,v]
  population_samples[[paste0('XVIIIHydroxylase_v.',v)]] <- exfit$XVIIIHydroxylase_v[,v]
  population_samples[[paste0('XIBetaHydroxylase_v.',v)]] <- exfit$XIBetaHydroxylase_v[,v]
}
for (v in 1:df$N_vars){
  for (hormone in 1:df$N_hormones){
    population_samples[[paste0('decay_v.',hormone,'.',v)]] <- exfit$decay_v[,hormone,v]
    population_samples[[paste0('init_v.',hormone,'.',v)]] <- exfit$init_v[,hormone,v]
    population_samples[[paste0('sigma_v.',hormone,'.',v)]] <- exfit$sigma_v[,hormone,v]
  }
}

population_samples$IIBetaHSD1_mean <- exfit$IIBetaHSD1_mean
population_samples$IIBetaHSD1_sd <- log(exfit$IIBetaHSD1_sd)
for (v in 1:df$N_vars){
  population_samples[[paste0('IIBetaHSD1_v.',v)]] <- exfit$IIBetaHSD1_v[,v]
}

for (h in 1:df$N_hormones){
  population_samples[[paste0('theta.',h)]] <- exfit$theta[,h]
}


#now approximate using multivariate normal
pop_par_mean = colMeans(t(matrix(unlist(population_samples), ncol = length(exfit$pulse_1_CCS_effect_mean), byrow = TRUE)))
pop_par_cov  = cov(t(matrix(unlist(population_samples), ncol = length(exfit$pulse_1_CCS_effect_mean), byrow = TRUE)))

save(pop_par_mean,file=paste0('FoldSaves/pop_par_mean',FoldNo,'.RData'))
save(pop_par_cov,file=paste0('FoldSaves/pop_par_cov',FoldNo,'.RData'))
