args <- commandArgs(trailingOnly = TRUE)
FoldNo <- as.numeric(sub("--FoldNo=", "", args[1]))

print(FoldNo)



library(rstan)

setwd("~/HormoneProject/Final")
load('df_full.RData')

LOAD_MODEL <- TRUE
if (LOAD_MODEL == TRUE){
  load('Prediction.RData')
} else {
  model <- stan_model('Prediction.stan')
  save(model,file='Prediction.RData')
}


load('test_individuals.RData')
test_inds <- test_individuals[FoldNo,]#select test individuals for this fold
#test_inds <- 6:213#for testing code quickly

#Restrict training data to those individuals not in the test
df$measurements <- df$measurements[test_inds,,]
df$measurements_indicator <- df$measurements_indicator[test_inds,,]
df$mdata <- df$mdata[test_inds,]
df$N_i <- length(test_inds)


LOAD_MODEL <- TRUE
if (LOAD_MODEL == TRUE){
  load('Prediction.RData')
} else {
  model <- stan_model('Prediction.stan')
  save(model,file='Prediction.RData')
}


load("~/HormoneProject/Test/ValidationInterpolate/Validation_InterpIndividual2.RData")
ex_init <- extract(samps)
Inits = list()
for (parameter in ls(ex_init)){
  #print(length(dim(ex_init[[parameter]])))
  if (length(dim(ex_init[[parameter]])) == 1){
    Inits[[parameter]] = array(ex_init[[parameter]][80],1)
  } else if (length(dim(ex_init[[parameter]])) == 2){
    Inits[[parameter]] = ex_init[[parameter]][80,]
  } else  if (length(dim(ex_init[[parameter]])) == 3){
    Inits[[parameter]] = ex_init[[parameter]][80,,]
  } else  if (length(dim(ex_init[[parameter]])) == 4){
    Inits[[parameter]] = ex_init[[parameter]][80,,,]
  }
  
  x <- Inits[[parameter]]
  
  # Case 1: Scalar with no dim → convert to length-1 vector
  if (is.null(dim(x)) && length(x) == 1) {
    Inits[[parameter]] <- array(x,1)
  }
  
  # Case 2: 1D vector with length > 1 → convert to 1 x N matrix
  else if (is.null(dim(x)) && length(x) == 7) {
    Inits[[parameter]] <- matrix(x, nrow = 1)
  }
  
  #Inits$pop_parameters = array(Inits$pop_parameters,279)  # Case 3: matrix 75×7 → array 75×1×7
  else if (!is.null(dim(x)) && all(dim(x) == c(75, 7))) {
    Inits[[parameter]] <- array(x, dim = c(75, 1, 7))
  }
}
Inits$eps_raw = array(0,c(df$T_max,1,df$N_hormones))
init_f <- function(){Inits}

load("interpolation_gaps.RData")
load("prediction_cutoffs.RData")
load(paste0('FoldSaves/pop_par_mean',FoldNo,'.RData'))
load(paste0('FoldSaves/pop_par_cov',FoldNo,'.RData'))
for (i in 1:df$N_i){
  df_i <- list(measurements = df$measurements[i,,,drop=FALSE],
                      measurements_indicator = df$measurements_indicator[i,,,drop=FALSE],
                      mdata = df$mdata[i,],
                      N_i = 1,
                      T_max = df$T_max,
                      N_vars = df$N_vars,
                      N_hormones = df$N_hormones,
                      N_pop_pars = length(pop_par_mean),
                      pop_par_mean = pop_par_mean,
                      pop_par_cov = chol(pop_par_cov))
  df_i_missingh <- df_i
  df_i_missingh$measurements_indicator[,,2:7] <- array(0,dim(df_i_missingh$measurements_indicator[,,2:7]))
  df_i_interp <- df_i
  df_i_interp$measurements_indicator[,interpolation_gaps[FoldNo,i,1]:interpolation_gaps[FoldNo,i,2],] <- array(0,dim(df_i_interp$measurements_indicator[,interpolation_gaps[FoldNo,i,1]:interpolation_gaps[FoldNo,i,2],,drop=FALSE]))
  df_i_predict <- df_i
  df_i_predict$measurements_indicator[,prediction_cutoffs[FoldNo,i]:df$T_max,] <- array(0,dim(df_i_predict$measurements_indicator[,prediction_cutoffs[FoldNo,i]:df$T_max,,drop=FALSE]))

  samps <- sampling(model,
                    data=df_i_missingh,
                    chains=2,
                    cores=2,
                    iter=1000,
                    warmup = 500,
                    thin = 1,
                    init=init_f
  )
  save(samps,file=paste0('FoldSaves/MissingHF',FoldNo,'I',i,'.RData'))
  samps <- sampling(model,
                    data=df_i_predict,
                    chains=2,
                    cores=2,
                    iter=1000,
                    warmup = 500,
                    thin = 1,
                    init=init_f
  )
  save(samps,file=paste0('FoldSaves/PredictF',FoldNo,'I',i,'.RData'))
  samps <- sampling(model,
                    data=df_i_interp,
                    chains=2,
                    cores=2,
                    iter=1000,
                    warmup = 500,
                    thin = 1,
                    init=init_f
  )
  save(samps,file=paste0('FoldSaves/InterpF',FoldNo,'I',i,'.RData'))
}

