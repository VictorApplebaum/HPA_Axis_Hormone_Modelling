library(rstan)

i <- 6
FoldNo <- 1
task <- 'Interp'

#in paper/chapter:
#i fold task
#------------
#21 5 Predict
#6 1 Interp
#9 11 MissingH


setwd("~/HormoneProject/Final")
load('df_full.RData')


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


load("interpolation_gaps.RData")
load("prediction_cutoffs.RData")
load(paste0('FoldSaves/pop_par_mean',FoldNo,'.RData'))
load(paste0('FoldSaves/pop_par_cov',FoldNo,'.RData'))

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


  load(paste0('FoldSaves/',task,'F',FoldNo,'I',i,'.RData'))
  N_quantiles = 199
  exfit <- extract(samps)
  mid = array(NA,c(dim(exfit$Measurement_guess)[1],df$T_max,df$N_i,df$N_hormones))
  quant = array(NA,c(N_quantiles,df$T_max,df$N_i,df$N_hormones))
    exfit <- extract(samps)
    for (hormone in 1:df$N_hormones){
      for (time in 1:df$T_max){
        for (iter in 1:dim(exfit$Measurement_guess)[1]){
          mid[iter,time,i,hormone] <- exfit$Measurement_guess[iter,time,1,hormone]
          #mid[iter,time,i,hormone] <- exp(log(exfit$Measurement_guess[iter,time,1,hormone])+rnorm(1,0,exfit$sigma[iter,1,hormone])) * normalisers2[hormone]
        }
        quant[,time,i,hormone] <- quantile(mid[,time,i,hormone],probs=(1:N_quantiles)/(N_quantiles+1),na.rm=TRUE)
      }
    }
  # Create a POSIXct time vector for the full series
  load('earliest_time.RData')
  make_time <- function(n_points) {
    earliest_time  + (0:(n_points-1)) * 20*60
  }
  
  
  hormone_labels <- c("F","E","THF","aTHF","18OHF","CCS","Aldo")
  # Select dataset based on task
  df_plot <- switch(task,
                    Predict  = df_i_predict,
                    MissingH = df_i_missingh,
                    Interp   = df_i_interp
  )
  
  plot_list <- lapply(1:7, function(hormone) {
    
    n_points <- length(quant[195,,i,hormone])
    time_vec <- make_time(n_points)
    
    pred_ind  <- df_plot$measurements_indicator[1,,hormone]
    pred_meas <- df_plot$measurements[1,,hormone]
    
    # Build observed series
    blue_y <- rep(NA, n_points)
    blue_y[pred_ind == 1] <- pred_meas[pred_ind == 1]
    
    red_y <- rep(NA, n_points)
    cond <- pred_ind == 0 & pred_meas != 0
    red_y[cond] <- pred_meas[cond]
    
    # Determine x-axis limits
    valid_idx <- which(!is.na(blue_y) | !is.na(red_y))
    
    if(length(valid_idx) > 0){
      xmin <- time_vec[min(valid_idx)]
      xmax <- time_vec[max(valid_idx)]
    } else {
      xmin <- min(time_vec)
      xmax <- max(time_vec)
    }
    
    # Plot
    ggplot() +
      geom_ribbon(
        aes(x = time_vec,
            ymin = quant[5,,i,hormone],
            ymax = quant[195,,i,hormone]),
        fill = "lightgreen"
      ) +
      geom_line(
        aes(x = time_vec, y = quant[100,,i,hormone]),
        color = "forestgreen",
        linewidth = 2
      ) +
      geom_line(aes(x = time_vec, y = blue_y),
                color = "blue", linewidth = 2) +
      geom_line(aes(x = time_vec, y = red_y),
                color = "red", linewidth = 2) +
      scale_x_datetime(limits = c(xmin, xmax), date_labels = "%H:%M") +
      xlab("Time of Day") +
      ylab("Hormone Concentration (nmol/L)") +
      ggtitle(hormone_labels[hormone]) +
      theme_minimal() +
      theme(text = element_text(size = 16))
  })
  
  library(gridExtra)
  
  full_plot <- grid.arrange(grobs = plot_list, ncol = 4)
  
  print(full_plot)
  ggsave(paste0(task,'.jpg'),
         plot=full_plot,
         width = 17,
         height = 10)  
  