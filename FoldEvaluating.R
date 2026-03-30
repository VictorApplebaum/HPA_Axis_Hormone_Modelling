library(parallel)
setwd("~/HormoneProject/Final")
N_folds = 20

#calculates calibration curves given a fold number and a task ('Predict','Interp', or 'MissingH')
calibration_curve_interp <- function(FoldNo,task){
  library(rstan)
  load('df_full.RData')
  load('test_individuals.RData')
  test_inds <- test_individuals[FoldNo,]#select test individuals for this fold
  #test_inds <- 6:213#for testing code quickly
  
  #Restrict training data to those individuals not in the test
  df$measurements <- df$measurements[test_inds,,]
  df$measurements_indicator <- df$measurements_indicator[test_inds,,]
  df$mdata <- df$mdata[test_inds,]
  df$N_i <- length(test_inds)
  
  N_quantiles = 199
  mid = array(NA,c(1000,df$T_max,df$N_i,df$N_hormones))
  quant = array(NA,c(N_quantiles,df$T_max,df$N_i,df$N_hormones))
  for (individual in 1:df$N_i){
    print(paste0('Calculating quantiles on individual ', individual, 'for task ', task))
    load(paste0('FoldSaves/',task,'F',FoldNo,'I',individual,'.RData'))
    exfit <- extract(samps)
    for (hormone in 1:df$N_hormones){
      for (time in 1:df$T_max){
        for (iter in 1:dim(exfit$Measurement_guess)[1]){
          mid[iter,time,individual,hormone] <- exfit$Measurement_guess[iter,time,1,hormone]
        }
        quant[,time,individual,hormone] <- quantile(mid[,time,individual,hormone],probs=(1:N_quantiles)/(N_quantiles+1),na.rm=TRUE)
      }
    }
  }
  
  load("interpolation_gaps.RData")
  load("prediction_cutoffs.RData")
  
  covered_points <- array(0,c(df$N_hormones+1,99))
  total_points <- array(0,c(df$N_hormones+1,99))
  for (quantile in 1:99){
    for (time in 1:df$T_max){
      for (individual in 1:df$N_i){
        for (hormone in 1:7){
          
          upper <- quant[100 + quantile, time, individual, hormone]
          lower <- quant[100 - quantile, time, individual, hormone]
          
          if (task == 'MissingH'){
            if (hormone %in% 2:7 & df$measurements_indicator[individual,time,hormone] == 1){
              if (df$measurements[individual,time,hormone] < upper & df$measurements[individual,time,hormone] > lower){
                covered_points[hormone,quantile] <- covered_points[hormone,quantile] + 1
              }
              total_points[hormone,quantile] <- total_points[hormone,quantile] + 1
            }
          }
          
          if (task == 'Predict'){
            if (df$measurements_indicator[individual,time,hormone] == 1 & time >= prediction_cutoffs[FoldNo,individual]){
              if (df$measurements[individual,time,hormone] < upper & df$measurements[individual,time,hormone] > lower){
                covered_points[hormone,quantile] <- covered_points[hormone,quantile] + 1
              }
              total_points[hormone,quantile] <- total_points[hormone,quantile] + 1
            }
          }

          if (task == 'Interp'){
            if (df$measurements_indicator[individual,time,hormone] == 1 & time >= interpolation_gaps[FoldNo,individual,1] & time <= interpolation_gaps[FoldNo,individual,2]){
              if (df$measurements[individual,time,hormone] < upper & df$measurements[individual,time,hormone] > lower){
                covered_points[hormone,quantile] <- covered_points[hormone,quantile] + 1
              }
              total_points[hormone,quantile] <- total_points[hormone,quantile] + 1
            }
          }
          
          
        }
        covered_points[8,quantile] <- sum(covered_points[1:7,quantile])
        total_points[8,quantile] <- sum(total_points[1:7,quantile])
      }
    }
  }
  coverage <- list(cov_p = covered_points,tot_p = total_points)
  save(coverage,file=paste0('FoldSaves/Coverage',task,'F',FoldNo,'.RData'))
}



#Now run in parallel
cl <- makeCluster(60)

tasks <- c("Predict", "Interp", "MissingH")
FoldNos <- 1:N_folds
job_grid <- expand.grid(FoldNo = FoldNos, task = tasks)

clusterExport(cl, c("calibration_curve_interp", "job_grid"))

results <- parLapply(cl, 1:nrow(job_grid), function(i) {
  calibration_curve_interp(
    FoldNo = job_grid$FoldNo[i],
    task   = job_grid$task[i]
  )
})

stopCluster(cl)


#Bring results together
task = 'Predict'
cov_p_all <- array(NA,c(8,99,N_folds))
tot_p_all <- array(NA,c(8,99,N_folds))
for (fold in 1:N_folds){
  load(paste0('FoldSaves/Coverage',task,'F',fold,'.RData'))
  cov_p_all[,,fold] <- coverage$cov_p
  tot_p_all[,,fold] <- coverage$tot_p
}
cov_ave <- array(NA,c(8,99))
for (quantile in 1:99){
  for (hormone in 1:8){
    cov_ave[hormone,quantile] = mean(cov_p_all[hormone,quantile,]/tot_p_all[hormone,quantile,])
  }
}

  
# x-axis
x <- (1:99) / 99

hormone_names <- c("F", "E", "THF", "aTHF", 
                   "18OHF", "CCS", "Aldo","Total")

# Build dataframe
df <- data.frame(
  x = rep(x, each = 8),
  coverage = as.vector(cov_ave),
  hormone = factor(rep(hormone_names, times = length(x)),
                   levels = hormone_names)
)

# Diagonal reference line
line_df <- data.frame(
  x = c(0, 1),
  y = c(0, 1)
)

ggplot(df, aes(x = x, y = coverage, color = hormone)) +
  geom_line(size = 3) +
  geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE, size = 3) +
  xlab("Quantile Size") +
  ylab("Proportion of Test Points Contained") +
  scale_color_discrete(name = "Hormone") +
  theme_minimal(base_size = 20)

ggsave("PredictCalibration.png", width = 10, height = 10)

SAD_predict <- array(NA,8)
for (h in 1:8){
  SAD_predict[h] <- sum(abs(cov_ave[h,] - x))
}


task = 'Interp'
cov_p_all <- array(NA,c(8,99,N_folds))
tot_p_all <- array(NA,c(8,99,N_folds))
for (fold in 1:N_folds){
  load(paste0('FoldSaves/Coverage',task,'F',fold,'.RData'))
  cov_p_all[,,fold] <- coverage$cov_p
  tot_p_all[,,fold] <- coverage$tot_p
}
cov_ave <- array(NA,c(8,99))
for (quantile in 1:99){
  for (hormone in 1:8){
    cov_ave[hormone,quantile] = mean(cov_p_all[hormone,quantile,]/tot_p_all[hormone,quantile,])
  }
}


# x-axis
x <- (1:99) / 99

hormone_names <- c("F", "E", "THF", "aTHF", 
                   "18OHF", "CCS", "Aldo","Total")

# Build dataframe
df <- data.frame(
  x = rep(x, each = 8),
  coverage = as.vector(cov_ave),
  hormone = factor(rep(hormone_names, times = length(x)),
                   levels = hormone_names)
)

# Diagonal reference line
line_df <- data.frame(
  x = c(0, 1),
  y = c(0, 1)
)

ggplot(df, aes(x = x, y = coverage, color = hormone)) +
  geom_line(size = 3) +
  geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE, size = 3) +
  xlab("Quantile Size") +
  ylab("Proportion of Test Points Contained") +
  scale_color_discrete(name = "Hormone") +
  theme_minimal(base_size = 20)
ggsave('InterpCalibration.png',
       width = 10,
       height = 10)

SAD_interp <- array(NA,8)
for (h in 1:8){
  SAD_interp[h] <- sum(abs(cov_ave[h,] - x))
}


task = 'MissingH'
cov_p_all <- array(NA,c(8,99,N_folds))
tot_p_all <- array(NA,c(8,99,N_folds))
for (fold in 1:N_folds){
  load(paste0('FoldSaves/Coverage',task,'F',fold,'.RData'))
  cov_p_all[,,fold] <- coverage$cov_p
  tot_p_all[,,fold] <- coverage$tot_p
}
cov_ave <- array(NA,c(8,99))
for (quantile in 1:99){
  for (hormone in 1:8){
    cov_ave[hormone,quantile] = mean(cov_p_all[hormone,quantile,]/tot_p_all[hormone,quantile,])
  }
}


# x-axis
x <- (1:99) / 99

hormone_names <- c("F", "E", "THF", "aTHF", 
                   "18OHF", "CCS", "Aldo","Total")

# Build dataframe
df <- data.frame(
  x = rep(x, each = 8),
  coverage = as.vector(cov_ave),
  hormone = factor(rep(hormone_names, times = length(x)),
                   levels = hormone_names)
)

# Diagonal reference line
line_df <- data.frame(
  x = c(0, 1),
  y = c(0, 1)
)

ggplot(df, aes(x = x, y = coverage, color = hormone)) +
  geom_line(size = 3) +
  geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE, size = 3) +
  xlab("Quantile Size") +
  ylab("Proportion of Test Points Contained") +
  scale_color_discrete(name = "Hormone") +
  theme_minimal(base_size = 20)
ggsave('MissingHCalibration.png',
       width = 10,
       height = 10)

SAD_missing <- array(NA,8)
for (h in 1:8){
  SAD_missing[h] <- sum(abs(cov_ave[h,] - x))
}


SAD_predict
SAD_interp
SAD_missing
