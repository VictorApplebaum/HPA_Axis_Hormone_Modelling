#library(nimble)
library(dplyr)
library(rstan)
library(ggplot2)
library(stargazer)
library(LaplacesDemon)

FitFromScratch <- TRUE



###############Load Data###################################

setwd("~/HormoneProject/Final")
data <- read.csv(file = "md_data_controls.csv")
mdata1 <- read.csv(file = "form_912_controls.csv")
mdata2 <- read.csv(file = "form_914_controls.csv")

#remove one individual with missing values
data <- data[data$MasterID != 465,]
mdata1 <- mdata1[mdata1$MasterID != 465,]
mdata2 <- mdata2[mdata2$MasterID != 465,]

N_used = 213
data$PID <- as.integer(factor(data$MasterID))#Re-number IDs from 1 to N_used
lookup <- unique(data[, c("MasterID", "PID")])#matrix for converting between PID and MasterID
data = data[data$PID %in% unique(data$PID)[1:N_used],]#select only relevent individuals
mdata1 = mdata1[lookup$PID[match(mdata1$'MasterID', lookup$MasterID)] %in% unique(data$PID)[1:N_used],]
mdata2 = mdata2[lookup$PID[match(mdata2$'MasterID', lookup$MasterID)] %in% unique(data$PID)[1:N_used],]
mdata1$PID <- lookup$PID[match(mdata1$'MasterID', lookup$MasterID)]#match mdata PIDs
mdata2$PID <- lookup$PID[match(mdata1$'MasterID', lookup$MasterID)]

#combine mdata1 and mdata2
mdata_collected <- cbind(mdata1,mdata2)

#now remove all uneccesary info from mdata
md_of_interest = c('GenderId','Age','ULT_SMOKER','SYSBP','DIABP','BMI')
mdata <- mdata_collected[ ,md_of_interest]

normalisers <- list('Gender_mean' = mean(mdata$GenderId),
                    'Age_mean' = mean(mdata$Age),
                    'Age_sd' = sd(mdata$Age),
                    'Smoke_mean' = mean(mdata$ULT_SMOKER),
                    'Sys_mean' = mean(mdata$SYSBP,na.rm=TRUE),
                    'Sys_sd' = sd(mdata$SYSBP,na.rm=TRUE),
                    'Dia_mean' = mean(mdata$DIABP,na.rm=TRUE),
                    'Dia_sd' = sd(mdata$DIABP,na.rm=TRUE),
                    'Bmi_mean' = mean(mdata$BMI),
                    'Bmi_sd' = sd(mdata$BMI)
                    )

mdata$GenderId = mdata$GenderId-1#normalisers$Gender_mean
mdata$Age = (mdata$Age -normalisers$Age_mean)/normalisers$Age_sd
mdata$ULT_SMOKER = mdata$ULT_SMOKER-1
mdata$SYSBP = (mdata$SYSBP -normalisers$Sys_mean)/normalisers$Sys_sd
mdata$DIABP = (mdata$DIABP -normalisers$Dia_mean)/normalisers$Dia_sd
mdata$BMI = (mdata$BMI -normalisers$Bmi_mean)/normalisers$Bmi_sd


#########################Dealing with trajectories#########################
#align all time points
group_20min <- as.numeric(
  floor(as.numeric( as.POSIXct(
    data$SampleTime,
    format = "%d/%m/%Y %H:%M",
    tz = "UTC"
  )) / (20 * 60))
)
group_20min <- group_20min-min(group_20min) + 1
data$time_index <- group_20min

#save earliest time
times <- as.POSIXct(data$SampleTime, format = "%d/%m/%Y %H:%M")
earliest_time <- min(times)
save(earliest_time,file='earliest_time.RData')

#now feature scale everything to be in zero-one range more or less
lin_rescale <- function(x){
  (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
}

hormones = c('Cortisol',
             'Cortisone',
             'THF',
             'aTHF',
             'X18OHF',
             'CCS',
             'Aldo'
)
N_hormones <- length(hormones)
measurements <- array(NA,dim=c(N_used,max(data$time_index),N_hormones))

for (point in 1:dim(data)[1]){
  measurements[data$PID[point],data$time_index[point],1] <- data$Cortisol[point]
  measurements[data$PID[point],data$time_index[point],2] <- data$Cortisone[point]
  measurements[data$PID[point],data$time_index[point],3] <- data$THF[point]
  measurements[data$PID[point],data$time_index[point],4] <- data$aTHF[point]
  measurements[data$PID[point],data$time_index[point],5] <- data$X18OHF[point]
  measurements[data$PID[point],data$time_index[point],6] <- data$CCS[point]
  measurements[data$PID[point],data$time_index[point],7] <- data$Aldo[point]
}

measurements_narm <- measurements
measurements_indicator <- array(1,dim=dim(measurements))
for (individual in 1:N_used){
  for (time in 1:dim(measurements_indicator)[2]){
    for (hormone in 1:N_hormones){
      if (is.na(measurements_narm[individual,time,hormone])){
        measurements_narm[individual,time,hormone] <- 0
        measurements_indicator[individual,time,hormone] <- 0
      }
    }
  }
}


#reparameterise
normalisers2 = array(NA,N_hormones)
for (hormone in 1:N_hormones){
  normalisers2[hormone] <- max(measurements_narm[,,hormone],na.rm=TRUE)
  measurements_narm[,,hormone] = measurements_narm[,,hormone]/normalisers2[hormone]
}

df <- list(measurements = measurements_narm,
           measurements_indicator = measurements_indicator,
           mdata = mdata,
           N_hormones = 7,
           N_i = N_used,
           T_max = dim(measurements)[2],
           N_vars = 6)
save(df,file='df_full.RData')



##############Randomizing Test and Training################
set.seed(1)
Folds <- 20
Fold_size <- round(N_used/5)
test_individuals <- array(NA,c(Folds,Fold_size))
for (fold in 1:Folds){
  test_individuals[fold,] <- sample(1:N_used,Fold_size)
}
save(test_individuals,file='test_individuals.RData')


############Cutoffs for testing#########
set.seed(1)
prediction_cutoffs <- array(NA,c(Folds,Fold_size))
for (fold in 1:Folds){
  row_has_one <- rowSums(df$measurements_indicator[1,,] == 1) > 0
  first_row <- which(row_has_one)[1]
  prediction_cutoffs[fold,] <- sample(first_row:df$T_max,Fold_size,replace=TRUE)
}
save(prediction_cutoffs,file='prediction_cutoffs.RData')

interpolation_gaps <- array(NA,c(Folds,Fold_size,2))
for (fold in 1:Folds){
  for (individual in 1:Fold_size){
    
    row_has_one <- rowSums(df$measurements_indicator[1,,] == 1) > 0
    first_row <- which(row_has_one)[1]
    last_row  <- tail(which(row_has_one), 1)
    interpolation_gaps[fold,individual,] = sort(sample(first_row:last_row,2))
  }
}
save(interpolation_gaps,file='interpolation_gaps.RData')


#######Setup Stan models##########
model <- stan_model('Model.stan')
save(model,file='Model.RData')
model <- stan_model('Prediction.stan')
save(model,file='Prediction.RData')

