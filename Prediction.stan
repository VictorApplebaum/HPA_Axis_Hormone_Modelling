data {
  int N_i;
  int N_vars;
  int N_hormones;
  int T_max;
  
  real measurements[N_i,T_max,N_hormones];
  real measurements_indicator[N_i,T_max,N_hormones];
  //vector[N_vars] mdata;
  real mdata[N_i,N_vars];

  
  int N_pop_pars;
  vector[N_pop_pars] pop_par_mean;
  matrix[N_pop_pars,N_pop_pars] pop_par_cov;
}

parameters {
  real<lower=0> pulse_1_length[N_i];
  real pulse_1_Cortisol_effect_unscaled[N_i];
  real pulse_1_CCS_effect_unscaled[N_i];
  
  
  real pulse_2_Cortisol_effect_unscaled[N_i];
  real pulse_2_CCS_effect_unscaled[N_i];
    
  real<lower=0> pulse_2_length[N_i];
  real<lower=0> pulse_2_init[N_i];

  real decay_unscaled[N_i,N_hormones];
  real sigma_unscaled[N_i,N_hormones];
  real init_unscaled[N_i,N_hormones];

  
  //Hormones
  real IIBetaHSD2_unscaled[N_i];
  real VBetaReductase_unscaled[N_i];
  real VAlphaReductase_unscaled[N_i];
  real XVIIIHydroxylase_unscaled[N_i];
  real XIBetaHydroxylase_unscaled[N_i];
  real IIBetaHSD1_unscaled[N_i];
  
  array[N_pop_pars] real pop_parameters;
    real eps_raw[T_max,N_i,N_hormones];


}


model {



  to_vector(pop_parameters) ~ multi_normal_cholesky(pop_par_mean,pop_par_cov);

  //Then assign names to each
  real pulse_1_length_mean;
  real pulse_1_length_sd;
  real pulse_1_Cortisol_effect_mean;
  real pulse_1_Cortisol_effect_sd ;
  real pulse_1_CCS_effect_mean;
  real pulse_1_CCS_effect_sd;
  
  real pulse_2_length_mean;
  real pulse_2_length_sd;
  real pulse_2_init_mean;
  real pulse_2_init_sd;
  real pulse_2_Cortisol_effect_mean;
  real pulse_2_Cortisol_effect_sd;
  real pulse_2_CCS_effect_mean;
  real pulse_2_CCS_effect_sd;
  
  
  real IIBetaHSD2_mean;
  real IIBetaHSD2_sd;
  real VBetaReductase_mean;
  real VBetaReductase_sd;
  real VAlphaReductase_mean;
  real VAlphaReductase_sd;
  real XVIIIHydroxylase_mean;
  real XVIIIHydroxylase_sd;
  real XIBetaHydroxylase_mean;
  real XIBetaHydroxylase_sd;
  real IIBetaHSD1_mean;
  real IIBetaHSD1_sd;
    pulse_1_length_mean = pop_parameters[1];
    pulse_1_length_sd = exp(pop_parameters[2]);
    pulse_1_Cortisol_effect_mean = pop_parameters[3];
    pulse_1_Cortisol_effect_sd = exp(pop_parameters[4]);
    pulse_1_CCS_effect_mean = pop_parameters[5];
    pulse_1_CCS_effect_sd= exp(pop_parameters[6]);
  
    pulse_2_length_mean = pop_parameters[7];
    pulse_2_length_sd = exp(pop_parameters[8]);
    pulse_2_init_mean = pop_parameters[9];
    pulse_2_init_sd = exp(pop_parameters[10]);
    pulse_2_Cortisol_effect_mean = pop_parameters[11];
    pulse_2_Cortisol_effect_sd = exp(pop_parameters[12]);
    pulse_2_CCS_effect_mean = pop_parameters[13];
    pulse_2_CCS_effect_sd = exp(pop_parameters[14]);
  
    IIBetaHSD2_mean = pop_parameters[15];
    IIBetaHSD2_sd = exp(pop_parameters[16]);
    VBetaReductase_mean = pop_parameters[17];
    VBetaReductase_sd = exp(pop_parameters[18]);
    VAlphaReductase_mean = pop_parameters[19];
    VAlphaReductase_sd = exp(pop_parameters[20]);
    XVIIIHydroxylase_mean = pop_parameters[21];
    XVIIIHydroxylase_sd = exp(pop_parameters[22]);
    XIBetaHydroxylase_mean = pop_parameters[23];
    XIBetaHydroxylase_sd = exp(pop_parameters[24]);
    IIBetaHSD1_mean = pop_parameters[265];
    IIBetaHSD1_sd = exp(pop_parameters[266]);
  

  
  array[N_hormones] real decay_mean;
  array[N_hormones] real decay_sd;
  array[N_hormones] real init_mean;
  array[N_hormones] real init_sd;
  array[N_hormones] real sigma_mean;
  array[N_hormones] real sigma_sd;
    for (hormone in 1:N_hormones){
      decay_mean[hormone] = pop_parameters[6*(hormone-1) + 25];
      decay_sd[hormone] = exp(pop_parameters[6*(hormone-1) + 26]);
      init_mean[hormone] = pop_parameters[6*(hormone-1) + 27];
      init_sd[hormone] = exp(pop_parameters[6*(hormone-1) + 28]);
      sigma_mean[hormone] = pop_parameters[6*(hormone-1) + 29];
      sigma_sd[hormone] = exp(pop_parameters[6*(hormone-1) + 30]);

    }

  
  array[N_vars] real pulse_1_length_v;
  array[N_vars] real pulse_1_Cortisol_effect_v;
  array[N_vars] real pulse_1_CCS_effect_v;
  array[N_vars] real pulse_2_length_v;
  array[N_vars] real pulse_2_init_v;
  array[N_vars] real pulse_2_Cortisol_effect_v;
  array[N_vars] real pulse_2_CCS_effect_v;
  array[N_vars] real IIBetaHSD2_v;
  array[N_vars] real VBetaReductase_v;
  array[N_vars] real VAlphaReductase_v;
  array[N_vars] real XVIIIHydroxylase_v;
  array[N_vars] real XIBetaHydroxylase_v;
  array[N_vars] real IIBetaHSD1_v;


    for (v in 1:N_vars){
      pulse_1_length_v[v] = pop_parameters[12*(v-1) + 67];
      pulse_1_Cortisol_effect_v[v] = pop_parameters[12*(v-1) + 68];
      pulse_1_CCS_effect_v[v] = pop_parameters[12*(v-1) + 69];
      pulse_2_length_v[v] = pop_parameters[12*(v-1) + 70];
      pulse_2_init_v[v] = pop_parameters[12*(v-1) + 71];
      pulse_2_Cortisol_effect_v[v] = pop_parameters[12*(v-1) + 72];
      pulse_2_CCS_effect_v[v] = pop_parameters[12*(v-1) + 73];
      IIBetaHSD2_v[v] = pop_parameters[12*(v-1) + 74];
      VBetaReductase_v[v] = pop_parameters[12*(v-1) + 75];
      VAlphaReductase_v[v] = pop_parameters[12*(v-1) + 76];
      XVIIIHydroxylase_v[v] = pop_parameters[12*(v-1) + 77];
      XIBetaHydroxylase_v[v] = pop_parameters[12*(v-1) + 78];
      IIBetaHSD1_v[v] = pop_parameters[v + 266];
    }
  

  array[N_hormones,N_vars] real decay_v;
  array[N_hormones,N_vars] real init_v;
  array[N_hormones,N_vars] real sigma_v;

    for (v in 1:N_vars){
      for (hormone in 1:N_hormones){
        decay_v[hormone,v] = pop_parameters[3*(hormone+(7*(v-1))-1)+139];
        init_v[hormone,v] = pop_parameters[3*(hormone+(7*(v-1))-1)+140];
        sigma_v[hormone,v] = pop_parameters[3*(hormone+(7*(v-1))-1)+141];
      }
    }


  
  //Individual level
  
  for (individual in 1:N_i){
    pulse_1_length[individual] ~ normal(pulse_1_length_mean + dot_product(to_vector(pulse_1_length_v), to_vector(mdata[individual,])) ,pulse_1_length_sd);
    pulse_1_Cortisol_effect_unscaled[individual] ~ normal(pulse_1_Cortisol_effect_mean + dot_product(to_vector(pulse_1_Cortisol_effect_v), to_vector(mdata[individual,])) ,pulse_1_Cortisol_effect_sd);
    pulse_1_CCS_effect_unscaled[individual] ~ normal(pulse_1_CCS_effect_mean + dot_product(to_vector(pulse_1_CCS_effect_v), to_vector(mdata[individual,])) ,pulse_1_CCS_effect_sd);
    
    pulse_2_length[individual] ~ normal(pulse_2_length_mean + dot_product(to_vector(pulse_2_length_v), to_vector(mdata[individual,])),pulse_2_length_sd);
    pulse_2_init[individual] ~ normal(pulse_2_init_mean + dot_product(to_vector(pulse_2_init_v), to_vector(mdata[individual,])) ,pulse_2_init_sd);
    pulse_2_Cortisol_effect_unscaled[individual] ~ normal(pulse_2_Cortisol_effect_mean + dot_product(to_vector(pulse_2_Cortisol_effect_v), to_vector(mdata[individual,])) ,pulse_2_Cortisol_effect_sd);
    pulse_2_CCS_effect_unscaled[individual] ~ normal(pulse_2_CCS_effect_mean + dot_product(to_vector(pulse_2_CCS_effect_v), to_vector(mdata[individual,])) ,pulse_2_CCS_effect_sd);
    for (hormone in 1:N_hormones){
      decay_unscaled[individual,hormone] ~ normal(0,1);//decay_mean[hormone] + dot_product(to_vector(decay_v[hormone,]), to_vector(mdata[individual,])) ,decay_sd[hormone]);
      init_unscaled[individual,hormone] ~ normal(init_mean[hormone] + dot_product(to_vector(init_v[hormone,]), to_vector(mdata[individual,])) ,init_sd[hormone]);
      sigma_unscaled[individual,hormone] ~ normal(sigma_mean[hormone] + dot_product(to_vector(sigma_v[hormone,]), to_vector(mdata[individual,])),sigma_sd[hormone]);
    }


  


    IIBetaHSD2_unscaled[individual] ~ normal(0,1);//IIBetaHSD2_mean + dot_product(to_vector(IIBetaHSD2_v), to_vector(mdata[individual,])),IIBetaHSD2_sd);
    VBetaReductase_unscaled[individual] ~ normal(0,1);//VBetaReductase_mean + dot_product(to_vector(VBetaReductase_v), to_vector(mdata[individual,])),VBetaReductase_sd);
    VAlphaReductase_unscaled[individual] ~ normal(0,1);//VAlphaReductase_mean + dot_product(to_vector(VAlphaReductase_v), to_vector(mdata[individual,])),VAlphaReductase_sd);
    XVIIIHydroxylase_unscaled[individual] ~ normal(0,1);//XVIIIHydroxylase_mean + dot_product(to_vector(XVIIIHydroxylase_v), to_vector(mdata[individual,])),XVIIIHydroxylase_sd);
    XIBetaHydroxylase_unscaled[individual] ~ normal(0,1);//XIBetaHydroxylase_mean + dot_product(to_vector(XIBetaHydroxylase_v), to_vector(mdata[individual,])),XIBetaHydroxylase_sd);
    IIBetaHSD1_unscaled[individual] ~ normal(0,1);//IIBetaHSD1_mean + dot_product(to_vector(IIBetaHSD1_v), to_vector(mdata[individual,])),IIBetaHSD1_sd);
  }
  


  real decay[N_i,N_hormones];
  real init[N_i,N_hormones];
  real sigma[N_i,N_hormones];


  real pulse_1_Cortisol_effect[N_i];
  real pulse_1_CCS_effect[N_i];
  real pulse_2_Cortisol_effect[N_i];
  real pulse_2_CCS_effect[N_i];
  real IIBetaHSD2[N_i];
  real VBetaReductase[N_i];
  real VAlphaReductase[N_i];
  real XVIIIHydroxylase[N_i] ;
  real XIBetaHydroxylase[N_i];
  real IIBetaHSD1[N_i];



  for (i in 1:N_i){
    for (hormone in 1:N_hormones){
      decay[i,hormone] = 1/(1+exp(-(decay_unscaled[i,hormone] * decay_sd[hormone] + decay_mean[hormone] + dot_product(to_vector(decay_v[hormone,]), to_vector(mdata[i,])))));
      init[i,hormone] = exp(init_unscaled[i,hormone]);
      sigma[i,hormone] = exp(sigma_unscaled[i,hormone]);// * sigma_sd[hormone] + sigma_mean[hormone] + dot_product(to_vector(sigma_v[hormone,]), to_vector(mdata[i,])));
    }
    pulse_1_Cortisol_effect[i] = exp(pulse_1_Cortisol_effect_unscaled[i]);
    pulse_1_CCS_effect[i] =  exp(pulse_1_CCS_effect_unscaled[i]);
    pulse_2_Cortisol_effect[i] = exp(pulse_2_Cortisol_effect_unscaled[i]);
    pulse_2_CCS_effect[i] = exp(pulse_2_CCS_effect_unscaled[i]);
    VBetaReductase[i] = exp(VBetaReductase_unscaled[i]* VBetaReductase_sd + VBetaReductase_mean + dot_product(to_vector(VBetaReductase_v), to_vector(mdata[i,])));
    VAlphaReductase[i] = exp(VAlphaReductase_unscaled[i]* VAlphaReductase_sd + VAlphaReductase_mean + dot_product(to_vector(VAlphaReductase_v), to_vector(mdata[i,])));
    XVIIIHydroxylase[i] = exp(XVIIIHydroxylase_unscaled[i]* XVIIIHydroxylase_sd + XVIIIHydroxylase_mean + dot_product(to_vector(XVIIIHydroxylase_v), to_vector(mdata[i,])));
    XIBetaHydroxylase[i] = exp(XIBetaHydroxylase_unscaled[i]* XIBetaHydroxylase_sd + XIBetaHydroxylase_mean + dot_product(to_vector(XIBetaHydroxylase_v), to_vector(mdata[i,])));
    IIBetaHSD1[i] = exp(IIBetaHSD1_unscaled[i]* IIBetaHSD1_sd + IIBetaHSD1_mean + dot_product(to_vector(IIBetaHSD1_v), to_vector(mdata[i,])));

    IIBetaHSD2[i] = exp(IIBetaHSD2_unscaled[i]* IIBetaHSD2_sd + IIBetaHSD2_mean + dot_product(to_vector(IIBetaHSD2_v), to_vector(mdata[i,])));

    
  }
  
  vector[N_hormones] theta = to_vector(pop_parameters[273:279]);

  ///////Model
  real Measurement_guess[T_max,N_i,N_hormones];
  real eps[T_max,N_i,N_hormones];

  ////Chain

  for (individual in 1:N_i){
    vector[N_hormones] prev_points =to_vector(init[individual, ]);
    for (time in 1:T_max){
      Measurement_guess[time,individual,1] = decay[individual,1] *prev_points[1] + IIBetaHSD1[individual] * prev_points[2];
      Measurement_guess[time,individual,6] = decay[individual,6] *prev_points[6];
      //add pulse 1
      if (time < (0 + 1) && time > 0 && time < (pulse_1_length[individual]+0)){
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_1_Cortisol_effect[individual] * (time - 0);
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_1_CCS_effect[individual] * (time - 0);
      } else if (time > 0 && time <= (pulse_1_length[individual]+0)){//other case where we are not at the first pulse point, but the pulse is still happenning
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_1_Cortisol_effect[individual];
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_1_CCS_effect[individual];
      } else if (time <= (pulse_1_length[individual]+0 + 1) && time > (pulse_1_length[individual]+0)){//other case where we are not at the final point
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_1_Cortisol_effect[individual] * (-time + (1+pulse_1_length[individual]+0));
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_1_CCS_effect[individual] * (-time + (1+pulse_1_length[individual]+0));
      }
      
      
      //add pulse 2
      //if we are on the first pulse point: add the increase multiplied by the time the pulse has been over
      //to do this we need time to be between the inittime and the inittime+1, while still before the pulse has ended
      if (time < (pulse_2_init[individual] + 1) && time > pulse_2_init[individual] && time < (pulse_2_length[individual]+pulse_2_init[individual])){
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_2_Cortisol_effect[individual] * (time - pulse_2_init[individual]);
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_2_CCS_effect[individual] * (time - pulse_2_init[individual]);
      } else if (time > pulse_2_init[individual] && time <= (pulse_2_length[individual]+pulse_2_init[individual])){//other case where we are not at the first pulse point, but the pulse is still happenning
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_2_Cortisol_effect[individual];
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_2_CCS_effect[individual];
      } else if (time <= (pulse_2_length[individual]+pulse_2_init[individual] + 1) && time > (pulse_2_length[individual]+pulse_2_init[individual])){//other case where we are not at the final point
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_2_Cortisol_effect[individual] * (-time + (1+pulse_2_length[individual]+pulse_2_init[individual]));
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_2_CCS_effect[individual] * (-time + (1+pulse_2_length[individual]+pulse_2_init[individual]));
      }
      
      for (hormone in {1,6}){
        eps_raw[time,individual,hormone] ~ normal(0,1);
        if (measurements_indicator[individual,time,hormone] == 1){
          eps[time,individual,hormone] = log(measurements[individual,time,hormone] / Measurement_guess[time,individual,hormone]);
          Measurement_guess[time,individual,hormone] = measurements[individual,time,hormone];
          target += normal_lpdf(eps[time,individual,hormone] | 0, (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone])));
        } else {
          eps[time,individual,hormone] = eps_raw[time,individual,hormone] * (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone]));
          Measurement_guess[time,individual,hormone] = Measurement_guess[time,individual,hormone] * exp(eps[time,individual,hormone]);
        }
      }
      
      Measurement_guess[time,individual,2] = decay[individual,2] *prev_points[2] + IIBetaHSD2[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,3] = decay[individual,3] *prev_points[3] + VBetaReductase[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,4] = decay[individual,4] *prev_points[4] + VAlphaReductase[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,5] = decay[individual,5] *prev_points[5] + XVIIIHydroxylase[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,7] = decay[individual,7] *prev_points[7] + XIBetaHydroxylase[individual] * Measurement_guess[time,individual,6];

  
      
      for (hormone in {2,3,4,5,7}){
        eps_raw[time,individual,hormone] ~ normal(0,1);
        if (measurements_indicator[individual,time,hormone] == 1){
          eps[time,individual,hormone] = log(measurements[individual,time,hormone] / Measurement_guess[time,individual,hormone]);
          Measurement_guess[time,individual,hormone] = measurements[individual,time,hormone];
          target += normal_lpdf(eps[time,individual,hormone] | 0, (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone])));
        } else {
          eps[time,individual,hormone] = eps_raw[time,individual,hormone] * (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone]));
          Measurement_guess[time,individual,hormone] = Measurement_guess[time,individual,hormone] * exp(eps[time,individual,hormone]);
        }
      }
      
      
      prev_points = to_vector(Measurement_guess[time, individual, ]);
    }
  }
}
generated quantities {
  


  //Then assign names to each
  real pulse_1_length_mean;
  real pulse_1_length_sd;
  real pulse_1_Cortisol_effect_mean;
  real pulse_1_Cortisol_effect_sd ;
  real pulse_1_CCS_effect_mean;
  real pulse_1_CCS_effect_sd;
  
  real pulse_2_length_mean;
  real pulse_2_length_sd;
  real pulse_2_init_mean;
  real pulse_2_init_sd;
  real pulse_2_Cortisol_effect_mean;
  real pulse_2_Cortisol_effect_sd;
  real pulse_2_CCS_effect_mean;
  real pulse_2_CCS_effect_sd;
  
  
  real IIBetaHSD2_mean;
  real IIBetaHSD2_sd;
  real VBetaReductase_mean;
  real VBetaReductase_sd;
  real VAlphaReductase_mean;
  real VAlphaReductase_sd;
  real XVIIIHydroxylase_mean;
  real XVIIIHydroxylase_sd;
  real XIBetaHydroxylase_mean;
  real XIBetaHydroxylase_sd;
  real IIBetaHSD1_mean;
  real IIBetaHSD1_sd;
    pulse_1_length_mean = pop_parameters[1];
    pulse_1_length_sd = exp(pop_parameters[2]);
    pulse_1_Cortisol_effect_mean = pop_parameters[3];
    pulse_1_Cortisol_effect_sd = exp(pop_parameters[4]);
    pulse_1_CCS_effect_mean = pop_parameters[5];
    pulse_1_CCS_effect_sd= exp(pop_parameters[6]);
  
    pulse_2_length_mean = pop_parameters[7];
    pulse_2_length_sd = exp(pop_parameters[8]);
    pulse_2_init_mean = pop_parameters[9];
    pulse_2_init_sd = exp(pop_parameters[10]);
    pulse_2_Cortisol_effect_mean = pop_parameters[11];
    pulse_2_Cortisol_effect_sd = exp(pop_parameters[12]);
    pulse_2_CCS_effect_mean = pop_parameters[13];
    pulse_2_CCS_effect_sd = exp(pop_parameters[14]);
  
    IIBetaHSD2_mean = pop_parameters[15];
    IIBetaHSD2_sd = exp(pop_parameters[16]);
    VBetaReductase_mean = pop_parameters[17];
    VBetaReductase_sd = exp(pop_parameters[18]);
    VAlphaReductase_mean = pop_parameters[19];
    VAlphaReductase_sd = exp(pop_parameters[20]);
    XVIIIHydroxylase_mean = pop_parameters[21];
    XVIIIHydroxylase_sd = exp(pop_parameters[22]);
    XIBetaHydroxylase_mean = pop_parameters[23];
    XIBetaHydroxylase_sd = exp(pop_parameters[24]);
    IIBetaHSD1_mean = pop_parameters[265];
    IIBetaHSD1_sd = exp(pop_parameters[266]);
  

  
  array[N_hormones] real decay_mean;
  array[N_hormones] real decay_sd;
  array[N_hormones] real init_mean;
  array[N_hormones] real init_sd;
  array[N_hormones] real sigma_mean;
  array[N_hormones] real sigma_sd;


    for (hormone in 1:N_hormones){
      decay_mean[hormone] = pop_parameters[6*(hormone-1) + 25];
      decay_sd[hormone] = exp(pop_parameters[6*(hormone-1) + 26]);
      init_mean[hormone] = pop_parameters[6*(hormone-1) + 27];
      init_sd[hormone] = exp(pop_parameters[6*(hormone-1) + 28]);
      sigma_mean[hormone] = pop_parameters[6*(hormone-1) + 29];
      sigma_sd[hormone] = exp(pop_parameters[6*(hormone-1) + 30]);

    }

  
  array[N_vars] real pulse_1_length_v;
  array[N_vars] real pulse_1_Cortisol_effect_v;
  array[N_vars] real pulse_1_CCS_effect_v;
  array[N_vars] real pulse_2_length_v;
  array[N_vars] real pulse_2_init_v;
  array[N_vars] real pulse_2_Cortisol_effect_v;
  array[N_vars] real pulse_2_CCS_effect_v;
  array[N_vars] real IIBetaHSD2_v;
  array[N_vars] real VBetaReductase_v;
  array[N_vars] real VAlphaReductase_v;
  array[N_vars] real XVIIIHydroxylase_v;
  array[N_vars] real XIBetaHydroxylase_v;
  array[N_vars] real IIBetaHSD1_v;


    for (v in 1:N_vars){
      pulse_1_length_v[v] = pop_parameters[12*(v-1) + 67];
      pulse_1_Cortisol_effect_v[v] = pop_parameters[12*(v-1) + 68];
      pulse_1_CCS_effect_v[v] = pop_parameters[12*(v-1) + 69];
      pulse_2_length_v[v] = pop_parameters[12*(v-1) + 70];
      pulse_2_init_v[v] = pop_parameters[12*(v-1) + 71];
      pulse_2_Cortisol_effect_v[v] = pop_parameters[12*(v-1) + 72];
      pulse_2_CCS_effect_v[v] = pop_parameters[12*(v-1) + 73];
      IIBetaHSD2_v[v] = pop_parameters[12*(v-1) + 74];
      VBetaReductase_v[v] = pop_parameters[12*(v-1) + 75];
      VAlphaReductase_v[v] = pop_parameters[12*(v-1) + 76];
      XVIIIHydroxylase_v[v] = pop_parameters[12*(v-1) + 77];
      XIBetaHydroxylase_v[v] = pop_parameters[12*(v-1) + 78];
      IIBetaHSD1_v[v] = pop_parameters[v + 266];
    }
  

  array[N_hormones,N_vars] real decay_v;
  array[N_hormones,N_vars] real init_v;
  array[N_hormones,N_vars] real sigma_v;

    for (v in 1:N_vars){
      for (hormone in 1:N_hormones){
        decay_v[hormone,v] = pop_parameters[3*(hormone+(7*(v-1))-1)+139];
        init_v[hormone,v] = pop_parameters[3*(hormone+(7*(v-1))-1)+140];
        sigma_v[hormone,v] = pop_parameters[3*(hormone+(7*(v-1))-1)+141];
      }
    }




  real decay[N_i,N_hormones];
  real init[N_i,N_hormones];
  real sigma[N_i,N_hormones];


  real pulse_1_Cortisol_effect[N_i];
  real pulse_1_CCS_effect[N_i];
  real pulse_2_Cortisol_effect[N_i];
  real pulse_2_CCS_effect[N_i];
  real IIBetaHSD2[N_i];
  real VBetaReductase[N_i];
  real VAlphaReductase[N_i];
  real XVIIIHydroxylase[N_i] ;
  real XIBetaHydroxylase[N_i];
  real IIBetaHSD1[N_i];

  for (i in 1:N_i){
    for (hormone in 1:N_hormones){
      decay[i,hormone] = 1/(1+exp(-(decay_unscaled[i,hormone] * decay_sd[hormone] + decay_mean[hormone] + dot_product(to_vector(decay_v[hormone,]), to_vector(mdata[i,])))));
      init[i,hormone] = exp(init_unscaled[i,hormone]);
      sigma[i,hormone] = exp(sigma_unscaled[i,hormone]);// * sigma_sd[hormone] + sigma_mean[hormone] + dot_product(to_vector(sigma_v[hormone,]), to_vector(mdata[i,])));
    }
    pulse_1_Cortisol_effect[i] = exp(pulse_1_Cortisol_effect_unscaled[i]);
    pulse_1_CCS_effect[i] =  exp(pulse_1_CCS_effect_unscaled[i]);
    pulse_2_Cortisol_effect[i] = exp(pulse_2_Cortisol_effect_unscaled[i]);
    pulse_2_CCS_effect[i] = exp(pulse_2_CCS_effect_unscaled[i]);
    VBetaReductase[i] = exp(VBetaReductase_unscaled[i]* VBetaReductase_sd + VBetaReductase_mean + dot_product(to_vector(VBetaReductase_v), to_vector(mdata[i,])));
    VAlphaReductase[i] = exp(VAlphaReductase_unscaled[i]* VAlphaReductase_sd + VAlphaReductase_mean + dot_product(to_vector(VAlphaReductase_v), to_vector(mdata[i,])));
    XVIIIHydroxylase[i] = exp(XVIIIHydroxylase_unscaled[i]* XVIIIHydroxylase_sd + XVIIIHydroxylase_mean + dot_product(to_vector(XVIIIHydroxylase_v), to_vector(mdata[i,])));
    XIBetaHydroxylase[i] = exp(XIBetaHydroxylase_unscaled[i]* XIBetaHydroxylase_sd + XIBetaHydroxylase_mean + dot_product(to_vector(XIBetaHydroxylase_v), to_vector(mdata[i,])));
    IIBetaHSD1[i] = exp(IIBetaHSD1_unscaled[i]* IIBetaHSD1_sd + IIBetaHSD1_mean + dot_product(to_vector(IIBetaHSD1_v), to_vector(mdata[i,])));

    IIBetaHSD2[i] = exp(IIBetaHSD2_unscaled[i]* IIBetaHSD2_sd + IIBetaHSD2_mean + dot_product(to_vector(IIBetaHSD2_v), to_vector(mdata[i,])));

    
  }

  vector[N_hormones] theta = to_vector(pop_parameters[273:279]);

  ///////Model
  real Measurement_guess[T_max,N_i,N_hormones];
  real eps[T_max,N_i,N_hormones];

  ////Chain

  for (individual in 1:N_i){
    vector[N_hormones] prev_points =to_vector(init[individual, ]);
    for (time in 1:T_max){
      Measurement_guess[time,individual,1] = decay[individual,1] *prev_points[1] + IIBetaHSD1[individual] * prev_points[2];
      Measurement_guess[time,individual,6] = decay[individual,6] *prev_points[6];
      //add pulse 1
      if (time < (0 + 1) && time > 0 && time < (pulse_1_length[individual]+0)){
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_1_Cortisol_effect[individual] * (time - 0);
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_1_CCS_effect[individual] * (time - 0);
      } else if (time > 0 && time <= (pulse_1_length[individual]+0)){//other case where we are not at the first pulse point, but the pulse is still happenning
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_1_Cortisol_effect[individual];
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_1_CCS_effect[individual];
      } else if (time <= (pulse_1_length[individual]+0 + 1) && time > (pulse_1_length[individual]+0)){//other case where we are not at the final point
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_1_Cortisol_effect[individual] * (-time + (1+pulse_1_length[individual]+0));
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_1_CCS_effect[individual] * (-time + (1+pulse_1_length[individual]+0));
      }
      
      
      //add pulse 2
      //if we are on the first pulse point: add the increase multiplied by the time the pulse has been over
      //to do this we need time to be between the inittime and the inittime+1, while still before the pulse has ended
      if (time < (pulse_2_init[individual] + 1) && time > pulse_2_init[individual] && time < (pulse_2_length[individual]+pulse_2_init[individual])){
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_2_Cortisol_effect[individual] * (time - pulse_2_init[individual]);
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_2_CCS_effect[individual] * (time - pulse_2_init[individual]);
      } else if (time > pulse_2_init[individual] && time <= (pulse_2_length[individual]+pulse_2_init[individual])){//other case where we are not at the first pulse point, but the pulse is still happenning
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_2_Cortisol_effect[individual];
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_2_CCS_effect[individual];
      } else if (time <= (pulse_2_length[individual]+pulse_2_init[individual] + 1) && time > (pulse_2_length[individual]+pulse_2_init[individual])){//other case where we are not at the final point
        Measurement_guess[time,individual,1] = Measurement_guess[time,individual,1] + pulse_2_Cortisol_effect[individual] * (-time + (1+pulse_2_length[individual]+pulse_2_init[individual]));
        Measurement_guess[time,individual,6] = Measurement_guess[time,individual,6] + pulse_2_CCS_effect[individual] * (-time + (1+pulse_2_length[individual]+pulse_2_init[individual]));
      }
      
      for (hormone in {1,6}){
        //eps_raw[time,individual,hormone] ~ normal(0,1);
        if (measurements_indicator[individual,time,hormone] == 1){
          eps[time,individual,hormone] = log(measurements[individual,time,hormone] / Measurement_guess[time,individual,hormone]);
          Measurement_guess[time,individual,hormone] = measurements[individual,time,hormone];
          //target += normal_lpdf(eps[time,individual,hormone] | 0, (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone])));
        } else {
          eps[time,individual,hormone] = eps_raw[time,individual,hormone] * (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone]));
          Measurement_guess[time,individual,hormone] = Measurement_guess[time,individual,hormone] * exp(eps[time,individual,hormone]);
        }
      }
      
      Measurement_guess[time,individual,2] = decay[individual,2] *prev_points[2] + IIBetaHSD2[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,3] = decay[individual,3] *prev_points[3] + VBetaReductase[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,4] = decay[individual,4] *prev_points[4] + VAlphaReductase[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,5] = decay[individual,5] *prev_points[5] + XVIIIHydroxylase[individual] * Measurement_guess[time,individual,1];
      Measurement_guess[time,individual,7] = decay[individual,7] *prev_points[7] + XIBetaHydroxylase[individual] * Measurement_guess[time,individual,6];

  
      
      for (hormone in {2,3,4,5,7}){
        //eps_raw[time,individual,hormone] ~ normal(0,1);
        if (measurements_indicator[individual,time,hormone] == 1){
          eps[time,individual,hormone] = log(measurements[individual,time,hormone] / Measurement_guess[time,individual,hormone]);
          Measurement_guess[time,individual,hormone] = measurements[individual,time,hormone];
          //target += normal_lpdf(eps[time,individual,hormone] | 0, (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone])));
        } else {
          eps[time,individual,hormone] = eps_raw[time,individual,hormone] * (sigma[individual,hormone]*exp(Measurement_guess[time,individual,hormone]*-theta[hormone]));
          Measurement_guess[time,individual,hormone] = Measurement_guess[time,individual,hormone] * exp(eps[time,individual,hormone]);
        }
      }
      
      
      prev_points = to_vector(Measurement_guess[time, individual, ]);
    }
  }
}


