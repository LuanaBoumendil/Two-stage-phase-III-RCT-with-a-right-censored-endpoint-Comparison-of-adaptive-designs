library(tidyverse)

data_generation <- function(N # sample size
                            , lambda0 # baseline hazard and scale parameter of the Weibull distribution in the control arm
                            , lambda1 # baseline hazard and scale parameter of the Weibull distribution in the experimental arm
                            , K # shape parameter of the Weibull distribution in the control arm 
                            , K0 # shape parameter of the Weibull distribution in the experimental arm
                            , end_accrual_after_nb_month  # end of accruals
                            , end_study_after_nb_month) { # end of the studdy
  tps_suivi <- NULL
  tt <- NULL
  trait <- rbinom(N, 1, 0.5) # treatment arm
  recrutement <- runif(N, 0, end_accrual_after_nb_month) # accrual times
  tps_suivi[trait == 0] <- rweibull(sum(trait == 0), shape = K0, scale = (1 / lambda0))
  tps_suivi[trait == 1] <- rweibull(sum(trait == 1), shape = K, scale = (1 / lambda1)) 
  time_since_study_begening <- tps_suivi + recrutement
  evt <- ifelse(time_since_study_begening > end_study_after_nb_month, 0, 1) # administrative censoring
  tt <- ifelse(evt == 1, tps_suivi, end_study_after_nb_month - recrutement) # follow-up
  df <- data.frame(trait, recrutement, tps_suivi, evt, tt, time_since_study_begening)
  return(df)
}

# The data simulation function for stage 2 is similar but the time of accrual goes from the time of interim analysis to the end of the accrual period
data_generation_Stage2 <- function(N, # the additional number of patients to be recruited in stage 2
                                   lambda0, 
                                   lambda1, 
                                   K, 
                                   K0, 
                                   end_accrual_after_nb_month, 
                                   end_study_after_nb_month,
                                   AI) { # the duration between the begening of the trial and the interim analysis time in months
  tps_suivi <- NULL
  tt <- NULL
  trait <- rbinom(N, 1, 0.5)
  recrutement <- runif(N, AI, end_accrual_after_nb_month)
  tps_suivi[trait == 0] <- rweibull(sum(trait == 0), shape = K0, scale = (1 / lambda0))
  tps_suivi[trait == 1] <- rweibull(sum(trait == 1), shape = K, scale = (1 / lambda1)) 
  time_since_study_begening <- tps_suivi + recrutement
  evt <- ifelse(time_since_study_begening > end_study_after_nb_month, 0, 1)
  tt <- ifelse(evt == 1, tps_suivi, end_study_after_nb_month - recrutement)
  df <- data.frame(trait, recrutement, tps_suivi, evt, tt, time_since_study_begening)
  return(df)
}
