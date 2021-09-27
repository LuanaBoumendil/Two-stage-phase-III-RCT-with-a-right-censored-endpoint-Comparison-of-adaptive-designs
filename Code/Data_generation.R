library(tidyverse)

data_generation <- function(N, lambda0, lambda1, K, K0, end_accrual_after_nb_month, end_study_after_nb_month) {
  tps_suivi <- NULL
  tt <- NULL
  trait <- rbinom(N, 1, 0.5)
  recrutement <- runif(N, 0, end_accrual_after_nb_month)
  tps_suivi[trait == 0] <- rweibull(sum(trait == 0), shape = K0, scale = (1 / lambda0))
  tps_suivi[trait == 1] <- rweibull(sum(trait == 1), shape = K, scale = (1 / lambda1)) 
  time_since_study_begening <- tps_suivi + recrutement
  evt <- ifelse(time_since_study_begening > end_study_after_nb_month, 0, 1)
  tt <- ifelse(evt == 1, tps_suivi, end_study_after_nb_month - recrutement)
  df <- data.frame(trait, recrutement, tps_suivi, evt, tt, time_since_study_begening)
  return(df)
}