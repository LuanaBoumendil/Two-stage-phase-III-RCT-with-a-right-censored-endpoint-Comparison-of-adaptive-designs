# Loading packages ####
library(tidyverse)
library(rpact)
library(survival)
library(nph)


# WHAT YOU HAVE PLANNED : Data simulation ####
end_accrual_after_nb_month <- 40
end_study_after_nb_month <- 60
tps_restant_apres_dernier_recrutement <- end_study - ta
lambda0 <- log(2) / 12
beta <- 0.2
alpha <- 0.05
K <- 1 # Weibull's shape parameter for the times of events in the experimental arm
median1 <- 24
#median1 <- 18
#median1 <- 15
lambda1 <- log(2)^ (1 / K) / median1
theta <- log(lambda1 / lambda0) # log(HR)
exp(theta) #HR


K <- 1

Get_Sample_Size_for_TAU0.5 <- function(alpha, beta, lambda0, lambda1,
                                       end_study_after_nb_month,
                                       end_accrual_after_nb_month) {
  tau <- 0.5
  schema <- getDesignGroupSequential(kMax = 2,
                                     alpha = alpha,
                                     sided = 2,
                                     beta = beta,
                                     informationRates = c(tau, 1),
                                     typeOfDesign = "asOF"      #O'Brien & Fleming type alpha spending
  )
  DesignPlan <- getSampleSizeSurvival(design = schema,
                                      typeOfComputation = "Schoenfeld",
                                      thetaH0 = 1,
                                      lambda1 = lambda1, #instantaneous risk in the experimental arm
                                      lambda2 = lambda0, #instantaneous risk in the control arm
                                      kappa = 1,
                                      allocationRatioPlanned = 1,
                                      accrualTime = c(0, end_accrual_after_nb_month),
                                      followUpTime = end_study_after_nb_month - end_accrual_after_nb_month)
  d <- ceiling(DesignPlan$maxNumberOfEvents)
  N <- ceiling(DesignPlan$maxNumberOfSubjects)
  if (N %% 2 != 0) {
    N <- N + 1
  }
  return(data.frame("Number of necessary events" = d,
                    "Number of patients to enroll" = N))
}

d <- as.numeric(Get_Sample_Size_for_TAU0.5(alpha = 0.05, beta = 0.2, lambda0, lambda1, 60, 40)[1])
N <- as.numeric(Get_Sample_Size_for_TAU0.5(alpha = 0.05, beta = 0.2, lambda0, lambda1, 60, 40)[2])

# SIMULATIONS ####
K <- 1
TAU <- c(seq(0.2, 0.95, 0.05))

Get_Analysis_Time <- function(TAU, alpha, beta, lambda0, lambda1, 
                              end_study_after_nb_month, end_accrual_after_nb_month) {
  IA <- c(rep(NA, length(TAU)))
  vect_early_stop <- c(rep(NA, length(TAU)))
  vect_final_stop <- c(rep(NA, length(TAU)))
  vect_c_alpha2 <- c(rep(NA, length(TAU)))
  for (b in seq_len(length(TAU))) {
    tau <- TAU[b]
    schema <- getDesignGroupSequential(kMax = 2,
                                       alpha = alpha,
                                       sided = 2,
                                       beta = beta,
                                       informationRates = c(tau, 1),
                                       typeOfDesign = "asOF"      #O'Brien & Fleming type alpha spending
    )
    DesignPlan <- getSampleSizeSurvival(design = schema,
                                        typeOfComputation = "Schoenfeld",
                                        thetaH0 = 1,
                                        lambda1 = lambda1, #instantaneous risk in the experimental arm
                                        lambda2 = lambda0, #instantaneous risk in the control arm
                                        kappa = 1,
                                        allocationRatioPlanned = 1,
                                        accrualTime = c(0, end_accrual_after_nb_month),
                                        followUpTime = end_study_after_nb_month - end_accrual_after_nb_month)
    IA[b] <- DesignPlan$analysisTime[1]
    vect_early_stop[b] <- DesignPlan$criticalValuesPValueScale[1]
    vect_final_stop[b] <- DesignPlan$criticalValuesPValueScale[2]
    vect_c_alpha2[b] <- exp(- qchisq((1 - DesignPlan$criticalValuesPValueScale[2]), df = 4, ncp = 0) / 2)
  }
  TAU <- c(TAU, 1)
  IA <- c(IA, 60)
  vect_early_stop <- c(vect_early_stop, NA)
  vect_final_stop <- c(vect_final_stop, alpha)
  vect_c_alpha2 <- c(vect_c_alpha2, exp(- qchisq((1 - vect_final_stop[length(TAU)]), df = 4, ncp = 0) / 2))
  return(data.frame(TAU, IA, vect_early_stop, vect_final_stop, vect_c_alpha2))
}
TAU <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                    end_study_after_nb_month = 60, end_accrual_after_nb_month = 40)$TAU)
IA <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                   end_study_after_nb_month = 60, end_accrual_after_nb_month = 40)$IA)
vect_early_stop <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                                end_study_after_nb_month = 60,
                                                end_accrual_after_nb_month = 40)$vect_early_stop)
vect_final_stop <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                                end_study_after_nb_month = 60,
                                                end_accrual_after_nb_month = 40)$vect_final_stop)
vect_c_alpha2 <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                              end_study_after_nb_month = 60,
                                              end_accrual_after_nb_month = 40)$vect_c_alpha2)

# WHAT REALLY HAPPENS ####

K <- 1
# To simulate under the alternative hypothesis
lambda1 <- log(2)^ (1 / K) / median1

# To simulate under the null for K=1 (proportional risks)
lambda1 <- lambda0 # If under the null hypothesis

# To simulate under the null for K!=1 (non proportional risks)
median0 <- 12 # median survival in the control arm
median1 <- median0 / (log(2)^ ((K - 1) / K))
lambda1 <- log(2)^ (1 / K) / median1

Get_reject_proportion <- function(nb_simu, N, d, lambda0, lambda1, K, K0, end_accrual_after_nb_month,
                                  end_study_after_nb_month, TAU, IA, vect_early_stop, vect_final_stop,
                                  vect_c_alpha2) {
  ta <- end_accrual_after_nb_month
  end_study <- end_study_after_nb_month
  efficacy <- c(rep(NA, length(TAU)))
  DESSEAUX <- c(rep(NA, length(TAU)))
  WASSMER <- c(rep(NA, length(TAU)))
  JORGENS <- c(rep(NA, length(TAU)))
  nb_simu <- 1
  for (b in seq_len(length(TAU))) {
    tau <- TAU[b]
    tps_interim <- IA[b]
    early_stop <- vect_early_stop[b]
    final_stop <- vect_final_stop[b]
    c_alpha2 <- vect_c_alpha2[b]
    w1 <- sqrt((tau * d) / d) # proportion of events expected to be observed at IA
    w2 <- sqrt(1 - w1^2)
    prop_rejet_WASSMER <- c(rep(NA, nb_simu))
    prop_efficacy <- c(rep(NA, nb_simu))
    prop_rejet_DESSEAUX <- c(rep(NA, nb_simu))
    prop_rejet_JORGENS <- c(rep(NA, nb_simu))
    set.seed(27011996)
    for (i in 1:nb_simu) {
      data <- data_generation(N, lambda0, lambda1, K, K0, ta, end_study)
      d1 <- dim(filter(filter(data, evt == 1), time_since_study_begening < tps_interim))[1] #number of events observed at IA
      d2 <- sum(data$evt)
      if (tau != 1 & d1 < d) {
        sous_base <- filter(data, recrutement < tps_interim) %>%
          mutate(evt = ifelse(time_since_study_begening >= tps_interim, 0, evt)) %>%
          mutate(tt = ifelse(time_since_study_begening >= tps_interim, tps_interim - recrutement, tt))
        LR_1 <- as.numeric(logrank.test(time = sous_base$tt,
                                        event = sous_base$evt,
                                        group = sous_base$trait,
                                        alternative = "two.sided")$test[3])
        p1 <- 2 * (1 - pnorm(abs(LR_1)))
        LR_2 <- as.numeric(logrank.test(time = data$tt,
                                        event = data$evt,
                                        group = data$trait,
                                        alternative = "two.sided")$test[3])
        Z2 <- ((sqrt(d2) * LR_2) - (sqrt(d1) * LR_1)) / (sqrt(d2 - d1))
        p2 <- 2 * (1 - pnorm(abs(Z2), mean = 0, sd = 1, log = FALSE)) # *2 because the test is two-sided
        
        if (p1 < early_stop) {
          early_stopping_D <- "efficacy"
          decision_D <- "reject H0"
          COMBI_D <- NA
        } else{
          early_stopping_D <- "NO"
          COMBI_D <- p1 * p2
          if (COMBI_D < c_alpha2) {
            decision_D <- "reject H0"
          }else{
            decision_D <- "don't reject H0"
          }
        }
        prop_rejet_DESSEAUX[i] <- ifelse(decision_D == "reject H0", 1, 0)
        prop_efficacy[i] <- ifelse(early_stopping_D == "efficacy", 1, 0)
        
        if (p1 < early_stop) {
          decision <- "reject H0"
        }else{
          COMBI <-  ((w1 * LR_1) + (w2 * Z2))
          p <- 2 * (1 - pnorm(abs(COMBI), 0, 1))
          decision <- ifelse(p < final_stop, "reject H0", "don't reject H0")
        }
        prop_rejet_WASSMER[i] <- ifelse(decision == "reject H0", 1, 0)
        
        # Expected number of events
        e_11 <- tau * d
        e_12 <-  round(dim(filter(data, recrutement < tps_interim))[1] / N * d) -  e_11
        e_22 <- pmax(0, d - e_11 - e_12)
        # Jorgens preplanned Weights
        w11 <- sqrt(e_11 / (e_11 + e_12 + e_22))
        if (e_12 <= 0) {
          w12 <- 0
        }else{
          w12 <- sqrt((e_12 / e_11) * w11^2)
        }
        if (e_22 != 0) {
          w22 <- sqrt(1 - (w11^2 + w12^2))
        }
        # Jorgens' decision rules
        if (p1 < early_stop) {
          p_final <- NA
          decision_J <- "reject H0"
        } else{
          data_22 <- filter(data, recrutement >= tps_interim)
          
          if (sum(((table(data_22$trait, data_22$evt) == 0) == T)) != 0) {
            COMBI <-  ((w1 * LR_1) + (w2 * Z2))
            p_final <- 2 * (1 - pnorm(abs(COMBI), 0, 1))
            decision_J <- ifelse(p_final < final_stop, "reject H0", "don't reject H0")
          }else{
            data_12 <- filter(data, recrutement < tps_interim)
            LR_12 <-  as.numeric(logrank.test(time = data_12$tt,
                                              event = data_12$evt,
                                              group = data_12$trait,
                                              alternative = "two.sided")$test[3])
            d_12 <- sum(data_12$evt)
            Z12 <- ((sqrt(d_12) * LR_12) - (sqrt(d1) * LR_1)) / (sqrt(d_12 - d1))
            if (e_22 > 0 & dim(table(filter(data, recrutement >= tps_interim)$trait)) == 2 &
                dim(table(filter(data, recrutement >= tps_interim)$evt)) == 2) {
              data_22 <- filter(data, recrutement >= tps_interim)
              LR_22 <- as.numeric(logrank.test(time = data_22$tt,
                                               event = data_22$evt,
                                               group = data_22$trait,
                                               alternative = "two.sided")$test[3])
              COMBI_J <-  (w11 * LR_1) + (w12 * Z12) + (w22 * LR_22)
            }else{
              COMBI_J <- (w11 * LR_1) + (w12 * Z12)
            }
            p_final <- 2 * (1 - pnorm(abs(COMBI_J), 0, 1))
            decision_J <- ifelse(p_final < final_stop, "reject H0", "don't reject H0")
          }
        }
        prop_rejet_JORGENS[i] <- ifelse(decision_J == "reject H0", 1, 0)
      }else{ #tau!=1
        LR_2 <- as.numeric(logrank.test(time = data$tt,
                                        event = data$evt,
                                        group = data$trait,
                                        alternative = "two.sided")$test[3])
        p2 <- 2 * (1 - pnorm(abs(LR_2)))
        if (p2 < alpha) {
          decision <- "reject H0"
        }else{
          decision <- "don't reject H0"
        }
        prop_rejet <- ifelse(decision == "reject H0", 1, 0)
        prop_efficacy[i] <- 0
        prop_rejet_DESSEAUX[i] <- prop_rejet
        prop_rejet_WASSMER[i] <- prop_rejet
        prop_rejet_JORGENS[i] <- prop_rejet
      }
    }
    efficacy[b] <- mean(na.omit(prop_efficacy))
    DESSEAUX[b] <- mean(prop_rejet_DESSEAUX)
    WASSMER[b] <- mean(prop_rejet_WASSMER)
    JORGENS[b] <- mean(prop_rejet_JORGENS)
  }
  return(data.frame(TAU, efficacy, DESSEAUX, WASSMER, JORGENS))
}

Get_reject_proportion(nb_simu = 100, N, d, lambda0, lambda1, K = 1, K0 = 1, 40,
                      60, TAU = TAU, IA = IA,
                      vect_early_stop,
                      vect_final_stop,
                      vect_c_alpha2)