# Loading packages ####
library(tidyverse)
library(survival)
library(rpact)

# Handicap values (found previously by simulation)####
Handicap0.5_N86_d66_K1 <- c(0.101, 0.094, 0.094, 0.093, 0.085, 0.085, 0.09,
                            0.083, 0.083, 0.082, 0.085, 0.049, 0.048, 0, 0,
                            0, 0)
Handicap0.67_N234_d192_K1 <- c(0.095, 0.099, 0.107, 0.116, 0.119, 0.125, 0.127,
                               0.127, 0.117, 0.117, 0.113, 0.095, 0.082, 0.073,
                               0.041, 0, 0)
Handicap0.8_N748_d633_K1 <- c(0.1, 0.109, 0.118, 0.127, 0.129, 0.131, 0.133,
                              0.133, 0.127, 0.12, 0.127, 0.118, 0.098, 0.093,
                              0.078, 0.04, 0)


# WHAT YOU HAVE PLANNED : Data simulation ####
end_accrual_after_nb_month <- 40
end_study_after_nb_month <- 60
tps_restant_apres_dernier_recrutement <- end_study - ta
lambda0 <- log(2) / 12
beta <- 0.2
alpha <- 0.05
K <- 1 # Weibull's shape parameter for the times of events in the experimental arm
#lambda1 <- log (2)^ (1/K) / 15
#lambda1 <- log (2)^ (1/K) / 18
lambda1 <- log(2)^ (1 / K) / 24
theta <- log(lambda1 / lambda0) # log(HR)
exp(theta) #HR


K <- 1

Get_Sample_Size_for_TAU0.5 <- function(alpha, beta, lambda0, lambda1,
                                       end_study_after_nb_month, end_accrual_after_nb_month) {
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

# SIMULATIONS
K <- 1
TAU <- c(seq(0.2, 0.95, 0.05))

Get_Analysis_Time <- function(TAU, alpha, beta, lambda0, lambda1, end_study_after_nb_month, end_accrual_after_nb_month) {
  IA <- c(rep(NA, length(TAU)))
  for (b in seq_len(length(TAU))) {
    tau <- TAU[b]
    schema <- getDesignGroupSequential(kMax = 2,
                                       alpha = alpha,
                                       sided = 2,
                                       beta = beta,
                                       informationRates = c(tau, 1),
                                       typeOfDesign = "asOF" #O'Brien & Fleming type alpha spending
    )
    DesignPlan <- getSampleSizeSurvival(design = schema,
                                        typeOfComputation = "Schoenfeld",
                                        thetaH0 = 1,
                                        lambda1 = lambda1,
                                        lambda2 = lambda0,
                                        kappa = 1,
                                        allocationRatioPlanned = 1,
                                        accrualTime = c(0, end_accrual_after_nb_month),
                                        followUpTime = end_study_after_nb_month - end_accrual_after_nb_month)
    IA[b] <- DesignPlan$analysisTime[1]
  }
  TAU <- c(TAU, 1)
  IA <- c(IA, 60)
  return(data.frame(TAU, IA))
}
TAU <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                    end_study_after_nb_month = 60, end_accrual_after_nb_month = 40)$TAU)
IA <- as.numeric(Get_Analysis_Time(TAU = c(seq(0.2, 0.95, 0.05)), alpha, beta, lambda0, lambda1,
                                   end_study_after_nb_month = 60, end_accrual_after_nb_month = 40)$IA)


# WHAT REALLY HAPPENS ####

lambda1 <- lambda0 # If under the null hypothesis
theta <- log(lambda1 / lambda0)
K0 <- 1
K <- 2 # Crossing hazard
nb_simu <- 1

Get_reject_proportion <- function(nb_simu, N, d, lambda0, lambda1, K, K0, end_accrual_after_nb_month,
                                  end_study_after_nb_month, TAU, IA, vect_handicap) {
  ta <- end_accrual_after_nb_month
  end_study <- end_study_after_nb_month
  var <- 4 # approximation
  early_stop <- c()
  power <- c()
  for (b in seq_len(length(TAU))) {
    tau <- TAU[b]
    tps_interim <- IA[b]
    handicap <- Handicap0.5_N86_d66_K1BIS[b]
    n <- c(tau * d, d)
    n0 <- var / (sqrt(var / (handicap * d))) ^ 2
    var_p <- var / (n + n0)
    Borne <- (var / n) * ((0 / var_p) - (theta0 * n0 / var) + qnorm(1 - 0.025) / sqrt(var_p))
    efficacy <- c()
    decision <- c()
    set.seed(27011996)
    for (i in 1:nb_simu) {
      data <- data_generation(N, lambda0, lambda1, K, K0, ta, end_study)
      d1 <- dim(filter(filter(data, evt == 1), time_since_study_begening < tps_interim))[1]
      if (tau != 1 & d1 < d) {
        fit1 <- coxph(Surv(tt, evt)~trait, data = filter(data, recrutement < tps_interim) %>%
                        mutate(evt = ifelse(time_since_study_begening >= tps_interim, 0, evt)) %>%
                        mutate(tt = ifelse(time_since_study_begening >= tps_interim, tps_interim - recrutement, tt)))
        fit2 <- coxph(Surv(tt, evt)~trait, data = data)
        Zstat1 <- summary(fit1)$coef[1]
        Zstat2 <- summary(fit2)$coef[1]
        if (abs(Zstat1) > Borne[1]) {
          efficacy[i] <- 1
          decision[i] <- 1
        } else {
          efficacy[i] <- 0
          decision[i] <- ifelse(abs(Zstat2) > Borne[2], 1, 0)
        }
      } else {
        fit2 <- coxph(Surv(tt, evt)~trait, data = data)
        Zstat2 <- summary(fit2)$coef[1]
        decision[i] <- ifelse(abs(Zstat2) > qnorm(1 - alpha / 2) / sqrt(d / var), 1, 0)
        efficacy[i] <- 0
      }
    }
    early_stop[b] <- mean(efficacy)
    power[b] <- mean(decision)
  }
  return(data.frame(TAU, vect_handicap, early_stop, power))
}


Get_reject_proportion(nb_simu = 10, N, d, lambda0, lambda1, K = 1, K0 = 1, 40,
                      60, TAU = TAU, IA = IA,
                      vect_handicap = Handicap0.5_N86_d66_K1)


lint(filename = "../../../M2 SMDS/Stage M2/Two-stage-phase-3-RCT-with-a-right-censored-endpoint-Comparison-of-sequential-designs/Code/Bayesian_approach_Grossman.R")
