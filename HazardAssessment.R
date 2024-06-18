# -*- coding: utf-8 -*-

# Author: Ziying Zhou
# Date: June 8, 2024
# Description: This script calculates joint probabilities, joint return periods, compound hazard index, and failure probabilities based on C-Vine Copula functions.

# Initialization Settings ####################################################

library(VineCopula)
library(CDVineCopulaConditional)

# Define the directories for input and output data
copula_dir <- "D:\\Project_CompoundHazard\\Data\\Copula"
probmargin_dir <- file.path(copula_dir, "ProbMargin")
RVM_dir <- file.path(copula_dir, "RVM")
bestcopula_dir <- file.path(copula_dir, "BestCopula")
returnperiod_dir <- file.path(copula_dir, "ReturnPeriod")
random_dir <- file.path(copula_dir, "Random")
failureprob_dir <- file.path(copula_dir, "FailureProb")
servtime_dir <- file.path(copula_dir, "ServiceTime")

# Path to the station list
stationlist_path <- "StationList.csv"

# Read the station list
dfStation <- read.csv(stationlist_path, stringsAsFactors = FALSE)

# Define the period of years
yearPeriod <- 2018 - 1979 + 1

# Function Definitions #######################################################

### JointProb_3D
JointProb_3D <- function(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS) {
  cdf_WR <- BiCopCDF(w, r, family = dist_WR, par = para1_WR, par2 = para2_WR, obj = NULL, check.pars = TRUE)
  cond_R.W <- cdf_WR / w
  cdf_WS <- BiCopCDF(w, s, family = dist_WS, par = para1_WS, par2 = para2_WS, obj = NULL, check.pars = TRUE)
  cond_S.W <- cdf_WS / w
  cdf_RS <- BiCopCDF(r, s, family = dist_RS, par = para1_RS, par2 = para2_RS, obj = NULL, check.pars = TRUE)
  cond_RS.W <- BiCopCDF(cond_R.W, cond_S.W, family = dist_RS, par = para1_RS, par2 = para2_RS, obj = NULL, check.pars = TRUE)
  cdf_3D <- cond_RS.W * w
  return(c(cdf_3D, cdf_WR, cdf_WS, cdf_RS))
}

### IndepExceed_OR
IndepExceed_OR <- function(w, r, s) {
  Exceed_3D_OR <- 1 - w * r * s
  Exceed_WR_OR <- 1 - w * r
  Exceed_WS_OR <- 1 - w * s
  Exceed_RS_OR <- 1 - r * s
  return(c(Exceed_3D_OR, Exceed_WR_OR, Exceed_WS_OR, Exceed_RS_OR))
}

### IndepExceed_AND
IndepExceed_AND <- function(w, r, s) {
  Exceed_3D_AND <- (1 - w) * (1 - r) * (1 - s)
  Exceed_WR_AND <- (1 - w) * (1 - r)
  Exceed_WS_AND <- (1 - w) * (1 - s)
  Exceed_RS_AND <- (1 - r) * (1 - s)
  return(c(Exceed_3D_AND, Exceed_WR_AND, Exceed_WS_AND, Exceed_RS_AND))
}

### JointExceed_OR
JointExceed_OR <- function(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS) {
  cdf_results <- JointProb_3D(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)
  cdf_3D <- cdf_results[1]
  cdf_WR <- cdf_results[2]
  cdf_WS <- cdf_results[3]
  cdf_RS <- cdf_results[4]
  Exceed_3D_OR <- 1 - cdf_3D
  Exceed_WR_OR <- 1 - cdf_WR
  Exceed_WS_OR <- 1 - cdf_WS
  Exceed_RS_OR <- 1 - cdf_RS
  return(c(Exceed_3D_OR, Exceed_WR_OR, Exceed_WS_OR, Exceed_RS_OR))
}

### JointExceed_AND
JointExceed_AND <- function(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS) {
  cdf_results <- JointProb_3D(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)
  cdf_3D <- cdf_results[1]
  cdf_WR <- cdf_results[2]
  cdf_WS <- cdf_results[3]
  cdf_RS <- cdf_results[4]
  Exceed_3D_AND <- 1 - w - r - s + cdf_WR + cdf_WS + cdf_RS - cdf_3D
  Exceed_WR_AND <- 1 - w - r + cdf_WR
  Exceed_WS_AND <- 1 - w - s + cdf_WS
  Exceed_RS_AND <- 1 - r - s + cdf_RS
  return(c(Exceed_3D_AND, Exceed_WR_AND, Exceed_WS_AND, Exceed_RS_AND))
}

### ServiceTime
ServiceTime <- function(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU, FP) {
  IndepExceed <- 1 - (1 - FP)^(MU / lifetime)
  JointExceed <- JointExceed_OR((1 - IndepExceed), (1 - IndepExceed), (1 - IndepExceed), dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[1]
  return(c(IndepExceed, JointExceed))
}

### RandomGenerate
RandomGenerate <- function(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, RVM) {
  set.seed(123)
  Sim <- RVineSim(1000, RVM)
  results_RG <- data.frame(w = numeric(), r = numeric(), s = numeric(), cdf_3D = numeric(), cdf_WR = numeric(), cdf_WS = numeric(), cdf_RS = numeric(), stringsAsFactors = FALSE)
  for (i in 1:nrow(Sim)) {
    w <- as.vector(Sim[i, 3])
    r <- as.vector(Sim[i, 1])
    s <- as.vector(Sim[i, 2])
    cdf_results <- JointProb_3D(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)
    results_RG <- rbind(results_RG, data.frame(w = w, r = r, s = s, cdf_3D = cdf_results[1], cdf_WR = cdf_results[2], cdf_WS = cdf_results[3], cdf_RS = cdf_results[4], stringsAsFactors = FALSE))
  }
  return(results_RG)
}

### UnivFailure
UnivFailure <- function(w, lifetime, MU) {
  Prob <- w
  fail <- 1 - (Prob)^(lifetime / MU)
  return(fail)
}

### IndepFailure_OR
IndepFailure_OR <- function(w, r, s, lifetime, MU) {
  cdf_3D <- IndepExceed_OR(w, r, s)[1]
  fail_3D <- 1 - (1 - cdf_3D)^(lifetime / MU)
  cdf_WR <- IndepExceed_OR(w, r, s)[2]
  fail_WR <- 1 - (1 - cdf_WR)^(lifetime / MU)
  cdf_WS <- IndepExceed_OR(w, r, s)[3]
  fail_WS <- 1 - (1 - cdf_WS)^(lifetime / MU)
  cdf_RS <- IndepExceed_OR(w, r, s)[4]
  fail_RS <- 1 - (1 - cdf_RS)^(lifetime / MU)
  return(c(fail_3D, fail_WR, fail_WS, fail_RS))
}

### IndepFailure_AND
IndepFailure_AND <- function(w, r, s, lifetime, MU) {
  cdf_3D <- IndepExceed_AND(w, r, s)[1]
  fail_3D <- 1 - (1 - cdf_3D)^(lifetime / MU)
  cdf_WR <- IndepExceed_AND(w, r, s)[2]
  fail_WR <- 1 - (1 - cdf_WR)^(lifetime / MU)
  cdf_WS <- IndepExceed_AND(w, r, s)[3]
  fail_WS <- 1 - (1 - cdf_WS)^(lifetime / MU)
  cdf_RS <- IndepExceed_AND(w, r, s)[4]
  fail_RS <- 1 - (1 - cdf_RS)^(lifetime / MU)
  return(c(fail_3D, fail_WR, fail_WS, fail_RS))
}

### JointFailure_OR
JointFailure_OR <- function(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU) {
  cdf_3D <- JointExceed_OR(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[1]
  fail_3D <- 1 - (1 - cdf_3D)^(lifetime / MU)
  cdf_WR <- JointExceed_OR(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[2]
  fail_WR <- 1 - (1 - cdf_WR)^(lifetime / MU)
  cdf_WS <- JointExceed_OR(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[3]
  fail_WS <- 1 - (1 - cdf_WS)^(lifetime / MU)
  cdf_RS <- JointExceed_OR(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[4]
  fail_RS <- 1 - (1 - cdf_RS)^(lifetime / MU)
  return(c(fail_3D, fail_WR, fail_WS, fail_RS))
}

### JointFailure_AND
JointFailure_AND <- function(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU) {
  cdf_3D <- JointExceed_AND(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[1]
  fail_3D <- 1 - (1 - cdf_3D)^(lifetime / MU)
  cdf_WR <- JointExceed_AND(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[2]
  fail_WR <- 1 - (1 - cdf_WR)^(lifetime / MU)
  cdf_WS <- JointExceed_AND(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[3]
  fail_WS <- 1 - (1 - cdf_WS)^(lifetime / MU)
  cdf_RS <- JointExceed_AND(w, r, s, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[4]
  fail_RS <- 1 - (1 - cdf_RS)^(lifetime / MU)
  return(c(fail_3D, fail_WR, fail_WS, fail_RS))
}

# Probability Calculations ####################################################

RP_OR_path <- file.path(returnperiod_dir, "ReturnPeriod_OR.csv")
RP_AND_path <- file.path(returnperiod_dir, "ReturnPeriod_AND.csv")

results_RP_OR <- data.frame(ID = character(),
                            ExceedProb = numeric(),
                            UnivRP = numeric(),
                            IndepRP_3D = numeric(),
                            IndepRP_WR = numeric(),
                            IndepRP_WS = numeric(),
                            IndepRP_RS = numeric(),
                            JointRP_3D = numeric(),
                            JointRP_WR = numeric(),
                            JointRP_WS = numeric(),
                            JointRP_RS = numeric(),
                            stringsAsFactors = FALSE)

results_RP_AND <- data.frame(ID = character(),
                             ExceedProb = numeric(),
                             UnivRP = numeric(),
                             IndepRP_3D = numeric(),
                             IndepRP_WR = numeric(),
                             IndepRP_WS = numeric(),
                             IndepRP_RS = numeric(),
                             JointRP_3D = numeric(),
                             JointRP_WR = numeric(),
                             JointRP_WS = numeric(),
                             JointRP_RS = numeric(),
                             stringsAsFactors = FALSE)

for(n in 1:nrow(dfStation)) {
  dfTemp <- dfStation[n, ]
  id <- as.character(dfTemp["ID"])
  station_id <- substr(as.character(dfTemp["ID"] + 100), 2, 3)
  
  probmargin_path <- file.path(probmargin_dir, paste0("ProbMargin", station_id, ".csv"))
  RVM_path <- file.path(RVM_dir, paste0("RVM", station_id, ".RData"))
  bestcopula_path <- file.path(bestcopula_dir, paste0("BestCopula", station_id, ".csv"))
  
  data_probmargin <- read.csv(probmargin_path)
  data_bestcopula <- read.csv(bestcopula_path)
  load(file = RVM_path)
  
  data_WR <- data_bestcopula[1, ]
  dist_WR <- data_WR[["Family"]]
  para1_WR <- data_WR[["Par1"]]
  para2_WR <- data_WR[["Par2"]]
  
  data_WS <- data_bestcopula[2, ]
  dist_WS <- data_WS[["Family"]]
  para1_WS <- data_WS[["Par1"]]
  para2_WS <- data_WS[["Par2"]]
  
  data_RS <- data_bestcopula[3, ]
  dist_RS <- data_RS[["Family"]]
  para1_RS <- data_RS[["Par1"]]
  para2_RS <- data_RS[["Par2"]]
  
  MU <- yearPeriod / nrow(data_probmargin)
  
  # Return Period Calculations ==============================================
  
  # Probability P to Return Period T
  P <- 0.95
  UnivRP <- MU / (1 - P)
  
  ## Scenario OR ------------------------------------------------
  
  ### IndepExceed_OR
  IndepExceed_3D_OR <- IndepExceed_OR(P, P, P)[1]
  IndepRP_3D_OR <- MU / IndepExceed_3D_OR
  
  IndepExceed_WR_OR <- IndepExceed_OR(P, P, P)[2]
  IndepRP_WR_OR <- MU / IndepExceed_WR_OR
  
  IndepExceed_WS_OR <- IndepExceed_OR(P, P, P)[3]
  IndepRP_WS_OR <- MU / IndepExceed_WS_OR
  
  IndepExceed_RS_OR <- IndepExceed_OR(P, P, P)[4]
  IndepRP_RS_OR <- MU / IndepExceed_RS_OR
  
  ### JointExceed_OR
  JointExceed_3D_OR <- JointExceed_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[1]
  JointRP_3D_OR <- MU / JointExceed_3D_OR
  
  JointExceed_WR_OR <- JointExceed_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[2]
  JointRP_WR_OR <- MU / JointExceed_WR_OR
  
  JointExceed_WS_OR <- JointExceed_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[3]
  JointRP_WS_OR <- MU / JointExceed_WS_OR
  
  JointExceed_RS_OR <- JointExceed_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[4]
  JointRP_RS_OR <- MU / JointExceed_RS_OR
  
  ### Export Results OR
  results_RP_OR <- rbind(results_RP_OR, data.frame(ID = id,
                                                   ExceedProb = 1 - P,
                                                   UnivRP = UnivRP,
                                                   IndepRP_3D = IndepRP_3D_OR,
                                                   IndepRP_WR = IndepRP_WR_OR,
                                                   IndepRP_WS = IndepRP_WS_OR,
                                                   IndepRP_RS = IndepRP_RS_OR,
                                                   JointRP_3D = JointRP_3D_OR,
                                                   JointRP_WR = JointRP_WR_OR,
                                                   JointRP_WS = JointRP_WS_OR,
                                                   JointRP_RS = JointRP_RS_OR))
  
  ## Scenario AND ------------------------------------------------
  
  ### IndepExceed_AND
  IndepExceed_3D_AND <- IndepExceed_AND(P, P, P)[1]
  IndepRP_3D_AND <- MU / IndepExceed_3D_AND
  
  IndepExceed_WR_AND <- IndepExceed_AND(P, P, P)[2]
  IndepRP_WR_AND <- MU / IndepExceed_WR_AND
  
  IndepExceed_WS_AND <- IndepExceed_AND(P, P, P)[3]
  IndepRP_WS_AND <- MU / IndepExceed_WS_AND
  
  IndepExceed_RS_AND <- IndepExceed_AND(P, P, P)[4]
  IndepRP_RS_AND <- MU / IndepExceed_RS_AND
  
  ### JointExceed_AND
  JointExceed_3D_AND <- JointExceed_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[1]
  JointRP_3D_AND <- MU / JointExceed_3D_AND
  
  JointExceed_WR_AND <- JointExceed_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[2]
  JointRP_WR_AND <- MU / JointExceed_WR_AND
  
  JointExceed_WS_AND <- JointExceed_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[3]
  JointRP_WS_AND <- MU / JointExceed_WS_AND
  
  JointExceed_RS_AND <- JointExceed_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS)[4]
  JointRP_RS_AND <- MU / JointExceed_RS_AND
  
  ### Export Results AND
  results_RP_AND <- rbind(results_RP_AND, data.frame(ID = id,
                                                     ExceedProb = 1 - P,
                                                     UnivRP = UnivRP,
                                                     IndepRP_3D = IndepRP_3D_AND,
                                                     IndepRP_WR = IndepRP_WR_AND,
                                                     IndepRP_WS = IndepRP_WS_AND,
                                                     IndepRP_RS = IndepRP_RS_AND,
                                                     JointRP_3D = JointRP_3D_AND,
                                                     JointRP_WR = JointRP_WR_AND,
                                                     JointRP_WS = JointRP_WS_AND,
                                                     JointRP_RS = JointRP_RS_AND))
  
  # Random Generation ========================================================
  
  RG_path <- file.path(random_dir, paste0("Random", station_id, ".csv"))
  results_RG <- RandomGenerate(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, RVM)
  write.csv(results_RG, RG_path, row.names = FALSE)
  print(RG_path)
  
  # Failure Probability Calculations =========================================
  
  listT_FP <- c(10, 20, 50, 100, 200)
  
  for (m in 1:length(listT_FP)) {
    T_FP <- listT_FP[m]
    lifetime_list <- seq(0, 100, by = 1)
    P <- -MU / T_FP + 1
    
    failureprob_sub_dir <- file.path(failureprob_dir, paste0("RP", substr(as.character(T_FP + 1000), 2, 4)))
    if (dir.exists(failureprob_sub_dir)) {
      # unlink(failureprob_sub_dir, recursive = TRUE)
    } else {
      dir.create(failureprob_sub_dir, recursive = TRUE)
    }
    print(sub_dir)
    
    FP_OR_path <- file.path(failureprob_sub_dir, paste0("FailureProb_OR_", station_id, ".csv"))
    FP_AND_path <- file.path(failureprob_sub_dir, paste0("FailureProb_AND_", station_id, ".csv"))
    
    results_FP_OR <- data.frame(lifetime = character(),
                                ExceedProb = numeric(),
                                UnivFP = numeric(),
                                IndepFP_3D = numeric(),
                                IndepFP_WR = numeric(),
                                IndepFP_WS = numeric(),
                                IndepFP_RS = numeric(),
                                JointFP_3D = numeric(),
                                JointFP_WR = numeric(),
                                JointFP_WS = numeric(),
                                JointFP_RS = numeric(),
                                stringsAsFactors = FALSE)
    
    results_FP_AND <- data.frame(lifetime = character(),
                                 ExceedProb = numeric(),
                                 UnivFP = numeric(),
                                 IndepFP_3D = numeric(),
                                 IndepFP_WR = numeric(),
                                 IndepFP_WS = numeric(),
                                 IndepFP_RS = numeric(),
                                 JointFP_3D = numeric(),
                                 JointFP_WR = numeric(),
                                 JointFP_WS = numeric(),
                                 JointFP_RS = numeric(),
                                 stringsAsFactors = FALSE)
    
    for (lifetime in lifetime_list) {
      
      ## Scenario OR ------------------------------------------------
      UnivFP <- UnivFailure(P, lifetime, MU)
      
      ### IndepFP_OR
      IndepFP_3D_OR <- IndepFailure_OR(P, P, P, lifetime, MU)[1]
      IndepFP_WR_OR <- IndepFailure_OR(P, P, P, lifetime, MU)[2]
      IndepFP_WS_OR <- IndepFailure_OR(P, P, P, lifetime, MU)[3]
      IndepFP_RS_OR <- IndepFailure_OR(P, P, P, lifetime, MU)[4]
      
      ### JointFP_OR
      JointFP_3D_OR <- JointFailure_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[1]
      JointFP_WR_OR <- JointFailure_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[2]
      JointFP_WS_OR <- JointFailure_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[3]
      JointFP_RS_OR <- JointFailure_OR(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[4]
      
      results_FP_OR <- rbind(results_FP_OR, data.frame(lifetime = lifetime,
                                                       ExceedProb = 1 - P,
                                                       UnivFP = UnivFP,
                                                       IndepFP_3D = IndepFP_3D_OR,
                                                       IndepFP_WR = IndepFP_WR_OR,
                                                       IndepFP_WS = IndepFP_WS_OR,
                                                       IndepFP_RS = IndepFP_RS_OR,
                                                       JointFP_3D = JointFP_3D_OR,
                                                       JointFP_WR = JointFP_WR_OR,
                                                       JointFP_WS = JointFP_WS_OR,
                                                       JointFP_RS = JointFP_RS_OR))
      
      ## Scenario AND ------------------------------------------------
      UnivFP <- UnivFailure(P, lifetime, MU)
      
      ### IndepFP_AND
      IndepFP_3D_AND <- IndepFailure_AND(P, P, P, lifetime, MU)[1]
      IndepFP_WR_AND <- IndepFailure_AND(P, P, P, lifetime, MU)[2]
      IndepFP_WS_AND <- IndepFailure_AND(P, P, P, lifetime, MU)[3]
      IndepFP_RS_AND <- IndepFailure_AND(P, P, P, lifetime, MU)[4]
      
      ### JointFP_AND
      JointFP_3D_AND <- JointFailure_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[1]
      JointFP_WR_AND <- JointFailure_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[2]
      JointFP_WS_AND <- JointFailure_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[3]
      JointFP_RS_AND <- JointFailure_AND(P, P, P, dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, lifetime, MU)[4]
      
      results_FP_AND <- rbind(results_FP_AND, data.frame(lifetime = lifetime,
                                                         ExceedProb = 1 - P,
                                                         UnivFP = UnivFP,
                                                         IndepFP_3D = IndepFP_3D_AND,
                                                         IndepFP_WR = IndepFP_WR_AND,
                                                         IndepFP_WS = IndepFP_WS_AND,
                                                         IndepFP_RS = IndepFP_RS_AND,
                                                         JointFP_3D = JointFP_3D_AND,
                                                         JointFP_WR = JointFP_WR_AND,
                                                         JointFP_WS = JointFP_WS_AND,
                                                         JointFP_RS = JointFP_RS_AND))
    }
    
    write.csv(results_FP_OR, FP_OR_path, row.names = FALSE)
    print(FP_OR_path)
    
    write.csv(results_FP_AND, FP_AND_path, row.names = FALSE)
    print(FP_AND_path)
  }
  
  # Service Time Calculations ===============================================
  
  ST_path <- file.path(servtime_dir, paste0("ServiceTime", station_id, ".csv"))
  
  listFP <- seq(from = 0.01, to = 0.99, by = 0.01)
  
  results_ST <- data.frame(FP = character(),
                           IndepExceed_LT010 = numeric(),
                           JointExceed_LT010 = numeric(),
                           IndepExceed_LT020 = numeric(),
                           JointExceed_LT020 = numeric(),
                           IndepExceed_LT050 = numeric(),
                           JointExceed_LT050 = numeric(),
                           IndepExceed_LT100 = numeric(),
                           JointExceed_LT100 = numeric(),
                           IndepExceed_LT200 = numeric(),
                           JointExceed_LT200 = numeric(),
                           stringsAsFactors = FALSE)
  
  for (j in 1:length(listFP)) {
    FP <- listFP[j]
    
    IndepExceed_LT010 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 10, MU, FP)[1]
    JointExceed_LT010 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 10, MU, FP)[2]
    
    IndepExceed_LT020 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 20, MU, FP)[1]
    JointExceed_LT020 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 20, MU, FP)[2]
    
    IndepExceed_LT050 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 50, MU, FP)[1]
    JointExceed_LT050 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 50, MU, FP)[2]
    
    IndepExceed_LT100 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 100, MU, FP)[1]
    JointExceed_LT100 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 100, MU, FP)[2]
    
    IndepExceed_LT200 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 200, MU, FP)[1]
    JointExceed_LT200 <- ServiceTime(dist_WR, para1_WR, para2_WR, dist_WS, para1_WS, para2_WS, dist_RS, para1_RS, para2_RS, 200, MU, FP)[2]
    
    results_ST <- rbind(results_ST, data.frame(FP = FP,
                                               IndepExceed_LT010 = IndepExceed_LT010,
                                               JointExceed_LT010 = JointExceed_LT010,
                                               IndepExceed_LT020 = IndepExceed_LT020,
                                               JointExceed_LT020 = JointExceed_LT020,
                                               IndepExceed_LT050 = IndepExceed_LT050,
                                               JointExceed_LT050 = JointExceed_LT050,
                                               IndepExceed_LT100 = IndepExceed_LT100,
                                               JointExceed_LT100 = JointExceed_LT100,
                                               IndepExceed_LT200 = IndepExceed_LT200,
                                               JointExceed_LT200 = JointExceed_LT200))
  }
  
  write.csv(results_ST, ST_path, row.names = FALSE)
  print(ST_path)
}

write.csv(results_RP_OR, RP_OR_path, row.names = FALSE)
print(RP_OR_path)

write.csv(results_RP_AND, RP_AND_path, row.names = FALSE)
print(RP_AND_path)