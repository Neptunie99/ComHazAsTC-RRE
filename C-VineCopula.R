# -*- coding: utf-8 -*-

# Author: Ziying Zhou
# Date: June 8, 2024
# Description: This script calculates the optimal C-Vine Copula function and its corresponding parameters.

# Initialization Settings ==================================================== #

library(VineCopula)
library(CDVineCopulaConditional)

# Define the directories for input and output data
copula_dir <- "E:\\Project_CompoundHazard\\Data\\Copula"
probmargin_dir <- file.path(copula_dir, "ProbMargin")
RVM_dir <- file.path(copula_dir, "RVM")
bestcopula_dir <- file.path(copula_dir, "BestCopula")

# Path to the station list
stationlist_path <- "StationList.csv"

# Read the station list
dfStation <- read.csv(stationlist_path, stringsAsFactors = FALSE)

# Calculate C-Vine Copula Functions ========================================== #

for(i in 1:nrow(dfStation)) {
  
  # Data Reading ------------------------------------------------------------
  dfTemp <- dfStation[i, ]

  station_id <- substr(as.character(dfTemp["ID"] + 100), 2, 3)
  
  probmargin_path <- file.path(probmargin_dir, paste0("ProbMargin", station_id, ".csv"))
  RVM_path <- file.path(RVM_dir, paste0("RVM", station_id, ".RData"))
  bestcopula_path <- file.path(bestcopula_dir, paste0("BestCopula", station_id, ".csv"))
  
  data_probmargin <- read.csv(probmargin_path)
  
  # Reorder the data (condition variable at the end)
  data_value <- data.frame(data_probmargin$Rain, data_probmargin$Surge, data_probmargin$Wind)
  
  # Function Fitting ---------------------------------------------------------
  # Define Vine structure
  matrix_vine <- matrix(c(1, 2, 3, 0, 2, 3, 0, 0, 3), nrow = 3)
  
  # Obtain Vine fitting
  RVM <- CDVineCondFit(data_value, Nx = 2, type = "CVine", selectioncrit = "AIC", rotations = FALSE, method = "mle", Matrix = matrix_vine)
  
  # Output Vine fitting results
  save(RVM, file = RVM_path)
  print(RVM)
  
  # Result Statistics --------------------------------------------------------
  # Capture print function output
  output <- capture.output(print(RVM))
  
  # Select specific lines from output
  selected_output <- output[c(3, 4, 7)]
  print(selected_output)
  
  pairs <- c("WR", "WS", "RS")
  copula_type <- c()
  
  for (line in selected_output) {
    # Split each line by the "(" symbol
    split_line <- strsplit(line, "\\(")[[1]]
    first_part <- split_line[1]
    
    # Process the first part to identify pairs and copula_type
    matches_first <- regmatches(first_part, regexec("^([0-9,;]+)\\s+(.*)", first_part))[[1]]
    if (length(matches_first) == 3) {
      copula_type <- c(copula_type, trimws(matches_first[3]))
    }
  }
  
  family_value <- c(RVM$family[3, 1], RVM$family[3, 2], RVM$family[2, 1])
  par1_value <- c(RVM$par[3, 1], RVM$par[3, 2], RVM$par[2, 1])
  par2_value <- c(RVM$par2[3, 1], RVM$par2[3, 2], RVM$par2[2, 1])
  tau_value <- c(RVM$tau[3, 1], RVM$tau[3, 2], RVM$tau[2, 1])
  emptau_value <- c(RVM$emptau[3, 1], RVM$emptau[3, 2], RVM$emptau[2, 1])
  
  # Correlation Calculations -------------------------------------------------
  cortest_WR <- cor.test(data_probmargin$Wind, data_probmargin$Rain, method = "kendall")
  cortest_WS <- cor.test(data_probmargin$Wind, data_probmargin$Surge, method = "kendall")
  
  result_WR <- c()
  for (j in 1:nrow(data_probmargin)) {
    u <- data_probmargin["Wind"][j, ]
    v <- data_probmargin["Rain"][j, ]
    p <- BiCopHfunc1(u, v, family = family_value[1], par = par1_value[1], par2 = par2_value[1])
    result_WR[j] <- p
  }
  data_probmargin$Cond_WR <- result_WR
  
  result_WS <- c()
  for (k in 1:nrow(data_probmargin)) {
    u <- data_probmargin["Wind"][k, ]
    v <- data_probmargin["Surge"][k, ]
    p <- BiCopHfunc1(u, v, family = family_value[2], par = par1_value[2], par2 = par2_value[2])
    result_WS[k] <- p
  }
  data_probmargin$Cond_WS <- result_WS
  
  cortest_RS <- cor.test(data_probmargin$Cond_WR, data_probmargin$Cond_WS, method = "kendall")
  
  p_value <- c(cortest_WR[["p.value"]], cortest_WS[["p.value"]], cortest_RS[["p.value"]])
  
  # AIC & BIC Calculations --------------------------------------------------
  copula_WR <- BiCopSelect(data_probmargin$Wind, data_probmargin$Rain, familyset = family_value[1], rotation = FALSE, method = "mle", selectioncrit = "AIC")
  copula_WS <- BiCopSelect(data_probmargin$Wind, data_probmargin$Surge, familyset = family_value[2], rotation = FALSE, method = "mle", selectioncrit = "AIC")
  copula_RS <- BiCopSelect(data_probmargin$Cond_WR, data_probmargin$Cond_WS, familyset = family_value[3], rotation = FALSE, method = "mle", selectioncrit = "AIC")
  
  aic_value <- c(copula_WR[["AIC"]], copula_WS[["AIC"]], copula_RS[["AIC"]])
  
  # Construct a data frame
  result_df <- data.frame(Pairs = pairs, Copula_Type = copula_type, Family = as.numeric(family_value), Par1 = as.numeric(par1_value), Par2 = as.numeric(par2_value), Tau = as.numeric(tau_value), Emptau = as.numeric(emptau_value), Pvalue = as.numeric(p_value), AIC = as.numeric(aic_value))
  
  # Write the data frame to a CSV file
  write.csv(result_df, bestcopula_path, row.names = FALSE, quote = FALSE)
  
  print(bestcopula_path)
  
}