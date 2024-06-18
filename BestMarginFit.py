# -*- coding: utf-8 -*-

# Author: Ziying Zhou
# Date: June 8, 2024
# Description: This script calculates the optimal marginal distribution for wind speed, rainfall, and storm surge.

# Initialization Settings =====================================================

import pandas as pd
from scipy.stats import kstest
from scipy.stats import burr, cauchy, expon, gamma, genextreme, genpareto, gumbel_l, logistic, lognorm, norm, weibull_max
import os

# Define the root directory and subdirectories for data
root = r"G:\Project_CompoundHazard\Data"
copula_dir = os.path.join(root, "Copula")
record_dir = os.path.join(copula_dir, "Record")
bestmargin_dir = os.path.join(copula_dir, "BestMargin")
probmargin_dir = os.path.join(copula_dir, "ProbMargin")

# Path to the station list
stationlist_path = r"StationList.csv"

# List of variables and distributions to evaluate
listVariable = ["Wind", "Rain", "Surge"]
listDistribution = ["burr", "cauchy", "expon", "gamma", "genextreme",
                    "genpareto", "gumbel_l", "logistic", "lognorm", "norm",
                    "weibull_max"]

# Define the marginal distribution functions
distributions_margin = {
    "burr": burr,    
    "cauchy": cauchy,
    "expon": expon,
    "gamma": gamma,
    "genextreme": genextreme,
    "genpareto": genpareto, 
    "gumbel_l": gumbel_l,
    "logistic": logistic,
    "lognorm": lognorm,  
    "norm": norm,
    "weibull_max": weibull_max
}

# Marginal Distribution Fitting ===============================================

# Read the station list
dfStation = pd.read_csv(stationlist_path)

# Loop through each station
for i in range(len(dfStation)):
  
    dfTemp = dfStation.iloc[i]  
    
    station_id = str(int(dfTemp["ID"] + 100))[-2:]
 
    record_path = os.path.join(record_dir, f"Record{station_id}.csv")
    bestmargin_path = os.path.join(bestmargin_dir, f"BestMargin{station_id}.csv")

    dfRecord = pd.read_csv(record_path)
       
    listBestMargin = []
       
    for column in ["Wind", "Rain", "Surge"]:
    
        column_data = dfRecord[column].dropna()
    
        for dist_name, dist_func in distributions_margin.items():
            
            try:
                # Fit the distribution to the data
                params = dist_func.fit(column_data)
                # Perform the KS test
                D, p_value = kstest(column_data, dist_name, args=params)
                
                listBestMargin.append({
                    "Variable": column,
                    "Distribution": dist_name,
                    "KS_D": D,
                    "KS_P": p_value,
                    "Parameters": params
                })
    
            except Exception as e:
                print(f"Error fitting {dist_name} for {column}: {e}")
    
    # Save the best margin results to a CSV file
    dfBestMargin = pd.DataFrame(listBestMargin)
    dfBestMargin.to_csv(bestmargin_path, index=False)
    
    print(bestmargin_path)
      
# Optimal Distribution Cumulative Distribution Function =======================

# Define the CDF functions for the distributions
distributions_margin_cdf = {
    "burr": burr.cdf,    
    "cauchy": cauchy.cdf,
    "expon": expon.cdf,
    "gamma": gamma.cdf,
    "genextreme": genextreme.cdf,
    "genpareto": genpareto.cdf, 
    "gumbel_l": gumbel_l.cdf,
    "logistic": logistic.cdf,
    "lognorm": lognorm.cdf,  
    "norm": norm.cdf,
    "weibull_max": weibull_max.cdf
}

# Read the station list again
dfStation = pd.read_csv(stationlist_path)

# Loop through each station
for i in range(len(dfStation)):
    
    dfTemp = dfStation.iloc[i]  

    station_id = str(int(dfTemp["ID"] + 100))[-2:]
    
    record_path = os.path.join(record_dir, f"Record{station_id}.csv")
    bestmargin_path = os.path.join(bestmargin_dir, f"BestMargin{station_id}.csv")
    probmargin_path = os.path.join(probmargin_dir, f"ProbMargin{station_id}.csv")

    dfRecord = pd.read_csv(record_path)
    dfBestMargin = pd.read_csv(bestmargin_path)

    listProbMargin = {}
        
    for variable in listVariable:
        
        ks_min = 1
        dist_best = ""
        
        for dist_info in dfBestMargin.itertuples():
            
            if dist_info.Variable != variable:
                continue
            
            dist_name = dist_info.Distribution 
            ks = dist_info.KS_D 
            params = eval(dist_info.Parameters)  
            
            # Calculate the CDF based on the best distribution
            cdf = distributions_margin_cdf[dist_name](dfRecord[variable], *params)
        
            if ks < ks_min:          
                listProbMargin[variable] = cdf
                ks_min = ks
                dist_best = dist_name
        
        print(variable, dist_best)
    
    # Save the probability margin results to a CSV file
    dfProbMargin = pd.DataFrame(listProbMargin)
    dfProbMargin.to_csv(probmargin_path, index=False)
    
    print(probmargin_path)
