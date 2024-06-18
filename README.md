# ComHazAsTC-RRE: Compound Hazard Assessment of Tropical Cyclones within RRE Framework

## Overview

This repository contains the source code and processed data used in the study titled "ComHazAsTC-RRE: Compound Hazard Assessment of Tropical Cyclones within RRE Framework". The study develops and applies the ComHazAsTC-RRE model to assess the compound hazards of wind, rainfall, and storm surge induced by tropical cyclones along China's coast.

## Contents

- `StationList.csv`: Contains basic information of coastal stations from the GTSM and GSOD datasets.
- `Record/`: Directory with records of wind, rain, and surge at coastal stations in China during tropical cyclone events.
- `BestMarginFit.py`: Python script for calculating the optimal marginal distribution for wind speed, rainfall, and storm surge.
- `C-VineCopula.R`: R script for calculating the optimal C-Vine Copula function and its corresponding parameters.
- `HazardAssessment.R`: R script for calculating joint probabilities, joint return periods, the Compound Hazard Index, and failure probabilities using C-Vine Copula functions.
- `README.md`: This readme file.

## Requirements

The following Python packages are required to run the scripts:

- pandas
- scipy

The following R packages are required to run the scripts:

- VineCopula
- CDVineCopulaConditional

## Data Availability 

The data used in this study are sourced from publicly accessible datasets:

- **[Global Surface Summary of the Day (GSOD)](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00516)**: Provides long-term time series of daily maximum wind speed and cumulative rainfall (24h) on a global scale.
- **[Global Tide and Surge Model (GTSM)](https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-water-level-change-timeseries-cmip6)**: Provides long-term time series of daily maximum storm surge on a global scale.

The processed data files are stored in the `Record/` directory.
