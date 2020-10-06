# MultiGrangerNetwork
# Copyright 2017, Andrey Skripnikov, All rights reserved.
Code to the paper titled "Estimation of Multi-Granger Network Causal Models" by A. Skripnikov

####################
##  Paper.functions.R
####################

Contains all the main functions for data generation, ADMM algorithm, separate and joint estimation procedures and more. Code is very well-documented.


####################
##  Simulated.Data.Two.Entities.R
####################

Executive file for the simulation study for the case two entities from Section 5. You can manipulate the parameters(number of variables per entity, number of time points, number of factors in L-factor model etc). Whole description can be found in the comments in the code.

####################
##  Real.Data.Four.States.R
####################

Executive file for the econometric time series data study from Section 6. Can pick one of four settings and manipulate multiple parameters, description can be found in the comments in the code.


####################
### PA.csv, IL.csv, OH.csv, MI.csv
####################

Files containing Federal Reserve (FRED) seasonally-adjusted monthly time series data - from https://fred.stlouisfed.org/tags/series -  on 18 economic indicator variables (all described in the paper) from 1976 up until 2016 for each state. Function 'Load.Data', defined in 'Paper.functions.R' file and executed in 'Real.Data.Four.States.R' file, deals with extracting our time series of interest out of these files.

