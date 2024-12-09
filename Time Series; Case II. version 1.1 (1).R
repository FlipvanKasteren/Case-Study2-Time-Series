library(ggplot2)
library(rstan)
library(vars)
library(gridExtra)
library(bootUR)
library(svars)
        
# 1. DATA EXPLORATION
# TASK 1.1: LOADING AND EXPLORING DATA
FRED_data <- read.csv('DataCase.csv')
FRED_data$sasdate <- as.Date(FRED_data$sasdate, format = "%m/%d/%Y") # Converting the sasdate column
Data1 <- subset(FRED_data, sasdate <= as.Date("2019-12-31")) # Excluding pandemic data

# The appendix gives the following description and transformation codes (tcode) for the variables:
# INDPRO    => IP Index, 5 (means: \Delta log(x_t))
# CPIULFSL  => CPI: All except food  (seasonally adjusted), 6 (means: \Delta^2 log(x_t))
# FEDFUNDS  => Effective Federal Funds Rate, 2 (means: \Delta t_x)
Data <- data.frame(Date = Data1$sasdate, INDPRO = Data1$INDPRO, CPIULFSL = Data1$CPIULFSL, FEDFUNDS = Data1$FEDFUNDS)

# Transform to Numeric Values (not characters, as this would mess up the plot)
Data$INDPRO <- as.numeric(ifelse(Data$INDPRO == "", NA, Data$INDPRO))
Data$CPIULFSL <- as.numeric(ifelse(Data$CPIULFSL == "", NA, Data$CPIULFSL))
Data$FEDFUNDS <- as.numeric(ifelse(Data$FEDFUNDS == "", NA, Data$FEDFUNDS))

ggplot(Data, aes(x = Date, y = INDPRO)) + geom_line(color = "blue") + labs(title = "Industrial Production Over Time", x = "Date", y = "Industrial Production Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = CPIULFSL)) + geom_line(color = "blue") + labs(title = "Consumer Price Index Over Time", x = "Date", y = "Consumer Price Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = FEDFUNDS)) + geom_line(color = "blue") + labs(title = "Federal Funds Over Time", x = "Date", y = "Effective Federal Funds Rate") + theme_minimal()

# TASK 1.2: MAKING THE THREE VARIABLES STATIONARY (according to the tcodes)
# For INDPRO we take the log differences:
Data$INDPRO_stationary <- c(NA, diff(log(Data$INDPRO), differences = 1)) #Growth Rate
# For CPIULFSL we take the second log differences:
Data$CPIULFSL_stationary <- c(NA, NA, diff(log(Data$CPIULFSL), differences = 2))
# For FEDFUNDS we take the first differences:
Data$FEDFUNDS_stationary <- c(NA, diff(Data$FEDFUNDS, differences = 1))

ggplot(Data, aes(x = Date, y = INDPRO_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Industrial Production Over Time", x = "Date", y = "Transformed Industrial Production Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = CPIULFSL_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Consumer Price Index Over Time", x = "Date", y = "Transformed Consumer Price Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = FEDFUNDS_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Federal Funds Over Time", x = "Date", y = "Transformed Federal Funds Rate") + theme_minimal()


# 2. ESTIMATING VAR MODELS AND SELECTING THE LAG ORDER
# TASK 2.1: SPECIFY AND TRIVIATE VAR MODEL
VAR_data <- na.omit(Data[, c("INDPRO_stationary", "CPIULFSL_stationary", "FEDFUNDS_stationary")])


## check if all time series are actually stationary
INDPRO_stationary <- VAR_data$INDPRO_stationary
CPIULFSL_stationary <- VAR_data$CPIULFSL_stationary
FEDFUNDS_stationary <- VAR_data$FEDFUNDS_stationary

# ADF test for INDPRO_stationary
adf_indpro <- adf(INDPRO_stationary)

# ADF test for CPIULFSL_stationary
adf_cpi <- adf(CPIULFSL_stationary)

# ADF test for FEDFUNDS_stationary
adf_fedfunds <- adf(FEDFUNDS_stationary)

print(adf_indpro)
print(adf_cpi)
print(adf_fedfunds)

# We have VAR(p) model, where p = lag order, and:
# y_{1,t} : INDPRO_stationary
# y_{2,t} : CPIULFSL_stationary
# y_{3,t} : FEDFUNDS_stationary
# y_{1,t} = c_1 + sum_{i=1}^p (\phi_{11}INPRO_{t-i} + \phi_{12}CPIULFSL_{t-i} + \phi_{13}FEDFUNDS_{t-i} + e_{1t})
# y_{2,t} = c_2 + sum_{i=1}^p (\phi_{21}INPRO_{t-i} + \phi_{22}CPIULFSL_{t-i} + \phi_{23}FEDFUNDS_{t-i} + e_{2t})
# y_{3,t} = c_3 + sum_{i=1}^p (\phi_{31}INPRO_{t-i} + \phi_{32}CPIULFSL_{t-i} + \phi_{33}FEDFUNDS_{t-i} + e_{3t})

# TASK 2.2: DETERMINING LAG ORDER
#POPULAR METHODs for LAG SELECTION
lag_selection <- VARselect(VAR_data, lag.max = 20, type = "const")
print(lag_selection)
# The 4 different criteria suggest the following number of lags (p):
# p <- 15 : AIC(n)
# p <- 5  : HQ(n)
# p <- 3 : SC(n)/BIC(n)
# p <- 15 : FPE(n) --> Final Prediction Error (find p that min FPE)

# Sequential tests on the nullity of coefficient matrices.
# Fit VAR models for different lags
max_lag <- 20  # Maximum lag length
var_models <- lapply(1:max_lag, function(p) VAR(VAR_data, p = p, type = "const"))

# Perform sequential likelihood ratio (LR) tests
p_values <- c()  # To store p-values
for (p in max_lag:2) {
  # Current model and previous model
  model_p <- var_models[[p]]
  model_p_minus_1 <- var_models[[p - 1]]
  
  # Extract log-likelihoods
  logL_p <- logLik(model_p)
  logL_p_minus_1 <- logLik(model_p_minus_1)
  
  # Calculate LR statistic
  LR_stat <- -2 * (as.numeric(logL_p_minus_1) - as.numeric(logL_p))
  
  # Degrees of freedom: number of parameters in one lag
  n_vars <- ncol(VAR_data)  # Number of variables in the VAR
  df <- n_vars^2
  
  # Calculate p-value
  p_value <- 1 - pchisq(LR_stat, df)
  p_values <- c(p_values, p_value)
}

# Combine results into a data frame
results <- data.frame(Lag = max_lag:2, P_Value = p_values)
print(results)

# Determine the optimal lag (where p-value > 0.05 for the first time)
optimal_lag <- min(results$Lag[results$P_Value > 0.05])
cat("Optimal lag selected:", optimal_lag, "\n")

## We decided to choose the BIC for the lag length selection, not the sequential test. So lag length = 3
#ESTIMATING VAR(p=3) MODEL
VAR_model <- VAR(VAR_data, p = 3, type = "const")
summary(VAR_model)

# TASK 2.3: VALIDATING VAR MODEL
# Check for stability of the VAR model
stability <- roots(VAR_model)
print(stability)

# Are all roots inside the unit circle?
if (all(Mod(stability) < 1)) {
  cat("The VAR model is stable.\n")
} else {
  cat("The VAR model is NOT stable.\n")
}

# Extract residuals
var_residuals <- residuals(VAR_model)

# Perform Ljung-Box test on residuals of each equation
for (i in colnames(var_residuals)) {
  cat("Ljung-Box Test for residuals of:", i, "\n")
  print(Box.test(var_residuals[, i], lag = 10, type = "Ljung-Box"))
  cat("\n")
}

# Ljung-box test suggests autocorrelation in residuals of INDPRO

# ACF for INDPRO_stationary residuals
acf(var_residuals[, "INDPRO_stationary"], main = "ACF of INDPRO Residuals", lag.max = 20)

# ACF for CPIULFSL_stationary residuals
acf(var_residuals[, "CPIULFSL_stationary"], main = "ACF of CPIULFSL Residuals", lag.max = 20)

# ACF for FEDFUNDS_stationary residuals
acf(var_residuals[, "FEDFUNDS_stationary"], main = "ACF of FEDFUNDS Residuals", lag.max = 20)

# Test for cross-correlations

# Calculate the covariance matrix of residuals
cov_matrix <- cov(var_residuals) 
print(cov_matrix)

# Cross-correlation between INDPRO and CPIULFSL residuals
ccf(var_residuals[, "INDPRO_stationary"], var_residuals[, "CPIULFSL_stationary"],
    main = "Cross-correlation between INDPRO and CPIULFSL residuals")

# Cross-correlation between INDPRO and FEDFUNDS residuals
ccf(var_residuals[, "INDPRO_stationary"], var_residuals[, "FEDFUNDS_stationary"],
    main = "Cross-correlation between INDPRO and FEDFUNDS residuals")

# Cross-correlation between CPIULFSL and FEDFUNDS residuals
ccf(var_residuals[, "CPIULFSL_stationary"], var_residuals[, "FEDFUNDS_stationary"],
    main = "Cross-correlation between CPIULFSL and FEDFUNDS residuals")

# 3. GRANGER CAUSALITY
# TASK 3.1: EXPLAIN GRANGER CAUSALITY
# - X Granger causes Y: Y can be better explained/predicted when we take X into account. 
# - Past values of X contain information that helps to predict Y beyond what can be 
#   explained by past values of Y alone.  
# --> X provides additional predictive power for Y

# = Not direct causal relationship

# TASK 3.2: PERFORMING CAUSALITY TESTS & INTERPRETING RESULTS
# H_0: The variable does not Granger-cause the other variable(s) in the model.
# H_1: The variable Granger-causes at least one other variable in the model.

#This tests granger-causality for ALL OTHER variables jointly in model

causality(VAR_model, cause = "INDPRO_stationary")
#--> INDPRO granger-causes at least one of the 2 other variables
causality(VAR_model, cause = "CPIULFSL_stationary")
#--> CPI does not granger-cause any of the 2 other variables
causality(VAR_model, cause = "FEDFUNDS_stationary")
#--> FEDFUNDS granger-causes at least one of the 2 other variables

# TASK 3.3: 6 Bivariate Granger-causality tests

# Define the data for the three variables
bivariate_data <- VAR_data[, c("INDPRO_stationary", "CPIULFSL_stationary", "FEDFUNDS_stationary")]

# Define all pairs of variables
variable_pairs <- list(
  c("INDPRO_stationary", "CPIULFSL_stationary"),
  c("INDPRO_stationary", "FEDFUNDS_stationary"),
  c("CPIULFSL_stationary", "FEDFUNDS_stationary")
)

# Loop over each pair of variables
for (pair in variable_pairs) {
  var1 <- pair[1]
  var2 <- pair[2]
  
  # Subset the data for the two variables
  subset_data <- bivariate_data[, c(var1, var2)]
  
  # Estimate the bivariate VAR model
  bivariate_VAR <- VAR(subset_data, p = 3, type = "const")
  
  # Perform Granger causality tests
  cat(paste0("\n--- Granger Causality Tests for: ", var1, " and ", var2, " ---\n"))
  
  # Test if var1 Granger-causes var2
  gc_test_1 <- causality(bivariate_VAR, cause = var1)
  print(paste0(var1, " Granger-causes ", var2, ":"))
  print(gc_test_1)
  
  # Test if var2 Granger-causes var1
  gc_test_2 <- causality(bivariate_VAR, cause = var2)
  print(paste0(var2, " Granger-causes ", var1, ":"))
  print(gc_test_2)
}

# RESULTS
# INDPRO Granger causes FEDFUNDS
# FEDFUNDS Granger causes CPI


#4 Impulse Response functions for reduced-form VAR's

# TASK 4.1: explaining IRF's
#Impulse Response Functions (IRFs) show how a 1-time unit shock to 1 variable affects itself and all other 
#variables in the VAR system 

#--> Shows whether response to shock are positive/negative
#--> Traces the path of the effect over time: how large is initial response?, 
#how long does the effect last?, does the shock dissapear /stabilize?


# IRFs for INDPRO_stationary
irf_indpro_to_indpro <- irf(VAR_model, impulse = "INDPRO_stationary", response = "INDPRO_stationary", n.ahead = 12, boot = TRUE, ci = 0.95)
irf_indpro_to_cpi <- irf(VAR_model, impulse = "INDPRO_stationary", response = "CPIULFSL_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_indpro_to_ffr <- irf(VAR_model, impulse = "INDPRO_stationary", response = "FEDFUNDS_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)

# IRFs for CPIULFSL_stationary
irf_cpi_to_indpro <- irf(VAR_model, impulse = "CPIULFSL_stationary", response = "INDPRO_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_cpi_to_cpi <- irf(VAR_model, impulse = "CPIULFSL_stationary", response = "CPIULFSL_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_cpi_to_ffr <- irf(VAR_model, impulse = "CPIULFSL_stationary", response = "FEDFUNDS_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)

# IRFs for FEDFUNDS_stationary
irf_ffr_to_indpro <- irf(VAR_model, impulse = "FEDFUNDS_stationary", response = "INDPRO_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_ffr_to_cpi <- irf(VAR_model, impulse = "FEDFUNDS_stationary", response = "CPIULFSL_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_ffr_to_ffr <- irf(VAR_model, impulse = "FEDFUNDS_stationary", response = "FEDFUNDS_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)

# Plot each IRF step by step
# INDPRO_stationary Impulse
plot(irf_indpro_to_indpro, main = "Response of INDPRO to a Shock in INDPRO")
plot(irf_indpro_to_cpi, main = "Response of CPI to a Shock in INDPRO")
plot(irf_indpro_to_ffr, main = "Response of FFR to a Shock in INDPRO")

# CPIULFSL_stationary Impulse
plot(irf_cpi_to_indpro, main = "Response of INDPRO to a Shock in CPI")
plot(irf_cpi_to_cpi, main = "Response of CPI to a Shock in CPI")
plot(irf_cpi_to_ffr, main = "Response of FFR to a Shock in CPI")

# FEDFUNDS_stationary Impulse
plot(irf_ffr_to_indpro, main = "Response of INDPRO to a Shock in FFR")
plot(irf_ffr_to_cpi, main = "Response of CPI to a Shock in FFR")
plot(irf_ffr_to_ffr, main = "Response of FFR to a Shock in FFR")

# TASK 4.3: Prize Puzzle

# Bernanke et al. (2005): The Prize Puzzle
# The price puzzle refers to the counter-intuitive result where an increase in monetary policy interest rates leads to an increase in the price level.

#--> OCCURS IN OUR ANALYSIS

# Relevance in Bernanke et al. (2005): This puzzle arises in certain VAR models when the model specification or data transformations 
# fail to account for all relevant information, particularly about inflation expectations.

# TASK 4.4: Confidence Intervals of IRF's

# Package computes the confidence intervals using a Bootstrapping procedure

#1. Estimate coefficients of VAR model to compute point estimates
#2. Obtain residuals from the fitted VAR equations
#3. Bootstrap: resample residuals to create new datasets
#4. For a large number of resamples, compute new IRF's using re-estimated VAR-coefficients
#5. After all bootstrap iterations, IRF estimates form an empirical distribution and CI's are derived as usual

# SIGNIFICANCE of confidence intervals?

# TASK 4.5: obtaining original IRF(in levels) from IRF corresponding to stationary series(in differences)

#Lutkepöhl (2005): to change IRF from differences to levels: compute cumulative sum of the IR's of the differences series
#for each time horizon --> IRF in levels

#For FEDFUNDS IRF (expressed in differences)
#--> Lutkepohl procedure

#For INDPRO IRF (expressed in log-differences):
#--> Adapt for log-differences

#For CPI (expressed in 2nd order log-differences):
#--> Adapt for 2nd order log-differences

#OBTAINING IRF(levels) from IRF(original)




                     


# TASK 4.6: 
# We have now obtained IRF's for the original time series 
# = IRF's for the "reduced-form" specification of the VAR

#--> Can we use these for structural analysis?

#NO: they are useful for understanding statistical relationships but not for causal
#.   interpretation or structural analysis

#WHY: the errors are correlated across equations 

#--> We require isolating these orthogonal shocks to identify the causal 
#.   impact of one variable on another


#--> Need to transform the reduced-form VAR into a Structural VAR

#5. Structural VARS (SVAR)

# TASK 5.1: Difference between reduced-form VAR and SVAR

#Reduced-form VAR = original VAR(3) model

#-->Each variable explained by lags of its own and other variables
#-->We measure correlations between variables but dont know why the correlation happens
#The errors are correlated because they capture combined effects of multiple unobserved factors

#Structural VAR: we go beyond correlation and try to identify causal relationships by imposing assumptions

#--> Goal is to decompose the reduced-form residuals into independent structural shocks 
#--> These shocks are uncorrelated and have a clear interpretation 
# enables you to use economic theory to set restrictions
#Reduced-form VAR shows what happens, SVAR shows why it happens

# TASK 5.2: SVAR Model Equation

# Ensure the variables are in the correct order
# Reorder the VAR data to match the recursive order: FEDFUNDS -> CPIULFSL -> INDPRO
VAR_data_reordered1 <- VAR_data[, c("CPIULFSL_stationary", "FEDFUNDS_stationary","INDPRO_stationary")]
VAR_data_reordered2 <- VAR_data[, c("INDPRO_stationary", "CPIULFSL_stationary","FEDFUNDS_stationary")]

# Fit the reduced-form VAR model with 3 lags (chosen previously)
VAR_model1 <- VAR(VAR_data_reordered1, p = 3, type = "const")

# Apply the Cholesky decomposition to identify the structural VAR
SVAR_model1 <- id.chol(VAR_model1)

# Display the identified structural VAR model
summary(SVAR_model1)

# Extract the structural shocks
structural_shocks <- residuals(SVAR_model1)
head(structural_shocks)

# Plot Impulse Response Functions (IRFs) for the Structural VAR
irf_svar <- irf(SVAR_model1, n.ahead = 12, boot = TRUE, ci = 0.95)

# Plot IRFs
plot(irf_svar, main = "Impulse Response Functions for SVAR (Cholesky)")

# Fit the reduced-form VAR model with 3 lags (chosen previously)
VAR_model2 <- VAR(VAR_data_reordered2, p = 3, type = "const")

# Apply the Cholesky decomposition to identify the structural VAR
SVAR_model2 <- id.chol(VAR_model2)

# Display the identified structural VAR model
summary(SVAR_model2)

# Extract the structural shocks
structural_shocks <- residuals(SVAR_model2)
head(structural_shocks)

# Plot Impulse Response Functions (IRFs) for the Structural VAR
irf_svar2 <- irf(SVAR_model2, n.ahead = 12, boot = TRUE, ci = 0.95)

# Plot IRFs
plot(irf_svar2, main = "Impulse Response Functions for SVAR (Cholesky)")
plot(irf_svar, main = "Impulse Response Functions for SVAR (Cholesky)")

names(irf_svar$irf)
names(irf_svar2$irf)
irf_svar2$irf <- irf_svar2$irf[names(irf_svar$irf)]
irf_svar2$Lower <- irf_svar2$Lower[names(irf_svar$irf)]
irf_svar2$Upper <- irf_svar2$Upper[names(irf_svar$irf)]
par(mfrow = c(1, 2))  # Set up side-by-side plots

# Plot Model 1
plot(irf_svar, main = "Impulse Response Functions for SVAR (Model 1)")

# Plot Model 2 (with reordered shocks)
plot(irf_svar2, main = "Impulse Response Functions for SVAR (Model 2)")
