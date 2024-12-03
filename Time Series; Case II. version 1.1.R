library(ggplot2)
library(rstan)
library(vars)
library(gridExtra)

# 1. DATA EXPLORATION
# TASK 1.1: LOADING AND EXPLORING DATA
FRED_data <- read.csv("C:/Users/flipv/OneDrive/Documents/Time series/case/Fred_Data.csv")
FRED_data$sasdate <- as.Date(FRED_data$sasdate, format = "%m/%d/%Y") # Converting the sasdate column
Data1 <- subset(FRED_data, sasdate <= as.Date("2019-12-31")) # Excluding pandemic data

# The appendix gives the following description and transformation codes (tcode) for the variables:
# INDPRO    => IP Index, 5 (means: \Delta log(x_t))
# CPIAUCSL  => CPI: All  (seasonally adjusted), 6 (means: \Delta^2 log(x_t))
# FEDFUNDS  => Effective Federal Funds Rate, 2 (means: \Delta t_x)
Data <- data.frame(Date = Data1$sasdate, INDPRO = Data1$INDPRO, CPIAUCSL = Data1$CPIAUCSL, FEDFUNDS = Data1$FEDFUNDS)

# Transform to Numeric Values (not characters, as this would mess up the plot)
Data$INDPRO <- as.numeric(ifelse(Data$INDPRO == "", NA, Data$INDPRO))
Data$CPIAUCSL <- as.numeric(ifelse(Data$CPIAUCSL == "", NA, Data$CPIAUCSL))
Data$FEDFUNDS <- as.numeric(ifelse(Data$FEDFUNDS == "", NA, Data$FEDFUNDS))

ggplot(Data, aes(x = Date, y = INDPRO)) + geom_line(color = "blue") + labs(title = "Industrial Production Over Time", x = "Date", y = "Industrial Production Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = CPIAUCSL)) + geom_line(color = "blue") + labs(title = "Consumer Price Index Over Time", x = "Date", y = "Consumer Price Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = FEDFUNDS)) + geom_line(color = "blue") + labs(title = "Federal Funds Over Time", x = "Date", y = "Effective Federal Funds Rate") + theme_minimal()

# TASK 1.2: MAKING THE THREE VARIABLES STATIONARY (according to the tcodes)
# For INDPRO we take the log differences:
Data$INDPRO_stationary <- c(NA, diff(log(Data$INDPRO), differences = 1)) #Growth Rate
# For CPIAUCSL we take the second log differences:
Data$CPIAUCSL_stationary <- c(NA, NA, diff(log(Data$CPIAUCSL), differences = 2))
# For FEDFUNDS we take the first differences:
Data$FEDFUNDS_stationary <- c(NA, diff(Data$FEDFUNDS, differences = 1))

ggplot(Data, aes(x = Date, y = INDPRO_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Industrial Production Over Time", x = "Date", y = "Transformed Industrial Production Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = CPIAUCSL_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Consumer Price Index Over Time", x = "Date", y = "Transformed Consumer Price Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = FEDFUNDS_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Federal Funds Over Time", x = "Date", y = "Transformed Federal Funds Rate") + theme_minimal()

# 2. ESTIMATING VAR MODELS AND SELECTING THE LAG ORDER
# TASK 2.1: SPECIFY AND TRIVIATE VAR MODEL
VAR_data <- na.omit(Data[, c("INDPRO_stationary", "CPIAUCSL_stationary", "FEDFUNDS_stationary")])

# We have VAR(p) model, where p = lag order, and:
# y_{1,t} : INDPRO_stationary
# y_{2,t} : CPIAUCSL_stationary
# y_{3,t} : FEDFUNDS_stationary
# y_{1,t} = c_1 + sum_{i=1}^p (\phi_{11}INPRO_{t-i} + \phi_{12}CPIAUCSL_{t-i} + \phi_{13}FEDFUNDS_{t-i} + e_{1t})
# y_{2,t} = c_2 + sum_{i=1}^p (\phi_{21}INPRO_{t-i} + \phi_{22}CPIAUCSL_{t-i} + \phi_{23}FEDFUNDS_{t-i} + e_{2t})
# y_{3,t} = c_3 + sum_{i=1}^p (\phi_{31}INPRO_{t-i} + \phi_{32}CPIAUCSL_{t-i} + \phi_{33}FEDFUNDS_{t-i} + e_{3t})


# TASK 2.2: DETERMINING LAG ORDER

#max 10 lags
lag_orders <- 1:10

# Table for AIC and BIC values for lags from 1 to 10
# Pre-allocate vectors for AIC and BIC values
aic_values <- rep(NA, length(lag_orders)) # AIC = -2ln(L) + 2k
bic_values <- rep(NA, length(lag_orders)) # BIC = -2ln(L) + ln(n)k => n <- #observations-lags

# Loop over lag orders and compute AIC and BIC

for (p in lag_orders) {
  var_model <- VAR(VAR_data, p = p, type = "const")
  aic_values[p] <- AIC(var_model)
  bic_values[p] <- BIC(var_model)
}

# Combine results into a data frame
model_selection <- data.frame(
  Lag_Order = lag_orders,
  AIC = aic_values,
  BIC = bic_values
)

# Print the model selection table
print(model_selection)

# According to BIC lags = 2 
VAR_model <- VAR(VAR_data, p = 2, type = "const")
summary(VAR_model)

# 
# TASK 2.3: VALIDATING VAR MODEL
#STATIONARITY OF VARIABLES = good

#NO SERIAL CORRELATION IN THE RESIDUALS
#PORTMANTEAU
serial_test <- serial.test(VAR_model, lags.pt = 12, type = "PT.asymptotic") # test residual autocorrelation
print(serial_test)
#BRUSCH-GODFREY
serial.test(var_model, lags.bg = 12, type = "BG")

#We reject H0 --> Serial correlation of residuals

#NORMALITY OF RESIDUALS
normality.test(var_model)

#We reject H0 --> NON-normality of residuals

#STABILITY OF VAR
stability_test <- stability(VAR_model)
plot(stability_test) 

stability <- roots(var_model)
print(stability)

# All eigenvalues lie within unit circle --> Stable VAR model

#To fix serial correlation and non-normality of residuals:

#Instead of using 10 lags, we use p = 2 lags:
VAR_model <- VAR(VAR_data, p = 2, type = "const")

#Now, see if assumptions hold for these lags

#STATIONARITY OF VARIABLES = good

#--> Add unit-root TEST

#NO SERIAL CORRELATION IN THE RESIDUALS
#PORTMANTEAU
serial_test <- serial.test(VAR_model, lags.pt = 2, type = "PT.asymptotic")
print(serial_test)
#BRUSCH-GODFREY
bg_test <- serial.test(VAR_model, lags.bg = 2, type = "BG")
print(bg_test)

#We reject H0 --> Serial correlation of residuals

plot(residuals(VAR_model))
hist(residuals(VAR_model))

qqnorm(residuals(VAR_model))
qqline(residuals(VAR_model))


#NORMALITY OF RESIDUALS
normality_test <- normality.test(VAR_model)
print(normality_test)

#We reject H0 --> NON-normality of residuals

#STABILITY OF VAR
stability_test <- stability(VAR_model)
plot(stability_test)

stability <- roots(VAR_model)
print(stability)

#--> Stability satisfied

# 3. GRANGER CAUSALITY
# TASK 3.1: EXPLAIN GRANGER CAUSALITY
# - X Granger causes Y: Y can be better explained/predicted when we take X into account. 
# - Past values of X contain information that helps to predict Y beyond what can be 
#   explained by past values of Y alone.  
# --> X provides additional predictive power for Y

# = Not direct causal relationship

#Example: interest rate Granger causes inflation

# TASK 3.2: PERFORMING CAUSALITY TESTS & INTERPRETING RESULTS
# H_0: The variable does not Granger-cause the other variable(s) in the model.
# H_1: The variable Granger-causes at least one other variable in the model.

#This tests granger-causality for ALL OTHER variables jointly in model

causality(VAR_model, cause = "INDPRO_stationary")
#--> INDPRO granger-causes at least one of the 2 other variables
causality(VAR_model, cause = "CPIAUCSL_stationary")
#--> CPI does not granger-cause any of the 2 other variables
causality(VAR_model, cause = "FEDFUNDS_stationary")
#--> FEDFUNDS granger-causes at least one of the 2 other variables

# TASK 3.3: 6 Bivariate Granger-causality tests

# Define the data for the three variables
bivariate_data <- VAR_data[, c("INDPRO_stationary", "CPIAUCSL_stationary", "FEDFUNDS_stationary")]

# Define all pairs of variables
variable_pairs <- list(
  c("INDPRO_stationary", "CPIAUCSL_stationary"),
  c("INDPRO_stationary", "FEDFUNDS_stationary"),
  c("CPIAUCSL_stationary", "FEDFUNDS_stationary")
)

# Loop over each pair of variables
for (pair in variable_pairs) {
  var1 <- pair[1]
  var2 <- pair[2]
  
  # Subset the data for the two variables
  subset_data <- bivariate_data[, c(var1, var2)]
  
  # Estimate the bivariate VAR model
  bivariate_VAR <- VAR(subset_data, p = 2, type = "const")
  
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
# FEDFUNDS Granger causes INDPRO
# FEDFUNDS Granger causes CPI

#4 Impulse Response functions for reduced-form VAR's

# TASK 4.1: explaining IRF's
#Impulse Response Functions (IRFs) show how a 1-time unit shock to 1 variable affects itself and all other 
#variables in the VAR system 

#--> Shows whether response to shock are positive/negative
#--> Traces the path of the effect over time: how large is initial response?, 
#how long does the effect last?, does the shock dissapear /stabilize?


# TASK 4.2: Compute IRF's 
VAR_model <- VAR(VAR_data, p = 2, type = "const")

# IRFs for INDPRO_stationary
irf_indpro_to_indpro <- irf(VAR_model, impulse = "INDPRO_stationary", response = "INDPRO_stationary", n.ahead = 12, boot = TRUE, ci = 0.95)
irf_indpro_to_cpi <- irf(VAR_model, impulse = "INDPRO_stationary", response = "CPIAUCSL_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_indpro_to_ffr <- irf(VAR_model, impulse = "INDPRO_stationary", response = "FEDFUNDS_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)

# IRFs for CPIAUCSL_stationary
irf_cpi_to_indpro <- irf(VAR_model, impulse = "CPIAUCSL_stationary", response = "INDPRO_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_cpi_to_cpi <- irf(VAR_model, impulse = "CPIAUCSL_stationary", response = "CPIAUCSL_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_cpi_to_ffr <- irf(VAR_model, impulse = "CPIAUCSL_stationary", response = "FEDFUNDS_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)

# IRFs for FEDFUNDS_stationary
irf_ffr_to_indpro <- irf(VAR_model, impulse = "FEDFUNDS_stationary", response = "INDPRO_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_ffr_to_cpi <- irf(VAR_model, impulse = "FEDFUNDS_stationary", response = "CPIAUCSL_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)
irf_ffr_to_ffr <- irf(VAR_model, impulse = "FEDFUNDS_stationary", response = "FEDFUNDS_stationary", n.ahead = 24, boot = TRUE, ci = 0.95)

# Plot each IRF step by step
# INDPRO_stationary Impulse
plot(irf_indpro_to_indpro, main = "Response of INDPRO to a Shock in INDPRO")
plot(irf_indpro_to_cpi, main = "Response of CPI to a Shock in INDPRO")
plot(irf_indpro_to_ffr, main = "Response of FFR to a Shock in INDPRO")

# CPIAUCSL_stationary Impulse
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

#Reduced-form VAR = original VAR(2) model

#-->Each variable explained by lags of its own and other variables
#-->We measure correlations between variables but dont know why the correlation happens
#The errors are correlated because they capture combined effects of multiple unobserved factors

#Structural VAR: we go beyond correlation and try to identify causal relationships by imposing assumptions

#--> Goal is to decompose the reduced-form residuals into independent structural shocks 
#--> These shocks are uncorrelated and have a clear interpretation 
# enables you to use economic theory to set restrictions
#Reduced-form VAR shows what happens, SVAR shows why it happens

# TASK 5.2: SVAR Model Equation
library(svars)

# Step 1: Estimate the reduced-form VAR model
VAR_model <- VAR(VAR_data, p = 1, type = "const")

# Step 2: Specify the short-term restrictions (economic reasoning)
B_matrix <- matrix(c(
  1,   0,   0,  # Inflation (CPIAUCSL) depends only on its own shock
  NA,   1,   0,  # Interest Rate (FEDFUNDS) depends on inflation
  NA,  NA,   1   # Industrial Production (INDPRO) depends on inflation and interest rates
), nrow = 3, byrow = TRUE)

# Step 3: Estimate the SVAR model
SVAR_model <- id.st(VAR_model, restriction_matrix = B_matrix)

# Step 4: Summarize the SVAR model
summary(SVAR_model)


# Compute IRFs
irf_svar <- irf(SVAR_model, n.ahead = 24, boot = TRUE, ci = 0.95)

# Plot IRFs
plot(irf_svar)
# Compute variance decomposition
fevd_svar <- fevd(SVAR_model, n.ahead = 24)

# Plot variance decomposition
plot(fevd_svar)


library(svars)

# Cholesky decomposition
SVAR_model <- id.chol(VAR_model)

# Summarize the SVAR model
summary(SVAR_model)

library(vars)
library(svars)

# Use Cholesky decomposition for SVAR identification
SVAR_model <- id.chol(VAR_model)

# Compute IRFs
# Specify the horizon (e.g., 24 steps ahead) and bootstrap for confidence intervals
irf_results_svar <- irf.svars(SVAR_model, n.ahead = 24, boot = TRUE, ci = 0.95)


# Plot the response of each variable to a shock in INDPRO_stationary
plot(irf_results_svar, impulse = "INDPRO_stationary", response = "INDPRO_stationary", main = "Response of INDPRO to INDPRO Shock")
plot(irf_results_svar, impulse = "INDPRO_stationary", response = "CPIAUCSL_stationary", main = "Response of CPIAUCSL to INDPRO Shock")
plot(irf_results_svar, impulse = "INDPRO_stationary", response = "FEDFUNDS_stationary", main = "Response of FEDFUNDS to INDPRO Shock")

# Plot the response of each variable to a shock in CPIAUCSL_stationary
plot(irf_results_svar, impulse = "CPIAUCSL_stationary", response = "INDPRO_stationary", main = "Response of INDPRO to CPIAUCSL Shock")
plot(irf_results_svar, impulse = "CPIAUCSL_stationary", response = "CPIAUCSL_stationary", main = "Response of CPIAUCSL to CPIAUCSL Shock")
plot(irf_results_svar, impulse = "CPIAUCSL_stationary", response = "FEDFUNDS_stationary", main = "Response of FEDFUNDS to CPIAUCSL Shock")

# Plot the response of each variable to a shock in FEDFUNDS_stationary
plot(irf_results_svar, impulse = "FEDFUNDS_stationary", response = "INDPRO_stationary", main = "Response of INDPRO to FEDFUNDS Shock")
plot(irf_results_svar, impulse = "FEDFUNDS_stationary", response = "CPIAUCSL_stationary", main = "Response of CPIAUCSL to FEDFUNDS Shock")
plot(irf_results_svar, impulse = "FEDFUNDS_stationary", response = "FEDFUNDS_stationary", main = "Response of FEDFUNDS to FEDFUNDS Shock")

print(SVAR_model$B)

# Variance decomposition
fevd_svar <- fevd(SVAR_model, n.ahead = 24)

# Plot variance decomposition
plot(fevd_svar)


library(svars)

# Define the B-matrix based on economic reasoning
B_matrix <- matrix(c(
  1,   0,   0,  # Inflation reacts only to its own shock
  NA,   1,  0,  # Interest rates react to inflation
  NA,  NA,   1   # Production reacts to inflation and interest rates
), nrow = 3, byrow = TRUE)

# Apply short-term restrictions to identify the SVAR
SVAR_model <- id.st(VAR_model, restriction_matrix = B_matrix)

# Compute new IRFs
irf_svar <- irf(SVAR_model, n.ahead = 24, boot = TRUE, ci = 0.95)

# Plot IRFs
plot(irf_svar, impulse = "CPIAUCSL_stationary", response = "INDPRO_stationary", main = "Response of INDPRO to CPIAUCSL Shock")


residuals_VAR <- residuals(VAR_model)
cov_matrix <- cov(residuals_VAR)
print(cov_matrix)
det(cov_matrix)  # Should be non-zero

SVAR_model <- id.chol(VAR_model)
summary(SVAR_model)

cor_matrix <- cor(na.omit(VAR_data))
print(cor_matrix)
eigenvalues <- eigen(cov(residuals(VAR_model)))$values
print(eigenvalues)

#testkajlsdhfkjahsdg