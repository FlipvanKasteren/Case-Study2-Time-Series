library(ggplot2)
library(rstan)
library(vars)

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

ggplot(Data, aes(x = Date, y = INDPRO)) + geom_line(color = "blue") + labs(title = "Industrial Production Over Time", x = "Date", y = "Industrial Production Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = CPIAUCSL)) + geom_line(color = "blue") + labs(title = "Consumer Price Index Over Time", x = "Date", y = "Consumer Price Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = FEDFUNDS)) + geom_line(color = "blue") + labs(title = "Federal Funds Over Time", x = "Date", y = "Effective Federal Funds Rate") + theme_minimal()

# TASK 1.2: MAKING THE THREE VARIABLES STATIONARY (according to the tcodes)
# For INDPRO we take the log differences:
Data$INDPRO_stationary <- c(NA, diff(log(Data$INDPRO), differences = 1))
# For CPIAUCSL we take the second log differences:
Data$CPIAUCSL_stationary <- c(NA, NA, diff(log(Data$CPIAUCSL), differences = 2))
# For FEDFUNDS we take the first differences:
Data$FEDFUNDS_stationary <- c(NA, diff(Data$FEDFUNDS, differences = 1))

ggplot(Data, aes(x = Date, y = INDPRO_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Industrial Production Over Time", x = "Date", y = "Transformed Industrial Production Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = CPIAUCSL_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Consumer Price Index Over Time", x = "Date", y = "Transformed Consumer Price Index") + theme_minimal()
ggplot(Data, aes(x = Date, y = FEDFUNDS_stationary)) + geom_line(color = "blue") + labs(title = "Transformed Federal Funds Over Time", x = "Date", y = "Transformed Federal Funds Rate") + theme_minimal()

## TESTING FOR STATIONARITY ADF??


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

lag_selection <- VARselect(VAR_data, lag.max = 10, type = "const")
print(lag_selection)
# The 4 different criteria suggest the following number of lags (p):
# p <- 10 : AIC(n)
# p <- 9  : HQ(n)
# p <- 3  : SC(n)
# p <- 10 : FPE(n)

# TASK 2.2: DETERMINING LAG ORDER
VAR_model <- VAR(VAR_data, p = 10, type = "const")
summary(VAR_model)
AIC(VAR_model) # AIC = -2ln(L) + 2k     => L <- max l.h. and k <- #parameters
BIC(VAR_model) # BIC = -2ln(L) + ln(n)k => n <- #observations

# Define the range of lag orders
lag_orders <- 1:10

# Initialize vectors to store AIC and BIC values
aic_values <- numeric(length(lag_orders))  # Correct initialization
bic_values <- numeric(length(lag_orders))  # Correct initialization

# Loop over lag orders and calculate AIC and BIC
for (p in lag_orders) {
  # Fit the VAR model
  var_model <- VAR(VAR_data, p = p, type = "const")
  
  # Store AIC and BIC values
  aic_values[p] <- AIC(var_model)
  bic_values[p] <- BIC(var_model)
}

# Create a data frame for model selection
model_selection <- data.frame(
  Lag_Order = lag_orders,
  AIC = aic_values,
  BIC = bic_values
)

# Print the result
print(model_selection)

# AIC prefers p = 10 and BIC prefers p = 2
# p = 10 has higher adjusted R-squared

# TASK 2.3: VALIDATING VAR MODEL
serial_test <- serial.test(VAR_model, lags.pt = 16, type = "PT.asymptotic") # test residual autocorrelation
# df = number of lags being tested
# p-value > 0.05: Fail to reject null
print(serial_test)

stability_test <- stability(VAR_model)
plot(stability_test) # checking for eigenvalues to lie in the unit circle
# Which assumptions need to be satisfied?

# 3. GRANGER CAUSALITY
# TASK 3.1: EXPLAIN GRANGER CAUSALITY
# - Tests whether one time series can provide useful information for forecasting another. 
# - If including past values of x in a model predicting y significantly improves the forecast 
# compared to a model that only uses past values of y, then x is said to Granger-cause y.
# - So, predictive relationship, but no true causation

# TASK 3.2: PERFORMING CAUSALITY TESTS & INTERPRETING RESULTS
# H_0: The variable cause variable does not Granger-cause the other variable(s) in the model.
# H_0: The coefficients on all lagged terms of cause variable in the equations of the other
# variable(s) are jointly zero.
causality(VAR_model, cause = "INDPRO_stationary")
causality(VAR_model, cause = "CPIAUCSL_stationary")
causality(VAR_model, cause = "FEDFUNDS_stationary")

#4 Impulsive response functions
# Impulse Response Functions (IRFs) are a crucial tool in Vector Autoregressions (VARs) for 
#analyzing the dynamic effects of a one-time shock (impulse) to one variable on the current 
#and future values (responses) of other variables in the system.

#What they measure: An IRF traces the effect of a one-standard-deviation shock to
#one variable on the entire system over time, holding everything else constant.
# Why they're important: IRFs provide insights into the propagation mechanisms of economic shocks,
#helping to interpret the relationships between variables.


# in vars package use of irf()
# irf_result <- irf(var_model, impulse = "variable1", response = "variable2", n.ahead = 10, boot = TRUE)

# Bernanke et al. (2005): The Price Puzzle
#The price puzzle refers to the counterintuitive result where an increase in monetary policy interest rates leads to an increase in the price level.

#Relevance in Bernanke et al. (2005): This puzzle arises in certain VAR models when the model specification or data transformations 
# fail to account for all relevant information, particularly about inflation expectations.

