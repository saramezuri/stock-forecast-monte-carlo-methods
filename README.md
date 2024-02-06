# "Forecasting Hermès Stock Price Using Monte Carlo Methods and Time Series Analysis"
## Sara Mezuri

## Introduction

Hermès is a global luxury fashion brand known for its high-quality products and timeless designs. The company has a strong presence in the market and is valued by investors for its consistent financial performance and strong brand reputation. However, investing in luxury fashion companies can be risky, and investors need to evaluate the potential risks and returns of their investment decisions. To address this challenge, this project aims to use Monte Carlo simulations and time series analysis to forecast the stock prices of Hermès and evaluate different investment strategies.

In this project, we analyze the historical stock prices of Hermès and use time series analysis techniques to identify patterns and trends in the data. We then use the time series model to forecast the future stock prices of Hermès for a given period. Additionally, we conduct Monte Carlo simulations to generate thousands of possible outcomes for the forecasted prices, taking into account various factors that affect the stock prices such as market volatility, industry trends, and economic indicators.

## Historical Stock Market Data 

We will collect monthly stock prices from Investing.com. (https://www.investing.com/equities/Hermès-international-historical-data)

![The data]()

# Time Series Analysis

Based on the graph provided, we can conclude that the time series is non-stationary, characterized by a rising trend but no noticeable seasonality.

Firstly, we start by transforming the time series. We want to transform our data into a stationary time series (a time series with mean and variance independent of the time).

![Hermès Stock Price Time Series]()

From the Time Series plot (fig. 1), we decide that we are dealing with an ARIMA(p, d, q) process. Our goal is to choose the best numerical values of p, d, and q. 

The next step is taking the logarithm of the time series. This is a common data transformation technique used in time series analysis.

In our data set, we notice an increase in the variance overtime, which makes it difficult to forecast the series accurately. Taking the logarithm of the series can help to stabilize the variance, by compressing the range of values for large observations and expanding the range for small observations.

![Hermès Stock Price Logged Time Series]()

Transformations such as logarithms can help to stabilize the variance of a time series. 
But when it comes to the mean of time series, differencing can help stabilize the mean of a time series by removing changes in the level of a time series and therefore eliminating (or reducing) trend and seasonality.

The differenced series is the change between consecutive observations in the original series and can be written as 

$$ X_t ' = X_{t+1} - X_t $$

where $X_t$ is a time series. 

At times, the differenced data might not appear to be stationary, and it might be required to perform a second differencing to achieve a stationary series: 

$$ X_t '' = X_{t+1} ' - X_t ' $$

The next step would be checking whether the first-order differenced time series is stationary or not. To do so, we can use the Augmented Dickey-Fuller (ADF) test. 

The null hypothesis of the Augmented Dickey-Fuller test is that the time series is non-stationary. A small test statistic value suggests that the null hypothesis of a non-stationary time series cannot be rejected. 
 
```{r, warning=FALSE}

library(tseries)

# Perform the ADF test

test_1 <- adf.test(diff1) 

# Print results

print(test_1)

cat("ADF Statistic: ", test_1$statistic, "\n")
cat("p-value: ", test_1$p.value, "\n")


```

```{r, warning = FALSE}

# We take the second-order difference of the logged time series

diff2 <- diff(diff1) 

test_2 <- adf.test(diff2)

# Print results

print(test_2)

cat("ADF Statistic: ", test_2$statistic, "\n")
cat("p-value: ", test_2$p.value, "\n")

```

By comparing the two ADF tests, we notice that both tests suggest that the time series is stationary since the $p$-value is $0.01 < 0.05$. But comparing the test statistics, we notice that the second test statistic is smaller than the first test statistic. 

When the test statistic is very small, it means that the time series is closer to being stationary. 

The first-order differencing removes most of the trend and seasonality of a time series. Therefore, we can conclude that the first-order difference time series is stationary, but still may have some degree of trend, that is not fully removed.  

This is why we take the second-order differencing, which gives us a better stationary time series. 

![Differenced Time Series]()

This conveys that the above is a second-order differenced model, hence $d = 2$, (see Fig.4). 

The autocorrelation function (ACF) and partial autocorrelation function (PACF) plots are commonly used tools in time series analysis. The ACF plot shows the correlation between a time series and its lagged values at different lags, while the PACF plot shows the correlation between the time series and its lagged values after removing the effects of the intervening lags. 

If the ACF plot decays exponentially, it suggests an auto-regressive $(AR(p))$ process.

$$ X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \ldots + \phi_p X_{t-p} + Z_t $$

While, if the PACF plot decays exponentially, this suggests a moving average $(MA(q))$ process.

$$ X_t = Z_t + \theta_1 Z_{t-1} + \theta_2 Z_{t-2} + \ldots + \theta_q Z_{t-q} $$

where $Z_t$ represents a white noise $(0, \sigma^2)$. 

![Autocorrelation and Partial Autocorrelation Functions]()

Based on the graphs of the ACF and PACF, (Fig.5), we conclude that the PACF is exponentially decaying. As a result, we are dealing with an MA process and $p = 0$. Since there is a significant spice at $lag 1$ in the ACF, but none beyond $lag 1$, this is an $MA(1)$ process

$$ X_t = Z_t + \theta_1 Z_{t-1} $$. 

Once the process has been identified, we proceed to fit a model, using the "arima" function on R and the Maximum Likelihood Method. 

The ARIMA (0, 2, 1)  model seems to be the best fit, considering that $AIC = -785.21$ is the smallest value. Additionally, $$Z_t\sim (0, 0.006274)$$, i.e normally distributed.

![Forecast of ARIMA(0, 2, 1) for 24 months]()

# Monte Carlo Simulation

The next step would be simulating the ARIMA model using Monte Carlo simulation. 

```{r}

# We simulate the distribution of future stock prices

mean_forecast <- mean(y.forecast) # the mean of the forecast model
mean_forecast
sd_forecast <- sd(y.forecast) # standard deviation of the forecast model
sd_forecast

sim_prices <- rnorm(1000, mean_forecast, sd_forecast) # this simulates 1000 stock prices

```
The next thing we want to do is calculate the expected returns for each simulated stock price.

For our first scenario, we choose the horizon, which is the length of time an investor plans to hold an investment, to be 2 years. 
Also, the risk-free rate, which is the rate of return an investor can expect to earn on a risk-free investment, we chose it to be $1%$. Thus, the investor is guaranteed to earn $1%$ of the return without taking any risks. 
The risk-free rate varies on several factors as economic conditions, duration of investment, and inflation. Therefore, by choosing different values for the horizon and the risk-free rate, we can calculate different scenarios.

```{r}

horizon <- 2 # investment horizon in years
rf_rate <- 0.01 # risk-free rate 
# The expected return is the average return of the simulated 
# stock price over the investment horizon of 2 years

# This is the return of the simulated stock prices
return <- (sim_prices - mean_forecast)/mean_forecast 

# the compounded annual return
annual_return <- exp(mean(log(1 + return)))^horizon - 1 - rf_rate 

```

After we defined the expected return of each simulated stock price, we repeat the same steps for $n = 1000$. 

```{r}

# We generate 1000 iterations
sim_iterations <- 1000
sim_results <- data.frame(matrix(NA, nrow = sim_iterations, ncol = 2))
colnames(sim_results) <- c("rf_rate")
for (i in 1:sim_iterations) {
  sim_price <- rnorm(1, mean = mean_forecast, sd = sd_forecast)
  sim_return <- (sim_price - mean_forecast)/mean_forecast
  sim_annual_return <- exp(mean(log(1 + sim_return)))^horizon - 1 - rf_rate
  sim_results[i, 1] <- sim_annual_return
}

```
![Expected Annual Return over 2 years with risk-free rate 1%]()

We can choose a different scenario. For instance, we can take a new risk-free rate to be $5$ % instead.
Then, we repeat the same process. (A comparison between different values of the risk-free rate is given below )

In conclusion, the distribution of expected annual returns informs us about the range of potential returns from an investment. We can learn about the risks and potential rewards of an investment by evaluating the shape and spread of the distribution.

The plot (see Fig.7) gives us an idea of the range of the potential return expected from the investment. A skewed distribution toward larger returns may imply a higher-risk investment, whereas a more symmetric distribution may indicate a lower-risk investment. 

The curve (Fig.7) seems to be somewhat wide and flat, suggesting a high degree of uncertainty in the expected returns. 

Furthermore, there are various other methods and measures that assist us to comprehend and assess better the risk and return of our investment. 

One of them would be the Sharpe Ratio, which is a measure of risk-adjusted return and it measures the excess return of an investment over the risk-free rate relative to its volatility. 

The Sharpe Ratio is given by the formula 

$$ Sharpe = \frac{R_p-R_f}{\sigma_p} $$

, where $R_p$ is the return of portofolio, $R_f$ is the risk-free rate and $\sigma_p$ is the standard deviation of the portofolio's return.

In the case when the risk-free rate is taken to be $1$%, the Sharpe Ratio is $-0.3178729$. A negative Sharpe Ratio indicates that on average, the investment's return might be lower than the risk-free rate. A Sharpe ratio of $-0.3178729$ indicates that the simulated stock prices are expected to generate a return that is $0.3178729$ standard deviations below the risk-free rate over a $2$-year investment horizon. 

This shows that the investment may not be a suitable decision when weighing the expected return versus risk, as the projected return is insufficient to compensate for the higher risk involved.

Another method for calculating the risk and loss of an investment would be Value at risk (VaR). The value at risk (VaR) statistic assesses the extent of potential financial losses inside a firm, portfolio, or position over a given time period. 

VaR modeling determines the entity's potential for loss as well as the likelihood that the defined loss will occur.
VaR is calculated by evaluating the amount of potential loss, the probability of the amount of loss occurring, and the time frame. 

Returning to the scenario we were studying previously, we want to calculate the $95$% VaR for a $2$-years investment horizon. 


```{r, fig.width=8, fig.height=6, fig.cap = "Expected Annual Return over 2 years with risk-free rate 0.5%"}

horizon <- 2 # investment horizon in years
rf_rate <- 0.01 # risk-free rate 
initial_investment <- 100000 # initial investment

# Calculate the investment return for each simulated stock price
returns <- (sim_prices - mean_forecast)/mean_forecast

# Sort the returns
sorted_returns <- sort(returns)

# Determine the percentile of the sorted returns
percentile <- qnorm(0.05, mean = mean(returns), sd = sd(returns))

selected_return <- sorted_returns[which(sorted_returns <= percentile)
                                  [length(which(sorted_returns <= percentile))]]

# Define VaR
var <- abs(initial_investment * selected_return)

cat("The VaR for a two-years investment horizon and a 95% confidence level is:", var, "\n")

```

The VaR (value at risk) of $18655.95$ implies a $5$% chance of the investment losing at least that amount over the specified investment horizon (in this case, two years) with $95$% confidence. For example, if you invest $100,000 for two years, there is a $5$% chance that the investment will lose at least $18655.95.

In conclusion, several elements are considered when measuring the risk of an investment, such as corporate performance, industry trends (in our case, fashion trends), market circumstances, regulatory environment, volatility, and so on. 

The past performance of the stock price, market trends, and volatility, as indicated in the ARIMA model, are taken into account to estimate the risk of investing in forecasting Hermès Stock Prices.  Furthermore, as a measure of the overall risk of investing in the company, the risk-free rate of return is taken into account when calculating the expected return and compounded annual return. The chart provided below (Fig.8) illustrates the varying impact that different risk-free rates have on each distribution.

To summarize, an investment with a wide normal curve, a negative Sharpe ratio and a positive VaR suggests an investment that is risky and most likely does not generate a positive return that compensates for the taken risk. 

![A Comparison of Expected Annual Returns over 2 years with Different Risk-Free Rates]

To summarize, an investment with a wide normal curve, a negative Sharpe ratio and a positive VaR suggests an investment that is risky and most likely does not generate a positive return that compensates for the taken risk. 

\newpage

# References

1. Rizzo, Maria L. (2019). Statistical Computing with R, Second Edition

2. Brockwell, Peter J., Davis, Richard A. (2016). Introduction to Time Series and Forecasting, Third Edition

3. Prabhakaran, S. (2021). ARIMA Model - Complete Guide to Time Series Forecasting in Python (https://www.machinelearningplus.com/time-series/arima-model-time-series-forecasting-python/)

4. CFI Team. (2023). Annual Return
(https://corporatefinanceinstitute.com/resources/capital-markets/annual-return/)

5. Stationarity Testing
(https://rpubs.com/richkt/269797)

6. Fernando, J. (2022). Sharpe Ratio Formula and Definition With Examples
(https://www.investopedia.com/terms/s/sharperatio.asp)

7. Kenton, W. (2023). Understanding Value at Risk (VaR) and How It Is Compted. 
(https://www.investopedia.com/terms/v/var.asp)
