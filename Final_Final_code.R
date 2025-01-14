library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

options(scipen=999)
set.seed(1201)
################################### Parameters #################################
S0 <- 100 # Underlying asset price
r <- 0.05 # risk-free rate
nsteps <- 12 # Time steps
npaths <- 500000 # Number of paths
T <- 1 # Maturity (years)
dt <- T/nsteps # Discritization (time)
K <- 110 # Strike price
sigma <- 0.1 # Volatility
discount <- exp(-r * dt) # For discounted payoff

#################################### Bermudan ##################################
# Longstaff function
simulate_bermudan_option <- function(npaths, nsteps, S0, r, sigma, K, discount) {
  # Time step size
  dt <- 1 / nsteps
  
  # Simulating stock price
  BM <- matrix(rnorm(npaths * nsteps, mean = 0, sd = sqrt(dt)), ncol = nsteps)
  S <- matrix(0, nrow = npaths, ncol = nsteps + 1)
  S[, 1] <- S0
  
  for (t in 1:nsteps) {
    S[, t + 1] <- S[, t] * exp((r - 0.5 * sigma^2) * dt + sigma * BM[, t])
  }
  
  # Payoff matrix for a Bermudan Call option
  payoff <- pmax(S - K, 0)
  
  # Initialize discounted cashflows and cashflow matrix
  disc_cashflow <- payoff[, nsteps + 1] * discount
  cashflow <- matrix(0, nrow = npaths, ncol = nsteps + 1)
  
  # Backward induction for Bermudan option pricing
  for (t in (nsteps + 1):1) {
    if (t == (nsteps + 1)) {
      cashflow[, t] <- payoff[, t]
      next
    }
    
    # Identify in-the-money paths
    in_the_money <- which(payoff[, t] > 0)
    continuation_value <- numeric(npaths)
    
    if (length(in_the_money) > 0) {
      # Fit regression to estimate continuation value
      X <- S[in_the_money, t]
      Y <- payoff[in_the_money, t + 1] * discount^(t)
      
      # Compute Laguerre basis functions
      L0 <- exp(-X / 2)
      L1 <- exp(-X / 2) * (1 - X)
      L2 <- exp(-X / 2) * (1 - 2 * X + X^2 / 2)
      L3 <- exp(-X / 2) * (exp(X) / (3 * 2)) * exp(-X) * (6 - 18 * X + 9 * X^2 - X^3)
      
      continuation_value[in_the_money] <- predict(lm(Y ~ L0 + L1 + L2 + L3),
                                                  newdata = data.frame(L0 = L0, L1 = L1, L2 = L2, L3 = L3))
      
      # Exercise or continue decision
      exercise <- payoff[, t] > continuation_value
      cashflow[exercise, t] <- payoff[exercise, t]
      cashflow[exercise, (t + 1):(nsteps + 1)] <- 0
      disc_cashflow <- ifelse(exercise, payoff[, t], disc_cashflow)
    }
    
    if (t != 1) {
      disc_cashflow <- disc_cashflow * discount
    }
  }
  
  cashflow[is.na(cashflow)] <- 0
  bermudan_price <- mean(disc_cashflow)
  
  return(list(bermudan_price = bermudan_price, cashflows = cashflow, disc_cashflow = disc_cashflow, S = S))
}

bermudan_price <- simulate_bermudan_option(npaths, nsteps, S0, r, sigma, K, discount)

################################### European ###################################
# Black-Scholes Function
black_scholes <- function(S0, K, T, r, sigma) {
  # Calculate d1 and d2
  d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  # N(d1) and N(d2)
  Nd1 <- pnorm(d1)  # Cumulative normal distribution
  Nd2 <- pnorm(d2)
  
  # call price 
  price <- S0 * Nd1 - K * exp(-r * T) * Nd2

  return(price)
}

european_price <- black_scholes(S0, K, T, r, sigma)



#################################### Binomial ##################################
binomial_model <- function(S0, K, r, T, sigma, nsteps, option_type = "call") {
  # Parameters
  dt <- T / nsteps
  u <- exp(sigma * sqrt(dt))    # Up factor
  d <- 1 / u                   # Down factor
  p <- (exp(r * dt) - d) / (u - d) # Risk-neutral probability
  discount <- exp(-r * dt) # For discounted payoff
  
  # Initialize stock price tree
  stock_tree <- matrix(0, nrow = nsteps + 1, ncol = nsteps + 1)
  for (i in 0:nsteps) {
    for (j in 0:i) {
      stock_tree[j + 1, i + 1] <- S0 * (u^j) * (d^(i - j))
    }
  }
  
  # Initialize option value tree
  option_tree <- matrix(0, nrow = nsteps + 1, ncol = nsteps + 1)
  
  # Terminal payoffs
  if (option_type == "call") {
    option_tree[, nsteps + 1] <- pmax(stock_tree[, nsteps + 1] - K, 0)
  } else {
    option_tree[, nsteps + 1] <- pmax(K - stock_tree[, nsteps + 1], 0)
  }
  
  # Backward induction for Bermudan option pricing
  for (i in (nsteps):1) {
    for (j in 0:(i - 1)) {
      continuation_value <- discount * (p * option_tree[j + 2, i + 1] + 
                                          (1 - p) * option_tree[j + 1, i + 1])
      intrinsic_value <- if (option_type == "call") {
        max(0, stock_tree[j + 1, i] - K)
      } else {
        max(0, K - stock_tree[j + 1, i])
      }
      option_tree[j + 1, i] <- max(intrinsic_value, continuation_value)
    }
  }
  
  # Initialize optimal stopping matrix
  optimal_stopping <- matrix(0, nrow = nsteps + 1, ncol = nsteps + 1)
  
  # Extract the optimal stopping points
  for (i in 1:nsteps) {  # Adjusted to match matrix dimensions
    for (j in 1:i) {
      continuation_value <- discount * (p * option_tree[j + 1, i + 1] + 
                                          (1 - p) * option_tree[j, i + 1])
      if (option_tree[j, i] > continuation_value) {
        optimal_stopping[j, i] <- 1
      }
    }
  }
  
  list(option_price = option_tree[1, 1], 
       stock_tree = stock_tree, 
       option_tree = option_tree,
       optimal_stopping = optimal_stopping)
}

binomial_price <- binomial_model(S0, K, r, T, sigma, nsteps, option_type = "call")

################################## Comparison ##################################
cat("Bermudan Call Option Price:", bermudan_price$bermudan_price)
cat("European Call Option Price:", european_price, "\n")
cat("Binomial Call Option Price:", binomial_price$option_price)

# European Early Exercise Value (simulated as Bermudan - European)
european_price - bermudan_price$bermudan_price
((european_price - bermudan_price$bermudan_price)/bermudan_price$bermudan_price) * 100

# Binomial Early Exercise Value (simulated as Bermudan - binomial)
binomial_price$option_price - bermudan_price$bermudan_price
((binomial_price$option_price - bermudan_price$bermudan_price)/bermudan_price$bermudan_price) * 100

# standard error for Bermudan
sd(bermudan_price$disc_cashflow) / sqrt(npaths)
