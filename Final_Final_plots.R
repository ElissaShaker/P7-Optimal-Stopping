start_time <- Sys.time()
################################ Determine npaths ##############################
# Simulation for each npaths
# Sequence of npaths
npaths_seq <- seq(1000, 1000000, by = 10000)
results <- numeric(length(npaths_seq))

# Simulation for each npaths
for (i in seq_along(npaths_seq)) {
  npaths <- npaths_seq[i]
  BM <- matrix(rnorm(npaths * nsteps, mean = 0, sd = sqrt(dt)), ncol = nsteps)
  S <- matrix(0, nrow = npaths, ncol = nsteps + 1)
  S[, 1] <- S0
  for (t in 1:nsteps) {
    S[, t + 1] <- S[, t] * exp((r - 0.5 * sigma^2) * dt + sigma * BM[, t])
  }
  # Example metric: Average terminal price
  results[i] <- mean(S[, nsteps + 1])
}
# Create a data frame for plotting
data <- data.frame(
  npaths = npaths_seq,
  mean_terminal_price = results
)

# ggplot2 plot
ggplot(data, aes(x = npaths, y = mean_terminal_price)) +
  geom_line(color = "blue", size = 1) +
  #geom_point(color = "red") +
  theme_minimal() +
  labs(
    title = "Average Stock Price at Maturity Time with Number of Paths",
    x = "Number of Paths",
    y = "Average Stock Price at Maturity Time"
  )

################################### Paths plot #################################
# Plotting a subset of simulated price paths, select 10 random paths.
df_paths <- as.data.frame(bermudan_price$S[sample(1:npaths, 20), ]) #(S[1:npaths, 10])
df_paths <- df_paths %>% 
  mutate(Path = row_number()) %>% 
  pivot_longer(-Path, names_to = "Time", values_to = "Price") %>% 
  mutate(Time = as.numeric(gsub("V", "", Time)) - 1 )

ggplot(df_paths, aes(x = Time, y = Price, color = factor(Path))) +
  geom_line() +
  labs(title = "Simulated Stock Price Paths",
       x = "Time Step",
       y = "Stock Price") +
  theme_minimal()+
  theme(legend.position = "none")

################################### Histogram ##################################
# Prepare data for histogram
# Replace all non-zero values in cashflow with 1 to count exercises
exercise_matrix <- ifelse(bermudan_price$cashflows != 0, 1, 0)

# Sum up exercises for each column (time step)
exercise_counts <- colSums(exercise_matrix)

# Create a data frame for visualization
exercise_df <- data.frame(
  Time = 0:nsteps,  # Time periods
  Exercises = exercise_counts
)

# Plot histogram
ggplot(exercise_df, aes(x = Time, y = Exercises)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Exercise Distribution Over Time",
       x = "Time Step",
       y = "Number of Exercises") +
  theme_minimal()


end_time <- Sys.time()
end_time - start_time


