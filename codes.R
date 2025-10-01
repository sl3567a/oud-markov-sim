library(tidyverse)

# Define original states and initialize the transition matrix
states <- c("Compulsive Use", "On MOUD", "Withdrawal", "Fatal Overdose", "Death")
P_statusquo <- matrix(0, nrow = 5, ncol = 5)
colnames(P_statusquo) <- rownames(P_statusquo) <- states

# Fill in transition probabilities for the status quo scenario
P_statusquo["Compulsive Use", ] <- c(0.9964504, 0.002675905, 0.00058, 0.0001253437, 0.0001683808)
P_statusquo["On MOUD", ] <- c(0.002448776*(1-0.0000713976), 0.9406888, 0.0567, 0.002448776*0.0000713976, 0.000162375)
P_statusquo["Withdrawal", ] <- c(0.00329*(1-0.0002506874), 0, 0.9965476, 0.00329*0.0002506874, 0.000162375)
P_statusquo["Fatal Overdose", ] <- c(0, 0, 0, 0, 1)
P_statusquo["Death", ] <- c(0, 0, 0, 0, 1)

# Copy a version for the intervention scenario
P_intervention <- P_statusquo

# Increase Compulsive Use → On MOUD transition probability
new_MOUD_prob <- 0.01
row <- P_intervention["Compulsive Use", ]
row_target <- "On MOUD"
other_states <- setdiff(states, row_target)
residual <- 1 - new_MOUD_prob
row[other_states] <- row[other_states] / sum(row[other_states]) * residual
row[row_target] <- new_MOUD_prob
P_intervention["Compulsive Use", ] <- row

# Reduce fatal relapse probability from On MOUD
relapse_prob_onmoud <- 0.002448776
old_od_prob_onmoud <- relapse_prob_onmoud * 0.0000713976
new_od_prob_onmoud <- relapse_prob_onmoud * 0.000007  # 更低的死亡率
row <- P_intervention["On MOUD", ]
delta <- old_od_prob_onmoud - new_od_prob_onmoud
row["Fatal Overdose"] <- new_od_prob_onmoud
row["Compulsive Use"] <- row["Compulsive Use"] + delta
P_intervention["On MOUD", ] <- row

# Reduce fatal relapse probability from Withdrawal
relapse_prob_withdrawal <- 0.00329
old_od_prob_withdrawal <- relapse_prob_withdrawal * 0.0002506874
new_od_prob_withdrawal <- relapse_prob_withdrawal * 0.0002
row <- P_intervention["Withdrawal", ]
delta <- old_od_prob_withdrawal - new_od_prob_withdrawal
row["Fatal Overdose"] <- new_od_prob_withdrawal
row["Compulsive Use"] <- row["Compulsive Use"] + delta
P_intervention["Withdrawal", ] <- row

# Final check: verify that row sums are 1 
row_sums_statusquo <- rowSums(P_statusquo)
row_sums_intervention <- rowSums(P_intervention)

simulate_markov <- function(P, initial, weeks) {
  pop_states <- matrix(0, nrow = weeks + 1, ncol = length(initial))
  colnames(pop_states) <- names(initial)
  pop_states[1, ] <- initial
  
  for (t in 1:weeks) {
    pop_states[t + 1, ] <- pop_states[t, ] %*% P
  }
  
  df <- as.data.frame(pop_states)
  df$Week <- 0:weeks
  df_long <- df %>%
    pivot_longer(-Week, names_to = "State", values_to = "Population")
  
  df_wide <- df_long %>%
    pivot_wider(names_from = State, values_from = Population) %>%
    arrange(Week) %>%
    mutate(
      Cumulative_Overdose_Deaths = cumsum(`Fatal Overdose`),
      Alive = `Compulsive Use` + `On MOUD` + `Withdrawal`
    )
  
  return(df_wide)
}

# Initial distribution and simulation duration
initial <- c("Compulsive Use" = 300000, "On MOUD" = 0, "Withdrawal" = 0,
             "Fatal Overdose" = 0, "Death" = 0)
weeks <- 52 * 5

# Simulate both status quo and intervention
df_statusquo <- simulate_markov(P_statusquo, initial, weeks)
df_intervention <- simulate_markov(P_intervention, initial, weeks)
df_statusquo$Scenario <- "Status Quo"
df_intervention$Scenario <- "Intervention"
df_compare <- bind_rows(df_statusquo, df_intervention)

write_csv(df_statusquo, "df_statusquo.csv")
write_csv(df_intervention, "df_intervention.csv")

ggplot(df_compare, aes(x = Week, y = Cumulative_Overdose_Deaths, color = Scenario)) +
  geom_line(size = 1.2) +
  labs(title = "Cumulative Overdose Deaths: Status Quo vs Intervention",
       y = "Cumulative Deaths") +
  theme_minimal()

ggplot(df_compare, aes(x = Week, y = Alive, color = Scenario)) +
  geom_line(size = 1.2) +
  labs(title = "Total Alive Population Over Time",
       y = "Number of Alive Individuals") +
  theme_minimal()

# Set initial values
initial <- c("Compulsive Use" = 300000, "On MOUD" = 0, "Withdrawal" = 0,
             "Fatal Overdose" = 0, "Death" = 0)
weeks <- 52 * 5

# Define parameter ranges for DSA）
dsa_params <- list(
  list(name = "CU_to_MOUD", base = 0.002675905, values = c(0.001338, 0.004014)),
  list(name = "OD_OnMOUD", base = 0.0000713976, values = c(0.000001, 0.0001)),
  list(name = "OD_Withdrawal", base = 0.0002506874, values = c(0.000001, 0.0001))
)

# Build transition matrix function
build_matrix_dsa <- function(cu_to_moud, od_onmoud, od_withdrawal) {
  states <- c("Compulsive Use", "On MOUD", "Withdrawal", "Fatal Overdose", "Death")
  P <- matrix(0, nrow = 5, ncol = 5)
  colnames(P) <- rownames(P) <- states
  
  # Fixed relapse probabilities
  relapse_on <- 0.002448776
  relapse_wd <- 0.00329
  
  # Compulsive Use row
  base_row <- c(Compulsive = 0.9964504, MOUD = 0.002675905, WD = 0.00058,
                OD = 0.0001253437, D = 0.0001683808)
  
  names(base_row) <- states
  base_row["On MOUD"] <- cu_to_moud
  other_states <- setdiff(states, "On MOUD")
  residual <- 1 - cu_to_moud
  base_row[other_states] <- base_row[other_states] / sum(base_row[other_states]) * residual
  row_cu <- base_row
  P["Compulsive Use", ] <- row_cu
  
  # On MOUD row
  row_on <- c(
    "Compulsive Use" = relapse_on * (1 - od_onmoud),
    "On MOUD" = 0.9406888,
    "Withdrawal" = 0.0567,
    "Fatal Overdose" = relapse_on * od_onmoud,
    "Death" = 0.000162375
  )
  row_on <- row_on / sum(row_on)
  P["On MOUD", ] <- row_on
  
  # Withdrawal row
  row_wd <- c(
    "Compulsive Use" = relapse_wd * (1 - od_withdrawal),
    "On MOUD" = 0,
    "Withdrawal" = 0.9965476,
    "Fatal Overdose" = relapse_wd * od_withdrawal,
    "Death" = 0.000162375
  )
  row_wd <- row_wd / sum(row_wd)
  P["Withdrawal", ] <- row_wd
  
  # Absorbing states
  P["Fatal Overdose", "Death"] <- 1
  P["Death", "Death"] <- 1
  
  return(P)
}


# Store DSA results
dsa_results <- data.frame()

for (p in dsa_params) {
  for (val in p$values) {
    
    # Set parameter for each round
    cu_moud <- if (p$name == "CU_to_MOUD") val else 0.002675905
    od_on <- if (p$name == "OD_OnMOUD") val else 0.0000713976
    od_wd <- if (p$name == "OD_Withdrawal") val else 0.0002506874
    
    P <- build_matrix_dsa(cu_moud, od_on, od_wd)
    df <- simulate_markov(P, initial, weeks)
    
    # Extract distribution at final week
    final <- df[nrow(df), ]
    
    dsa_results <- rbind(
      dsa_results,
      data.frame(
        Parameter = p$name,
        Value = val,
        OverdoseDeaths = final$Cumulative_Overdose_Deaths,
        Alive = final$Alive,
        CompulsiveUse = final$`Compulsive Use`,
        OnMOUD = final$`On MOUD`,
        Withdrawal = final$Withdrawal
      )
    )
  }
}

beta_params <- function(mean, se) {
  var <- se^2
  temp <- mean * (1 - mean) / var - 1
  alpha <- mean * temp
  beta <- (1 - mean) * temp
  return(c(alpha = alpha, beta = beta))
}

beta_params(0.002675905, 0.0005)

# Updated transition matrix constructor
build_transition_matrix <- function(cu_to_moud, onmoud_to_od, withdrawal_to_od) {
  states <- c("Compulsive Use", "On MOUD", "Withdrawal", "Fatal Overdose", "Death")
  P <- matrix(0, nrow = 5, ncol = 5)
  colnames(P) <- rownames(P) <- states
  
  # Compulsive Use 
  base_row <- c(
    "Compulsive Use" = 0.9964504,
    "On MOUD" = 0.002675905,
    "Withdrawal" = 0.00058,
    "Fatal Overdose" = 0.0001253437,
    "Death" = 0.0001683808
  )
  row_target <- "On MOUD"
  other_states <- setdiff(states, row_target)
  residual <- 1 - cu_to_moud
  row <- numeric(length(states))
  names(row) <- states
  row[row_target] <- cu_to_moud
  row[other_states] <- base_row[other_states] / sum(base_row[other_states]) * residual
  P["Compulsive Use", ] <- row
  
  # On MOUD 
  relapse_prob_onmoud <- 0.002448776
  row <- numeric(length(states))
  names(row) <- states
  row["Compulsive Use"] <- relapse_prob_onmoud * (1 - onmoud_to_od)
  row["Fatal Overdose"] <- relapse_prob_onmoud * onmoud_to_od
  row["On MOUD"] <- 0.9406888
  row["Withdrawal"] <- 0.0567
  row["Death"] <- 0.000162375
  row["On MOUD"] <- 1 - sum(row[setdiff(states, "On MOUD")])
  P["On MOUD", ] <- row
  
  # Withdrawal
  relapse_prob_withdrawal <- 0.00329
  row <- numeric(length(states))
  names(row) <- states
  row["Compulsive Use"] <- relapse_prob_withdrawal * (1 - withdrawal_to_od)
  row["Fatal Overdose"] <- relapse_prob_withdrawal * withdrawal_to_od
  row["Withdrawal"] <- 0.9965476
  row["Death"] <- 0.000162375
  row["Withdrawal"] <- 1 - sum(row[setdiff(states, "Withdrawal")])
  P["Withdrawal", ] <- row
  
  # Absorbing states
  P["Fatal Overdose", ] <- c(0, 0, 0, 0, 1)
  P["Death", ] <- c(0, 0, 0, 0, 1)
  
  return(P)
}

# Markov simulation function
simulate_markov <- function(P, initial, weeks) {
  pop_states <- matrix(0, nrow = weeks + 1, ncol = length(initial))
  colnames(pop_states) <- names(initial)
  pop_states[1, ] <- initial
  for (t in 1:weeks) {
    pop_states[t + 1, ] <- pop_states[t, ] %*% P
  }
  df <- as.data.frame(pop_states)
  df$Week <- 0:weeks
  df$Cumulative_Overdose_Deaths <- cumsum(df$`Fatal Overdose`)
  df$Alive <- df$`Compulsive Use` + df$`On MOUD` + df$`Withdrawal`
  return(df)
}

# Initial setup
initial <- c("Compulsive Use" = 300000, "On MOUD" = 0, "Withdrawal" = 0, "Fatal Overdose" = 0, "Death" = 0)
weeks <- 52 * 5
n_sim <- 1000

# PSA parameters
alpha_cu <- 28.56; beta_cu <- 10645.42
alpha_onod <- 0.0306; beta_onod <- 174835.67
alpha_wdod <- 0.6802; beta_wdod <- 824759.19

# Store simulation results
results <- data.frame()

for (i in 1:n_sim) {
  cu_to_moud <- rbeta(1, alpha_cu, beta_cu)
  onmoud_to_od <- rbeta(1, alpha_onod, beta_onod)
  withdrawal_to_od <- rbeta(1, alpha_wdod, beta_wdod)
  
  # Build transition matrices for two scenarios
  P_stat <- build_transition_matrix(cu_to_moud, onmoud_to_od, withdrawal_to_od)
  P_int <- build_transition_matrix(0.01, 0.000007, 0.0002)
  
  # Run simulation
  df_stat <- simulate_markov(P_stat, initial, weeks)
  df_int <- simulate_markov(P_int, initial, weeks)
  
  # Extract and store results
  results <- rbind(
    results,
    # 5-year outcomes
    data.frame(
      Sim = i,
      Scenario = "Status Quo",
      Time = "Year 5",
      OverdoseDeaths = tail(df_stat$Cumulative_Overdose_Deaths, 1),
      Alive = tail(df_stat$Alive, 1),
      CompulsiveUse = tail(df_stat$`Compulsive Use`, 1),
      OnMOUD = tail(df_stat$`On MOUD`, 1),
      Withdrawal = tail(df_stat$Withdrawal, 1)
    ),
    data.frame(
      Sim = i,
      Scenario = "Intervention",
      Time = "Year 5",
      OverdoseDeaths = tail(df_int$Cumulative_Overdose_Deaths, 1),
      Alive = tail(df_int$Alive, 1),
      CompulsiveUse = tail(df_int$`Compulsive Use`, 1),
      OnMOUD = tail(df_int$`On MOUD`, 1),
      Withdrawal = tail(df_int$Withdrawal, 1)
    ),
    
    # 1-year outcomes (week 52)
    data.frame(
      Sim = i,
      Scenario = "Status Quo",
      Time = "Year 1",
      OverdoseDeaths = df_stat$Cumulative_Overdose_Deaths[53],
      Alive = df_stat$Alive[53],
      CompulsiveUse = df_stat$`Compulsive Use`[53],
      OnMOUD = df_stat$`On MOUD`[53],
      Withdrawal = df_stat$Withdrawal[53]
    ),
    data.frame(
      Sim = i,
      Scenario = "Intervention",
      Time = "Year 1",
      OverdoseDeaths = df_int$Cumulative_Overdose_Deaths[53],
      Alive = df_int$Alive[53],
      CompulsiveUse = df_int$`Compulsive Use`[53],
      OnMOUD = df_int$`On MOUD`[53],
      Withdrawal = df_int$Withdrawal[53]
    )
  )
}

# Summarize simulation results
summary_results <- results %>%
  group_by(Scenario, Time) %>%
  summarise(
    Mean_OverdoseDeaths = mean(OverdoseDeaths),
    UI_Low_OverdoseDeaths = quantile(OverdoseDeaths, 0.025),
    UI_High_OverdoseDeaths = quantile(OverdoseDeaths, 0.975),
    
    Mean_Alive = mean(Alive),
    UI_Low_Alive = quantile(Alive, 0.025),
    UI_High_Alive = quantile(Alive, 0.975),
    
    Mean_CompulsiveUse = mean(CompulsiveUse),
    UI_Low_CompulsiveUse = quantile(CompulsiveUse, 0.025),
    UI_High_CompulsiveUse = quantile(CompulsiveUse, 0.975),
    
    Mean_OnMOUD = mean(OnMOUD),
    UI_Low_OnMOUD = quantile(OnMOUD, 0.025),
    UI_High_OnMOUD = quantile(OnMOUD, 0.975),
    
    Mean_Withdrawal = mean(Withdrawal),
    UI_Low_Withdrawal = quantile(Withdrawal, 0.025),
    UI_High_Withdrawal = quantile(Withdrawal, 0.975)
  )

print(summary_results)

