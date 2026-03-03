#SCENARIO 1 Targeting entire Scalloped Hammerhead Sharks

library(deSolve)  # ODE solver
library(ggplot2)  # plotting
library(dplyr)    # pipes and data manipulation
library(tidyr)    # pivot_longer
library(grid)     # for unit() in arrow plots

make_arrows <- function(df, k, x, y) {
  df %>%
    mutate(row = row_number()) %>%
    filter((row - 1) %% (k + 2) < 2) %>%
    mutate(
      xend = lead({{ x }}),
      yend = lead({{ y }})
    ) %>%
    filter(row %% 2 == 1)
}

bioecon_logistic <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    # Numerical guard: keep states non-negative (avoids small negative drift)
    S <- max(S, 0)
    E <- max(E, 0)
    
    H  <- q * E * S
    pi <- (p * q * S - c) * E
    
    dS <- r * S * (1 - S / K) - H
    dE <- theta * pi
    
    list(c(dS, dE))
  })
}

# Parameter set (units are abstract; keep the same symbols as prior tutorials)
parms_L <- c(r = 0.5, K = 60, q = 0.007, p = 1, c = 0.5, theta = 0.02)

# Grids for plotting
S_grid <- seq(0, parms_L["K"], length.out = 500)

# Nullclines
E_null_S <- (parms_L["r"]/parms_L["q"]) * (1 - S_grid/parms_L["K"])  # dS=0, S>0 branch
S_null_E <- parms_L["c"]/(parms_L["p"]*parms_L["q"])                # dE=0, E>0 branch

df_null <- data.frame(S = S_grid, E = E_null_S)

ggplot(df_null, aes(x = S, y = E)) +
  # dS/dt = 0 interior branch
  geom_line(linewidth = 1, linetype = "dashed") +
  # dS/dt = 0 extinction branch
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
  # dE/dt = 0 interior branch: S = c/(pq)
  geom_vline(xintercept = S_null_E, linetype = "dashed", linewidth = 1) +
  # dE/dt = 0 boundary branch: E = 0
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 1) +
  labs(title = "Nullclines in the (S, E) plane (Logistic) Scenario 1",
       x = "Stock (S)", y = "Effort (E)") +
  coord_cartesian(xlim = c(0, parms_L["K"])) +
  theme_minimal()

times <- seq(0, 100, by = 0.1)

# A small set of initial conditions (extend in class)
inits <- list(
  c(S = 30,  E = 32),
  c(S = 30,  E = 3.2)
)

outs <- lapply(seq_along(inits), function(i) {
  out <- ode(y = inits[[i]],
             times = times,
             func = bioecon_logistic,
             parms = parms_L)
  out <- as.data.frame(out)
  out$ic <- paste0("S0=", inits[[i]]["S"], ", E0=", inits[[i]]["E"])
  out
})

traj_L <- bind_rows(outs)

traj_L_long <- traj_L %>%
  pivot_longer(cols = c("S","E"), names_to = "var", values_to = "val")

ggplot(traj_L_long, aes(x = time, y = val, color = ic)) +
  geom_line() +
  facet_wrap(~var, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "ODE-system Scenario 1", y = "Size", x = "Time", color = "Initial condition")

bioecon_allee <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    # Numerical guard: keep states non-negative (avoids small negative drift)
    S <- max(S, 0)
    E <- max(E, 0)
    
    H  <- q * E * S
    pi <- (p * q * S - c) * E
    
    dS <- r * S * (1 - S / K) * (S / A - 1) - H
    dE <- theta * pi
    
    list(c(dS, dE))
  })
}

# Parameters (Allee threshold inside (0, K))
parms_A <- c(r = 0.5, K = 60, A = 20, q = 0.007, p = 1, c = 0.5, theta = 0.1)

# Grid for stock (avoid S=0 since we focus on interior nullclines)
S_grid <- seq(1e-6, parms_A["K"], length.out = 800)

# Interior nullcline: dS/dt = 0 for S>0
# E = (r/q) (1 - S/K) (S/A - 1)
E_null_S <- (parms_A["r"]/parms_A["q"]) *
  (1 - S_grid/parms_A["K"]) *
  (S_grid/parms_A["A"] - 1)

df_null <- data.frame(S = S_grid, E = E_null_S)

# Interior nullcline: dE/dt = 0 for E>0  ->  S = c/(pq)
S_null_E <- as.numeric(parms_A["c"]/(parms_A["p"]*parms_A["q"]))

# Plot only the interior nullclines
ggplot(df_null, aes(x = S, y = E)) +
  geom_line(linewidth = 1, linetype = "dashed") +                 # dS/dt=0 (interior)
  geom_vline(xintercept = S_null_E, linetype = "dashed") +        # dE/dt=0 (interior)
  coord_cartesian(ylim = c(0, max(0, max(E_null_S, na.rm = TRUE)))) +
  labs(
    title = "Interior nullclines in the (S, E) plane (Allee effect) Scenario 1",
    subtitle = "Only interior branches are shown: S>0 and E>0",
    x = "Stock (S)", y = "Effort (E)"
  ) +
  theme_minimal()

times <- seq(0, 2000, by = 0.1)

inits <- list(
  c(S = 70,  E = 32),
  c(S = 70,  E = 3.2)
)

outs_A <- lapply(seq_along(inits), function(i) {
  out <- ode(y = inits[[i]],
             times = times,
             func = bioecon_allee,
             parms = parms_A)
  out <- as.data.frame(out)
  out$ic <- paste0("S0=", inits[[i]]["S"], ", E0=", inits[[i]]["E"])
  out
})

traj_A <- bind_rows(outs_A)

traj_A_long <- traj_A %>%
  pivot_longer(cols = c("S","E"), names_to = "var", values_to = "val")

ggplot(traj_A_long, aes(x = time, y = val, color = ic)) +
  geom_line() +
  facet_wrap(~var, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Open access trajectories (Allee) Scenario 1", y = "Size", x = "Time",color= "Initial condition")

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

# Social–ecological model: (C, S) system
# State variables:
#   C = number of cooperators
#   S = stock biomass
#
# Parameters:
#   N = community size (constant), D = N - C
#   ec, ed = effort levels of cooperators/defectors
#   r, K = logistic growth parameters
#   a = social pressure (conversion of defectors by cooperators)
#   b = temptation to defect (erodes cooperation)

ses_fishery <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    # Guards: keep within meaningful bounds (prevents numerical drift)
    C <- min(max(C, 0), N)
    S <- max(S, 0)
    
    D <- N - C
    
    # Ecological dynamics (logistic growth minus harvest by both groups)
    dS <- r * S * (1 - S / K) - ed * D * S - ec * C * S
    
    # Social dynamics (norm enforcement - temptation to defect)
    dC <- a * C * (N - C) - C * b * (ed / ec)
    
    list(c(dC, dS))
  })
}

# Baseline parameters (from the exercise statement)
parms_ses <- c(
  a = 0.001,
  b = 0.01,
  r = 0.5,
  ec = 0.002,
  ed = 0.004,
  N = 100,
  K = 60
)

# Initial conditions
state0_ses <- c(C = 50, S = 50)

# Time grid
times <- seq(0, 500, by = 0.1)

# Simulate
out_ses <- ode(y = state0_ses, times = times, func = ses_fishery, parms = parms_ses)
df_ses  <- as.data.frame(out_ses)

# Time series plot
df_long <- df_ses %>%
  pivot_longer(cols = c("C", "S"), names_to = "var", values_to = "val")

ggplot(df_long, aes(x = time, y = val, color = var)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Scenario 1 Social–ecological system: cooperators and fish stock",
       x = "Time", y = NULL)

