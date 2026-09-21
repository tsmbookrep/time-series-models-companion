# Estimate the score-driven dynamic Bernoulli model for Boat Race outcomes.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(project_root, "Data")
results_dir <- file.path(project_root, "Results")

if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the readxl package before running this script.")
if (!requireNamespace("nloptr", quietly = TRUE)) stop("Install the nloptr package before running this script.")
source(file.path(script_dir, "Bernoulli_Logl.R"))

# This part loads the Boat Race outcomes time series and fills inherited missing values.
Y <- suppressWarnings(as.numeric(readxl::read_excel(
  file.path(data_dir, "BoatRace46.xlsx"), col_names = FALSE
)[[1]]))
for (i in seq_along(Y)) if (is.na(Y[i]) && i > 1) Y[i] <- Y[i - 1]
Y <- Y[!is.na(Y)] # Drop a nonnumeric spreadsheet heading, if present.
if (!all(Y %in% c(0, 1))) stop("Boat Race outcomes must be coded as 0 or 1.")

# Set to 1 for the simple Markov-chain specification; 0 selects the score-driven AR(1) model.
tsmc <- 0

if (tsmc == 1) {
  start <- c(omega = 0.1, kappa = 0.1)
  C <- rbind(c(1, 0), c(1, 1))
  lower <- c(-Inf, -1); upper <- c(Inf, 1)
} else {
  start <- c(omega = 0, phi = 0.99, kappa = 0.1)
  C <- rbind(c(1, 0, 0), c(1, 0, 1), c(1, 1, 0), c(1, 1, -1))
  lower <- c(-Inf, -1, -Inf); upper <- c(Inf, 1, Inf)
}

# The MATLAB constraints impose 0 <= C %*% par <= 1.
objective <- function(par) Bernoulli_Logl(par, Y, tsmc)$Logl

# SLSQP uses numerical objective derivatives and the exact Jacobian of the linear constraints.
objective_with_gradient <- function(par) {
  step <- 1e-6 * pmax(1, abs(par))
  gradient <- vapply(seq_along(par), function(j) {
    plus <- par; minus <- par
    plus[j] <- plus[j] + step[j]; minus[j] <- minus[j] - step[j]
    (objective(plus) - objective(minus)) / (2 * step[j])
  }, numeric(1))
  list(objective = objective(par), gradient = gradient)
}
inequality_with_jacobian <- function(par) {
  list(
    constraints = c(drop(C %*% par) - 1, -drop(C %*% par)),
    jacobian = rbind(C, -C)
  )
}

fit <- nloptr::nloptr(
  x0 = start, eval_f = objective_with_gradient, eval_g_ineq = inequality_with_jacobian,
  lb = lower, ub = upper,
  opts = list(algorithm = "NLOPT_LD_SLSQP", maxeval = 67000, xtol_rel = 1e-8, print_level = 1)
)
x_s <- fit$solution
fitted <- Bernoulli_Logl(x_s, Y, tsmc)

dir.create(results_dir, showWarnings = FALSE)
write.csv(data.frame(estimate = x_s), file.path(results_dir, "estPar_R.csv"), row.names = FALSE)
print(data.frame(estimate = x_s))

# Plot Boat Race winners and the fitted dynamic probability.
years <- seq(1946, by = 1, length.out = length(Y))
plot(years, Y, pch = 16, cex = 1.1, xlab = "Year", ylab = "Result and probability",
     main = "Dynamic Bernoulli probability", ylim = c(0, 1))
lines(years, fitted$mu, col = "#d95319", lwd = 1.5)
legend("topright", c("Winner (0 Oxf, 1 Cam)", "Dyn. Prob. (AR1)"),
       col = c("black", "#d95319"), pch = c(16, NA), lty = c(NA, 1), bty = "n")
