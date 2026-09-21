# Estimate the dynamic-scale Beta-t-EGARCH model on an equity price series.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(project_root, "Data")
results_dir <- file.path(project_root, "Results")

if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the readxl package before running this script.")
source(file.path(script_dir, "Beta_univ_t_Egarch_Logl.R"))
source(file.path(script_dir, "Beta_univ_t_Egarch_Logl_fact.R"))
source(file.path(script_dir, "Info_Crit.R"))

# This code sets up the dynamic-scale Beta-t-EGARCH model. The active input is the NSX price series.
prices <- as.numeric(readxl::read_excel(file.path(data_dir, "NSX Close---.xlsx"), col_names = FALSE)[[1]])
Y <- diff(log(prices))
t <- length(Y)

# If set to 1, scale dynamics are integrated; otherwise they follow a stationary AR(1) recursion.
I_s <- 0
# If set to 1, standardize the score by its information matrix.
Inf_s <- 1
# If set to 1, add the leverage term kappa_l * sign(-y_t) * (u_t + 1).
lev_s <- 1
Constr_s <- 1
kk_s <- 0.01

# Plot the unconditional distribution and autocorrelation diagnostics.
par(mfrow = c(2, 2))
hist(Y, breaks = 100, probability = TRUE, main = "Unconditional distribution of the data", xlab = "Returns")
curve(dnorm(x, mean(Y), sd(Y)), add = TRUE, col = "darkgreen", lwd = 2)
acf(Y, lag.max = 100, main = "ACF data")
acf(Y^2, lag.max = 100, main = "ACF squared data")
qqnorm(Y, main = "Q-Q plot against Gaussian distribution"); qqline(Y)
par(mfrow = c(1, 1))

# Construct the same starting values and parameter bounds used in the MATLAB script.
if (I_s == 1) {
  if (lev_s == 1) {
    rS <- rep(kk_s, 5); rS[c(1, 2, 4)] <- c(mean(Y), log(sd(Y)), log(8))
    lower <- c(-Inf, -Inf, 0, -Inf, -Inf); upper <- rep(Inf, 5)
  } else {
    rS <- rep(kk_s, 4); rS[c(1, 2, 4)] <- c(mean(Y), log(sd(Y)), log(8))
    lower <- c(-Inf, -Inf, 0, -Inf); upper <- rep(Inf, 4)
  }
} else {
  if (lev_s == 1) {
    rS <- rep(kk_s, 6); rS[c(1, 2, 3, 5)] <- c(mean(Y), log(sd(Y)), .9, log(8))
    lower <- c(-Inf, -Inf, 0, 0, -Inf, -Inf); upper <- c(Inf, Inf, 1, Inf, Inf, Inf)
  } else {
    rS <- rep(kk_s, 5); rS[c(1, 2, 3, 5)] <- c(mean(Y), log(sd(Y)), .9, log(8))
    lower <- c(-Inf, -Inf, 0, 0, -Inf); upper <- c(Inf, Inf, 1, Inf, Inf)
  }
}

# L-BFGS-B is R's bounded optimizer counterpart to constrained MATLAB optimization.
f_s <- function(x) Beta_univ_t_Egarch_Logl(x, Y, Inf_s, I_s)
fit_optim <- optim(rS, f_s, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 67000))
x_s <- fit_optim$par
dir.create(results_dir, showWarnings = FALSE)
write.csv(data.frame(estimate = x_s), file.path(results_dir, "estPar_R.csv"), row.names = FALSE)

# Obtain filtered scale, residuals, scores, and information criteria.
fitted <- Beta_univ_t_Egarch_Logl_fact(x_s, Y, Inf_s, I_s)
criteria <- Info_Crit(fitted$Logl, x_s, t)
vega_s <- exp(if (I_s == 1) x_s[4] else x_s[5])
mu_s <- x_s[1]

# Plot fitted scale against returns and absolute returns.
par(mfrow = c(2, 1))
plot(Y - mu_s, type = "l", main = "Scale fit: returns", ylab = "Returns / scale")
lines(exp(fitted$lam), col = "red", lwd = 2)
legend("topright", c("Y - mu", "Fitted conditional scale"), col = c("black", "red"), lty = 1)
plot(abs(Y - mu_s), type = "l", main = "Scale fit: absolute returns", ylab = "Absolute returns / scale")
lines(exp(fitted$lam), col = "red", lwd = 2)
par(mfrow = c(1, 1))

# Compare standardized residuals with fitted Student t and Gaussian densities.
hist(fitted$res, breaks = 100, probability = TRUE, main = "Fitted residual distribution", xlab = "Standardized residuals")
grid <- seq(min(fitted$res), max(fitted$res), length.out = 500)
lines(grid, dt(grid, df = vega_s), col = "red", lwd = 2)
lines(grid, dnorm(grid), col = "darkgreen", lwd = 2)
legend("topright", c("Fitted t", "Standard Gaussian"), col = c("red", "darkgreen"), lty = 1)

# Residual, squared-residual, and score ACF diagnostics.
par(mfrow = c(3, 1))
acf(fitted$res, lag.max = 100, main = "ACF residuals")
acf(fitted$res^2, lag.max = 100, main = "ACF squared residuals")
acf(fitted$u, lag.max = 100, main = "ACF fitted scores")
par(mfrow = c(1, 1))

# The analytic Student t CDF is R's counterpart to the numerical MATLAB PIT integration.
PIT_s <- pt(sort(fitted$res), df = vega_s)
plot(PIT_s, type = "l", main = "Ordered PIT DCS model versus uniform", xlab = "Residual rank", ylab = "PIT")
lines(seq(0, 1, length.out = length(PIT_s)), lty = 2)
