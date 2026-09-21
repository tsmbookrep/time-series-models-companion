# Information criteria used by the DCSt model.

Info_Crit <- function(Logl, par, n) {
  BIC <- -2 * Logl + length(par) * log(n)
  AIC <- 2 * length(par) - 2 * Logl
  list(BIC = BIC, AIC = AIC)
}
