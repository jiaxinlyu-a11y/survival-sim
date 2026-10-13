library(survival)

# Generate data: event time = piecewise exponential (hazard h1 before tau, h2 after)
gen_data2 <- function(n, rate_old, hr1, hr2, tau, cens_rate) {
  group <- rep(c("old", "new"), each = n)
  h1 <- ifelse(group == "old", rate_old, rate_old * hr1)   # hazard before tau
  h2 <- ifelse(group == "old", rate_old, rate_old * hr2)   # hazard after tau
  
  e <- rexp(2 * n, 1)    # random "risk budget" for each patient
  t_event <- ifelse(e <= h1 * tau, e / h1, tau + (e - h1 * tau) / h2)
  
  t_cens <- rexp(2 * n, cens_rate)
  data.frame(
    group = factor(group, levels = c("old", "new")),
    time  = pmin(t_event, t_cens),
    event = as.integer(t_event <= t_cens)
  )
}

# Cox HR and p-value, log-rank p-value
run_methods <- function(d) {
  cox <- summary(coxph(Surv(time, event) ~ group, data = d))
  lr  <- survdiff(Surv(time, event) ~ group, data = d)
  c(HR        = as.numeric(cox$coefficients[1, "exp(coef)"]),
    p_cox     = as.numeric(cox$coefficients[1, "Pr(>|z|)"]),
    p_logrank = 1 - pchisq(lr$chisq, df = 1))
}

# RMST difference (new - old) up to tau, and its p-value
rmst_diff <- function(d, tau = 24) {
  fit  <- survfit(Surv(time, event) ~ group, data = d)
  tb   <- summary(fit, rmean = tau)$table
  diff <- tb["group=new", "rmean"] - tb["group=old", "rmean"]
  se   <- sqrt(tb["group=new", "se(rmean)"]^2 + tb["group=old", "se(rmean)"]^2)
  c(rmst_diff = diff, p_rmst = 2 * (1 - pnorm(abs(diff / se))))
}

# True survival curve: hazard h1 before tau, h2 after
S <- function(t, h1, h2, tau = 6) {
  ifelse(t <= tau, exp(-h1 * t), exp(-h1 * tau - h2 * (t - tau)))
}

# True RMST: area under the true curve from 0 to L
rmst_true <- function(h1, h2, L = 24) {
  curve_fn <- function(t) S(t, h1, h2)
  integrate(curve_fn, 0, L)$value
}

# Censoring rate that makes a target share of the old group censored
cens_rate_for <- function(target, rate_old = 0.05) {
  target * rate_old / (1 - target)
}

# Monte Carlo standard error of a proportion
mc_se <- function(p, nsim) sqrt(p * (1 - p) / nsim)

# hr1 (before month 6) and hr2 (after month 6) for each scenario
hr_table <- list(
  null         = c(1, 1),
  proportional = c(0.7, 0.7),
  delayed      = c(1, 0.5),
  crossing     = c(1.5, 0.5)
)