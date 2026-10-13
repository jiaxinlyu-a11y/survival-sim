source("R/functions.R")

# All scenarios: 4 shapes x 2 censoring levels x 3 sample sizes = 24 rows
grid <- expand.grid(
  shape = c("null", "proportional", "delayed", "crossing"),
  cens  = c(0.2, 0.5),
  n     = c(100, 300, 600),
  stringsAsFactors = FALSE
)

# Run one scenario nsim times, return one row of results
sim_cell <- function(shape, cens, n, nsim) {
  hr <- hr_table[[shape]]
  res <- replicate(nsim, {
    d <- gen_data2(n, 0.05, hr1 = hr[1], hr2 = hr[2], tau = 6,
                   cens_rate = cens_rate_for(cens))
    c(run_methods(d), rmst_diff(d))
  })
  truth <- rmst_true(0.05 * hr[1], 0.05 * hr[2]) - rmst_true(0.05, 0.05)
  
  rej_cox  <- mean(res["p_cox", ] < 0.05)
  rej_lr   <- mean(res["p_logrank", ] < 0.05)
  rej_rmst <- mean(res["p_rmst", ] < 0.05)
  
  data.frame(shape = shape, cens = cens, n = n, nsim = nsim,
             true_rmst_diff = truth,
             mean_rmst_diff = mean(res["rmst_diff", ]),
             bias           = mean(res["rmst_diff", ]) - truth,
             mcse_bias      = sd(res["rmst_diff", ]) / sqrt(nsim),
             reject_cox     = rej_cox,
             reject_logrank = rej_lr,
             reject_rmst    = rej_rmst,
             mcse_cox       = mc_se(rej_cox, nsim),
             mcse_logrank   = mc_se(rej_lr, nsim),
             mcse_rmst      = mc_se(rej_rmst, nsim))
}

set.seed(2026)
nsim <- 1000
results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  sim_cell(grid$shape[i], grid$cens[i], grid$n[i], nsim)
}))
write.csv(results, "results2.csv", row.names = FALSE)