# Survival analysis simulation: Cox, log-rank and RMST

A simulation study (ADEMP framework) comparing three ways to detect a treatment
difference in time-to-event data: the Cox hazard ratio, the log-rank test and the
RMST difference (restricted mean survival time up to 24 months).
Scenarios: proportional hazards, delayed effect, crossing curves, and a null
scenario, each with two censoring levels and three sample sizes
(24 scenarios x 1000 replications).

## Files

- `R/functions.R`: data generation, the three methods, true RMST
- `R/run_simulation.R`: runs all 24 scenarios and writes `results2.csv`
- `report.qmd`: figures and summary of the results

## How to rerun

1. Install the R packages: survival, dplyr, tidyr, ggplot2.
2. From the project folder, run the simulation: `Rscript R/run_simulation.R`
3. Render the report: `quarto render report.qmd`

The random seed is fixed (2026), so the results are reproducible.

## Main findings

- Type I error of all three methods is close to the nominal 5%.
- Under proportional and delayed effects, Cox and log-rank are more powerful than RMST.
- Under crossing curves, Cox and log-rank detect that the hazards differ, while RMST
  stays near 5% because the true 24-month RMST difference is tiny (0.097849).
- Heavy censoring (50%) lowers power for all methods.
- No method is best everywhere: they answer different questions.