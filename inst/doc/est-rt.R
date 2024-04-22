## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE, warning=FALSE-----------------------------------------
suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(patchwork)
  library(ern)
})

## ----load sample data---------------------------------------------------------
# Loading sample SARS-CoV-2 wastewater data
ww.conc = ern::ww.data

head(ww.conc)

## ----define fec shed and gi---------------------------------------------------
# Define SARS-CoV-2 fecal shedding and generation interval distributions
dist.fec = ern::def_dist(
  dist = "gamma",
  mean = 12.90215,
  mean_sd = 1.136829,
  shape = 1.759937,
  shape_sd = 0.2665988,
  max = 33
)

dist.gi  = ern::def_dist(
  dist     = "gamma",
  mean     = 6.84,
  mean_sd  = 0.7486,
  shape    = 2.39,
  shape_sd = 0.3573,
  max      = 15
)

## ----plot dist fecal, fig.width = 6-------------------------------------------
plot_dist(dist.fec) + labs(title = paste0("Mean fecal shedding distribution (", dist.fec$dist, ")"))

## ----plot dist gi, fig.width = 6----------------------------------------------
plot_dist(dist.gi) + labs(title = paste0("Mean generation interval distribution (", dist.gi$dist, ")"))

## ----initializing smoothing and Rt params-------------------------------------
# Initializing scaling factor
scaling.factor = 1

# Initializing smoothing parameters
prm.smooth = list(
  align  = 'center', # smoothing alignment
  method = 'loess',  # smoothing method
  span   = 0.30,     # smoothing span (used for loess smoothing only)
  floor  = 5         # minimum smoothed concentration value (optional, loess smoothing only)
)

# Initialzing Rt settings
prm.R = list(
  iter   = 20,   # number of iterations in Rt ensemble
  CI     = 0.95, # confidence interval
  window = 10,   # Time window for Rt calculations
  config.EpiEstim = NULL # optional EpiEstim configuration for Rt calculations
)

## ----estimate rt--------------------------------------------------------------
r.estim = ern::estimate_R_ww(
  ww.conc        = ww.conc,
  dist.fec       = dist.fec,
  dist.gi        = dist.gi,
  scaling.factor = scaling.factor,
  prm.smooth     = prm.smooth,
  prm.R = prm.R,
  silent = TRUE # suppress output messages
)

## ----plot rt, fig.width = 6, fig.height = 6, warning = FALSE------------------
g = ern::plot_diagnostic_ww(r.estim)
plot(g)

## ----load-data----------------------------------------------------------------
dat <- (ern::cl.data 
    |> dplyr::filter(
      pt == "qc",
      dplyr::between(date, as.Date("2021-07-01"), as.Date("2021-09-01"))
))

## ----load-dists---------------------------------------------------------------
# define reporting delay 
dist.repdelay = ern::def_dist(
  dist = 'gamma',
  mean = 5,
  mean_sd = 1,
  sd = 1,
  sd_sd = 0.1,
  max = 10
)

# define reporting fraction 
dist.repfrac = ern::def_dist(
  dist = "unif",
  min = 0.1,
  max = 0.3
)

# define incubation period
dist.incub = ern::def_dist(
  dist     = "gamma",
  mean     = 3.49,
  mean_sd  = 0.1477,
  shape    = 8.5,
  shape_sd = 1.8945,
  max      = 8
)

# define generation interval
dist.gi = ern::def_dist(
  dist     = "gamma",
  mean     = 6.84,
  mean_sd  = 0.7486,
  shape    = 2.39,
  shape_sd = 0.3573,
  max      = 15
)

## ----plot-dist-repdelay, fig.width = 6----------------------------------------
plot_dist(dist.repdelay) + labs(title = paste0("Mean reporting delay distribution (", dist.repdelay$dist, ")"))

## ----plot-dist-incub, fig.width = 6-------------------------------------------
plot_dist(dist.incub) + labs(title = paste0("Mean incubation period distribution (", dist.incub$dist, ")"))

## ----plot-dist-gi, fig.width = 6----------------------------------------------
plot_dist(dist.gi) + labs(title = paste0("Mean generation interval distribution (", dist.gi$dist, ")"))

## ----init-prms----------------------------------------------------------------
# settings for daily report inference
prm.daily = list(
  method = "renewal",
  popsize = 1e7, # population size
   # Here, low value for `burn` and `iter` 
   # to have a fast compilation of the vignette.
   # For real-world applications, both `burn` and `iter`
   # should be significantly increased (e.g., 10,000).
   # Also, the number of chains should be at least 3 
   # (instead of 1 here) for real-world applications.
  burn = 100,
  iter = 200,
  chains = 1,
  prior_R0_shape = 2,
  prior_R0_rate = 0.6,
  prior_alpha_shape = 1,
  prior_alpha_rate = 1
)

# settings for checks of daily inferred reports
prm.daily.check = list(
  agg.reldiff.tol = 200
)

# smoothing settings for daily inferred reports
prm.smooth = list(
  method = "rollmean",
  window = 3,
  align = 'center'
)

# Rt settings
prm.R = list(
  iter = 10, # number of iterations in Rt ensemble
  CI = 0.95, # 95% confidence interval
  window = 7, # time window for each Rt estimate
  config.EpiEstim = NULL
)

## ----est-rt-------------------------------------------------------------------
r.estim = estimate_R_cl(
  cl.data      = dat,
  dist.repdelay = dist.repdelay,
  dist.repfrac  = dist.repfrac,
  dist.incub    = dist.incub,
  dist.gi       = dist.gi,
  prm.daily     = prm.daily,
  prm.daily.check = prm.daily.check,
  prm.smooth    = prm.smooth,
  prm.R         = prm.R,
  silent        = TRUE # suppress output messages
)

## ----plot-rt, fig.width = 6, fig.height = 6, warning = FALSE------------------
g = plot_diagnostic_cl(r.estim)
plot(g)

