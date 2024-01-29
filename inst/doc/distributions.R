## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE, message = FALSE, warning = FALSE--------------------
library(ern)

library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)

## ----manica-table, echo = FALSE, fig.align="center", fig.cap='**Table 2 from (Manica _et al._ 2022)**. Original caption: "Estimates for the incubation period, diagnostic delay, intrinsic and realized generation time, and household serial intervals of the SARS-CoV-2 Omicron variant. Reported parameters of shape and scale for the incubation period and intrinsic generation time refer to a gamma distribution. Estimates of the incubation period are dervied from the analysis of 80 participants to a single superspreading event in Norway. Data taken from Brandal _et al_."'----
knitr::include_graphics("figs/manica-table.png",
                        dpi = 144)

## ----optim--------------------------------------------------------------------
# Fit a standard deviation for a gamma distribution
# given a mean and 95% quantiles
fit_sd <- function(mean, ci){
  
  # Score function for optimization
  # s = shape
  # m = assumed mean
  # ci_target = the 95% credible interval bounds, assumed to be 2.5% and 97.5% quantiles
  score_fun <- function(s, m, ci_target) {
    shape <- s
    scale <- m/s
    ci <- qgamma(c(0.025, 0.975), shape = shape, scale = scale)
    return(sum((ci-ci_target)^2))
  }

  # Optimize over shape
  opt_out <- optim(
    par = 1, 
    fn = score_fun,
    m = mean,
    ci_target = ci,
    lower = 0.0001, upper = 1000,
    method = "Brent"
  )
  
  # Get best shape parameter found and return associated standard deviation
  shape_best <- opt_out$par
  sd_best <- shape_to_sd(shape_best, mean)
  
  return(list(
    sd = sd_best,
    score = opt_out$value
  ))
}

# Helper function to compute a gamma's standard deviation given the shape and mean
shape_to_sd <- function(shape, mean) mean/sqrt(shape)

## ----make_fit_table, echo = FALSE---------------------------------------------
#' Make table summarizing parameters for each family of distributions
#'
#' @param mean.inputs list with elements `mean` (scalar) and `ci` (numeric vector of length 2) 
#' @param shape.inputs same format as `mean.inputs`
#'
#' @return tibble 
make_summary_table <- function(mean.inputs, shape.inputs){

  fit.mean <- fit_sd(mean = mean.inputs$mean, ci = mean.inputs$ci)
  fit.shape <- fit_sd(mean = shape.inputs$mean, ci = shape.inputs$ci)
  
  tb <- tribble(
    ~parameter,           ~`value type`,        ~origin,    ~`value`,
    "distribution mean",  "mean",               "assumed",  mean.inputs$mean,
    NA,                   "standard deviation", "fitted",   fit.mean$sd,
    "distribution shape", "mean",               "assumed",  shape.inputs$mean,
    NA,                   "standard deviation", "fitted",   fit.shape$sd,
  )
}

## ----display_table, echo = FALSE----------------------------------------------
display_table <- function(tb, caption = "", digits = 4){
  knitr::kable(
    tb |> mutate(
    across(where(is.character), ~replace_na(.x, ""))
  ), 
    digits = digits,
    caption = caption
  )
}

## ----table-gen-int, echo = FALSE----------------------------------------------
summary.gen.int <- make_summary_table(
  mean.inputs = list(
    mean = 6.84,
    ci = c(5.72, 8.60)
  ),
  shape.inputs = list(
    mean = 2.39,
    ci = c(2.01, 3.34)
  )
)

display_table(summary.gen.int,
              caption = "Parameters for the family of generation interval distributions")

## ----table-inc-per, echo = FALSE----------------------------------------------
summary.inc.per <- make_summary_table(
  mean.inputs = list(
    mean = 3.49,
    ci = c(3.19, 3.77)
  ),
  shape.inputs = list(
    mean = 8.50,
    ci = c(6.14, 13.20)
  )
)

display_table(summary.inc.per,
              caption = "Parameters for the family of incubation period distributions")

## ----plot_mean_dist, echo = FALSE---------------------------------------------
plot_mean_dist <- function(summary, max, title){
  summary <- fill(summary, parameter)
  mean = (summary 
    |> filter(parameter == "distribution mean",
               `value type` == "mean")
    |> pull(value)
  )
  shape = (summary 
    |> filter(parameter == "distribution shape",
               `value type` == "mean")
    |> pull(value)
  )
  
  y <- get_discrete_dist(list(
    dist = "gamma",
    mean = mean,
    mean_sd = NA,
    shape = shape,
    shape_sd = NA,
    max = max
  ))
  
  df <- tibble(x = 1:max, y = y)
  
  (ggplot(df,
          aes(x = x, y = y))
    + geom_col()
    + scale_x_continuous(breaks = 1:max,
                         expand = expansion(mult = 0))
    + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    + labs(x = "days",
           title = title)
    + theme_bw()
    + theme(
      axis.title.y = element_blank(),
      text = element_text(size = 12)
    )
  )
}

## ----plot-mean-gen-int, echo = FALSE, fig.width = 4, fig.height = 2.5, fig.cap = "Generation interval distribution for the mean value of the mean and shape parameter"----
plot_mean_dist(summary.gen.int, max = 21,
               title = "Mean generation interval distribution")

## ----plot-mean-inc-per, echo = FALSE, fig.width = 4, fig.height = 2.5, fig.cap = "Incubation period distribution for the mean value of the mean and shape parameter"----
plot_mean_dist(summary.inc.per, max = 8, 
               title = "Mean incubation period distribution")

## ----echo=FALSE, fig.width = 6, fig.height = 4--------------------------------
shed.i = data.frame(
  x = seq(from = 1, to = 12, by = 2),
  y = c(5, 7.30103, 7.083333, 6.677121, 6.30, 5.9)
)

plot(x = shed.i$x,
     y = shed.i$y,
     "b",
     main = "Fecal shedding kinetics: I compartment.\nMean days of shedding: 12",
     xlab = "days",
     ylab = "log viral load")

shed.a = data.frame(
  x = seq(from = 1, to = 10, by = 2),
  y = c(5,7.30103,7.083333,6.677121, 6.30)
)
plot(x = shed.a$x,
     y = shed.a$y,
     "b",
     main = "Fecal shedding kinetics: A compartment.\nMean days of shedding: 10",
     xlab = "days",
     ylab = "log viral load")

shed.j = data.frame(
  x = seq(from = 1, to = 8, by = 2),
  y = c(5,7.30103,7.083333,6.677121)
)
plot(x = shed.j$x,
     y = shed.j$y,
     "b",
     main = "Fecal shedding kinetics: J compartment.\nMean days of shedding: 8",
     xlab = "days",
     ylab = "log viral load")

shed.z = data.frame(
  x = seq(from = 1, to = 24, by = 4),
  y = c(5.40103,4.60103,3.90103,3.10103,2.10103,1.2)
)
plot(x = shed.z$x,
     y = shed.z$y,
     "b",
     main = "Fecal shedding kinetics: Z compartment.\nMean days of shedding: 24",
     xlab = "days",
     ylab = "log viral load")

## -----------------------------------------------------------------------------
shed.symp = data.frame(
  x = c(seq(from = 1, to = 12, by = 2),seq(from = 13, to = 36, by = 4)),
  y = c(5, 7.30103, 7.083333, 6.677121, 6.30, 5.9, 5.40103,4.60103,3.90103,3.10103,2.10103,1.2)
)

shed.asymp = data.frame(
  x = c(seq(from = 1, to = 10, by = 2),seq(from = 11, to = 34, by = 4)),
  y = c(5, 7.30103, 7.083333, 6.677121, 6.30, 5.40103,4.60103,3.90103,3.10103,2.10103,1.2)
)

shed.severe = data.frame(
  x = seq(from = 1, to = 8, by = 2),
  y = c(5,7.30103,7.083333,6.677121)
)



## -----------------------------------------------------------------------------
fec.l = list(
  shed.symp,
  shed.asymp,
  shed.severe
)

df = dplyr::bind_rows(fec.l) |>
  dplyr::group_by(x) |>
  dplyr::summarise(y = mean(y, na.rm=TRUE))


## ----echo=FALSE, fig.width = 6, fig.height = 4--------------------------------
plot(x = df$x,
     y = df$y,
     main = "Average fecal shedding kinetics distribution",
     xlab = "days",
     ylab = "log viral load", typ='b')

## ----warning=FALSE, message=FALSE---------------------------------------------
get_gamma_fec <- function(df){
  log_likelihood = function(mean,shape){
    a = shape
    s = mean/shape
    return(-sum(df$y*(dgamma(df$x, shape = a, scale = s, log = TRUE))))
  }

  mle.res = stats4::mle(minuslogl = log_likelihood,
                                    start = list(mean = 14, shape = 2),
                                    nobs = nrow(df))
  
  res = stats4::summary(mle.res)
  ci = stats4::confint(object = mle.res, level = 0.68)
  
  f = list(
    dist = "gamma",
    mean  = res@coef[[1]],
    shape = res@coef[[2]],
    # we assume this parameter will be sampled with 
    # a normal distribution (hence 68%CI ~ 2 std dev):
    mean_sd  = (ci[1,2] - ci[1,1])/2, 
    shape_sd = (ci[2,2] - ci[2,1])/2,
    max = max(df$x)
  )

  return(f)
}

fec = get_gamma_fec(df)

fec

## ----echo=FALSE, fig.width = 6, fig.height = 4--------------------------------
y = dgamma(df$x, shape = fec$shape, scale = fec$mean/fec$shape)
plot(x = df$x,
     y = y/sum(y), 
     'l',
     main = "Fitted fecal shedding kinetics distribution",
     xlab = "days",
     ylab = "log viral load (normalized)")
points(x = df$x, y = df$y/sum(df$y))

## ----echo=FALSE---------------------------------------------------------------
summary.fec.shed <- tribble(
  ~parameter,             ~`value type`,        ~origin,  ~`value`,
    "distribution mean",  "mean",               "fitted", fec$mean,
    NA,                   "standard deviation", "fitted", fec$mean_sd,
    "distribution shape", "mean",               "fitted", fec$shape,
    NA,                   "standard deviation", "fitted", fec$shape_sd,
)

display_table(summary.fec.shed,
              caption = "Parameters for the family of fecal shedding distributions")

## ----display_table_all, echo = FALSE------------------------------------------
display_table_all <- function(df.list, caption){
  df <- (bind_rows(lapply(names(df.list), function(x){
    df <- (df.list[[x]] |> mutate(distribution = NA) |> relocate(distribution, .before = everything()))
    df$distribution[1] <- x
    
    df
  }))
    |> select(-origin)
  )
  
  display_table(df, caption = caption)
}

## ----table-all, echo = FALSE--------------------------------------------------
df.list <- list(
  `generation interval` = summary.gen.int,
  `incubation period` = summary.inc.per,
  `fecal shedding` = summary.fec.shed
)

display_table_all(
  df.list = df.list,
  caption = "Summary of parameters for each family of distributions as implemented in `ern`"
)

