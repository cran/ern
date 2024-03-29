% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agg_to_daily.R
\name{agg_to_daily}
\alias{agg_to_daily}
\title{Infer daily counts from aggregates}
\usage{
agg_to_daily(cl.input, dist.gi, popsize, prm.daily, silent = FALSE)
}
\arguments{
\item{cl.input}{Data frame. Must have variables:
\itemize{
\item \code{date}: calendar date of report
\item \code{value}: count of reported cases
}}

\item{dist.gi}{List. Parameters for the generation interval distribution in the same format as returned by \code{\link[=def_dist_generation_interval]{def_dist_generation_interval()}}.}

\item{popsize}{Integer. Population size to use in MCMC simulation to infer daily observations from aggregated input data.}

\item{prm.daily}{List. Parameters for daily report inference via MCMC. Elements include:
\itemize{
\item \code{burn}: Numeric. Length of burn-in period (number of days).
\item \code{iter}: Numeric. Number of iterations after burn-in period (number of days).
\item \code{chains}: Numeric. Number of chains to simulate.
\item \code{first.agg.period}: length of aggregation period for first aggregated observation (number of days); if NULL, assume same aggregation period as observed for second observation (gap between first and second observations)
\item \code{prior_R0_shape}: Shape of the (hyper-)parameter for the prior Gamma distribution for R0.
\item \code{prior_R0_rate}: Rate of the (hyper-)parameter for the prior Gamma distribution for R0.
\item \code{prior_alpha_shape}: Shape of the (hyper-)parameter for the prior Gamma distribution for alpha.
\item \code{prior_alpha_rate}: Rate of the (hyper-)parameter for the prior Gamma distribution for alpha.
}}

\item{silent}{Logical. Flag to suppress all output messages, warnings, and progress bars.}
}
\value{
Data frame with individual realizations of daily reported cases
}
\description{
Infer daily counts from aggregates
}
\examples{

# Importing data attached to the `ern` package
# and selecting the Omicron wave in Ontario, Canada.
# This is *weekly* incidence.
data(cl.input)
data = cl.input[cl.input$pt == 'on' & 
                  cl.input$date > as.Date('2021-11-30') & 
                  cl.input$date < as.Date('2021-12-31'),] 
head(data)
dist.gi = def_dist_generation_interval('sarscov2')

a = agg_to_daily(
cl.input = data, 
dist.gi = dist.gi, 
  popsize = 14e6, 
  prm.daily = list(
  # MCMC parameters.
  # small values for computation speed for this example.
  # Increase for better accuracy
  burn = 100,
  iter = 100,
  chains = 1,
  # - - - - - 
  prior_R0_shape = 2,
  prior_R0_rate = 0.6,
  prior_alpha_shape = 1,
  prior_alpha_rate = 1
))
# This is a Bayesian inference, so we 
# have a posterior distribution of  
# daily incidences. Here we just plot
# one single draw:
a1 = a[a$id==1,]
plot(x = a1$t, y = a1$value, typ = 'o',
     xlab = 'days', ylab = 'daily incidence',
     main = 'Posterior daily incidence infered from weekly incidence')
      

}
