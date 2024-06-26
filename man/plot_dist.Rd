% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_dist}
\alias{plot_dist}
\title{Plot a distribution}
\usage{
plot_dist(d)
}
\arguments{
\item{d}{List that defines the distribution (as returned by \code{def_dist_incubation_period()} for example)}
}
\value{
A ggplot object.
}
\description{
Plot a distribution
}
\examples{
# Define a `ern` distribution:
gi  = ern::def_dist(
  dist     = "gamma",
  mean     = 6.84,
  mean_sd  = 0.7486,
  shape    = 2.39,
  shape_sd = 0.3573,
  max      = 15
  )

# Plot can be customized like any `ggplot` object:
g = plot_dist(gi) + ggplot2::labs(subtitle = 'your subtitle')
plot(g)


}
