% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulated_dataset.R
\name{simulated_dataset}
\alias{simulated_dataset}
\title{Simulated Dataset Generator}
\usage{
simulated_dataset(n, problem = c("meta-analysis", "multiple outcomes"))
}
\arguments{
\item{n}{Sample size for each study/outcome.}

\item{problem}{A character string specified what problem the simulated dataset is
generated for. \code{problem} can be set to "meta-analysis" or "multiple outcomes".}
}
\value{
A list object of the ingredients from the simulated dataset. The elements of
this list depends on value of \code{problem}.

For \code{problem = "meta-analysis"},
\item{Xlist}{a list object with \eqn{k}th element denoting the baseline covariate matrix of \eqn{k}th study }
\item{Ylist}{a list object with \eqn{k}th element denoting the response vector of \eqn{k}th study}
\item{Trtlist}{a list object with \eqn{k}th element denoting the treatment vector of \eqn{k}th study and coded as 0 or 1}
\item{B}{the coefficient matrix containing \eqn{\delta_0}, \eqn{\bm{\delta}},
\eqn{\theta_0} and \eqn{\bm{\theta}}}

For \code{problem = "multiple outcomes"},
\item{X}{a matrix object denoting the baseline covariate matrix}
\item{Ylist}{a list object with \eqn{k}th element denoting the response vector of \eqn{k}th outcome}
\item{Trt}{a vector denoting the treatment and coded as 0 or 1)}
\item{B}{the coefficient matrix containing \eqn{\delta_0}, \eqn{\bm{\delta}},
\eqn{\theta_0} and \eqn{\bm{\theta}}}
}
\description{
Generate a simulated dataset, which could be used to demonstrate the
features of the mpersonalized package.
}
\details{
In the simulated dataset, outcomes are generated from the model
\deqn{Y = \delta_0 + \bm{X} \bm{\delta} + A (\theta_0 + \bm{X}\bm{\theta})+\epsilon,}
where \eqn{\bm{X}} is the baseline covariates and \eqn{A} is the treatment indicator coded
as 0,1. For different outcomes or studies, values of \eqn{\delta_0}, \eqn{\bm{\delta}},
\eqn{\theta_0} and \eqn{\bm{\theta}} are also different so as to represent the heterogeneity
in real problems.

The number of different studies/outcomes is set to be 6 and total number of candidate
covariates is 50. Treatment indicator \eqn{A} is generated with equal probability of
0 or 1.

This function randomly generates the coefficients for each study/outcome and then
generates the baseline covariates and error term for each subject. Depending on the
value of \code{problem}, generation of baseline covariates are slightly different.
For \code{problem = "meta-analysis"}, baseline covariates are generated independently
for each study; for \code{problem = "multiple outcomes"}, baseline covariates are
the same across different outcomes.
}
\examples{
set.seed(123)
sim_dat  = simulated_dataset(n = 200, problem = "meta-analysis")
str(sim_dat$Xlist)
str(sim_dat$Ylist)
str(sim_dat$Trtlist)
set.seed(NULL)

}
