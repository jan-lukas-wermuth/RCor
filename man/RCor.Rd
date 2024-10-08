% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCor.R
\name{RCor}
\alias{RCor}
\title{Kendall, Spearman and Pearson correlation and their generalizations for non-continuous data}
\usage{
RCor(X, Y, alpha = 0.1, method = "gamma", IID = TRUE, Fisher = TRUE)
}
\arguments{
\item{X}{a n x 1 numeric vector, matrix or data frame.}

\item{Y}{a n x 1 numeric vector, matrix or data frame.}

\item{alpha}{a numeric value specifying the significance level. The confidence level will be 1 - alpha.}

\item{method}{a character string specifying the correlation coefficient to be used for the independence test. Possible values are "tau", "tau_b", "tau_b_mod", "gamma", "rho", "rho_b" and "r". The recommendation for data with ties is "gamma". Specifying "tau_b_mod" only yields the independence test for IID data.}

\item{IID}{logical indicator determining whether the inference shall be conducted under iid (default) or time series assumptions (see CITATION for a precise description of the assumptions)}

\item{Fisher}{logical indicator determining whether the CIs shall be computed by using the Fisher transformation.}
}
\value{
The value of the chosen correlation coefficient along with its confidence interval and an independence test as well as an uncorrelatedness test.
}
\description{
\code{RCor()} computes the specified correlation with corresponding confidence intervals and P-values for the associated independence test either in the iid or in the time series case.
}
\examples{
X <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
Y <- c(1, 1, 1, 1, 2, 1, 2, 2, 2)
RCor(X, Y)
}
