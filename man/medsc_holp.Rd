% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medsc_holp.R
\name{medsc_holp}
\alias{medsc_holp}
\title{Screen mediators by Ridge-HOLP method}
\usage{
medsc_holp(y, x, M, COV.S = NULL, d = NULL, r = 1)
}
\arguments{
\item{y}{n-dimensional vector of continuous outcomes.}

\item{x}{n-dimensional vector of exposures}

\item{M}{n by p matrix of mediators}

\item{COV.S}{a \code{data.frame} or \code{matrix} of covariates}

\item{d}{desired number of chosen mediators. Default value is \eqn{d = \cdot n/\log(n)}.}

\item{r}{scalar for the penalty parameter. Default value is 1.}
}
\value{
A vector for indexes of  the selected mediators.
}
\description{
Screen mediators by Ridge-HOLP method
}
\examples{
data(ExampleData) # Load ExampleData from the data from package
y <- ExampleData$y #   'y' is a vector of continuous outcomes in ExampleData
x <- ExampleData$x #   'x' is a vector of exposures in ExampleData
M <- ExampleData$M #   'M' is a matrix  of mediators in ExampleData

chosen_ind <- medsc_holp(y, x, M)
print(chosen_ind) # print indexes for the selected mediators

}
