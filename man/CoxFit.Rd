\name{CoxFit}
\alias{CoxFit}
\title{Fitting a proportional hazard relative survival model}
\description{Fitting a proportional hazard relative survival model}
\usage{CoxFit(X, nAlive, nDied, nLost, expSurv)}
\arguments{
\item{X}{The input design matrix.}
\item{nAlive}{Number of people at risk.}
\item{nDied}{Number of people who died.}
\item{nLost}{Number of people who are censored.}
\item{expSurv}{The expected survival rate.}
}
\value{A CoxFit class object.}
