\name{CoxModel_Year}
\alias{CoxModel_Year}
\title{Fitting a proportional hazard relative survival model with year as a covariate}
\description{Fitting a proportional hazard relative survival model with year as a covariate}
\usage{CoxModel_Year(formula, data, subset, ...)}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted.  If the input data is from SEER*Stat, the formula can be ~Year.  Otherwise, the full formula should be specified as: ~Alive_at_Start + Died + Lost_to_Followup + Expected_Survival_Interval + Interval + Year}
\item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{\dots}{additional arguments to be passed to the low level regression fitting functions.}
}
\value{
A list with attributes:
\item{coefficients}{a named vector of coefficients and standard errors}
\item{converged}{convergence status}
\item{predicted}{the fitted relative survival rates}
\item{xbeta}{the linear predictor}
\item{ll}{log likelihood}
\item{aic}{AIC}
\item{bic}{BIC}
}
