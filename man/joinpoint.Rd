\name{joinpoint}
\alias{joinpoint}
\alias{print.joinpoint}
\alias{plot.joinpoint}
\alias{summary.joinpoint}
\title{Fitting a join point relative survival model}
\description{Fitting a joinpoint relative survival model}
\usage{joinpoint(formula, data, subset, numJPoints = 0, ...)}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted.  If the input data is from SEER*Stat, the formula can be ~Year.  Otherwise, the full formula should be specified as: ~Alive_at_Start + Died + Lost_to_Followup + Expected_Survival_Interval + Interval + Year}
\item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{numJPoints}{The maximum number of join points allowed.  The default is zero, which is equivalent to a proportional hazard relative survival model.}
\item{\dots}{additional arguments to be passed to the low level regression fitting functions.}
}
\value{
An object of class "joinpoint" will be returned with attributes:
\item{coefficients}{a named vector of coefficients and standard errors}
\item{jp}{the estimates of the join points}
\item{converged}{convergence status}
\item{predicted}{the fitted relative survival rates}
\item{xbeta}{the linear predictor}
\item{ll}{log likelihood}
\item{aic}{AIC}
\item{bic}{BIC}
}
\examples{
#load the SEER-9 prostate cancer dataset.
data(prostate);
# Fit the survival join point model with zero join points, i.e., fit the proportional hazard relative survival model.
fit1 = joinpoint(~Year, data=prostate, numJPoints = 0);
# Fit the survival join point model on a subset of the dataset.
fit2 = joinpoint(~Year, data=prostate, Site_and_age_for_total_stage == 18 & Sex == 1, numJPoints = 0);
# Fit the survival join point model for non-SEER*Stat datasets.
fit3 = joinpoint(~Alive_at_Start + Died + Lost_to_Followup + Expected_Survival_Interval
		+ Interval + Year, prostate, numJPoints = 0);
plot(fit1);
}
\references{
Yu, B., Huang, L., Tiwari, R. C., Feuer, E. J. and Johnson, K. A. (2009), Modeling population-based cancer survival trends by using join point models for grouped survival data.  Journal of the Royal Statistical Soceity: Series A, 172, 405-425.
}
