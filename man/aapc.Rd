\name{aapc}
\alias{aapc}
\title{Trend summary measures for joint point relative survival model}
\description{Get the trend summary measures for joint point relative survival model.  Meausres include annual percentage changes of hazard, annual percentage changes of cumulative relative survival, annual changes of cumulative relative survival.}
\usage{aapc(fit, type="HAZ_AC(CS)", interval=5)}
\arguments{
\item{fit}{Object of class "joinpoint".}
\item{type}{Type of trend summary measure.  Supported measures are: HAZ_APC(HR) - annual percentage changes of hazard, HAZ_APC(CS) - annual percentage changes of cumulative relative surviva, HAZ_AC(CS) - annual changes of cumulative relative survival.  The default is HAZ_AC(CS).}
\item{interval}{Years after diagnosis.  Only needed for HAZ_APC(CS) and HAZ_AC(CS).}
}
\value{The estimates and standard erros of the trend summary measure.}
\examples{
data(prostate);
# Fit the survival join point model with zero join points, i.e., fit the proportional hazard relative survival model.
fit1 = joinpoint(~Year, data=prostate, numJPoints = 0);
# Get the estimates and standard errors of the annual changes of cumulative relative survival.
haz_ac = aapc(fit1, type="HAZ_AC(CS)", interval=5);
}
