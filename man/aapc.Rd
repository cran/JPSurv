\name{aapc}
\alias{aapc}
\title{Trend summary measures for joinpoint relative survival model}
\description{Get the trend summary measures for joinpoint relative survival model quickly using default interval=5 (if available).  
A function that returns the trend summary measures including the Absolute Change in Survival (ACS), the Relative Change in Survival (RCS) .
and the Relative Change in Hazard (RCH) for each joinpoint segment}
\usage{aapc(fit, type="AbsChgSur",interval=5)}
\arguments{
\item{fit}{Joinpoint object with predicted values and joinpoint model selections.}
\item{type}{Type of trend summary measure.  
Supported measures are: 
RelChgHaz - Hazard of cancer death, 
AbsChgSur - Absolute change in survival.
RelChgSur - Relative change in survival. The default is AbsChgSur.}
\item{interval}{The number of years since diagnosis (follow-up years).}
}
\value{A data frame with columns \code{start.year}, \code{end.year}, \code{estimate},
      \code{std.error}, \code{lowCI}, \code{upCI} containing the 
estimates, standard errors, and confidence interval of the trend summary measure.}

\examples{
data("fit2", package="JPSurv")

# Get the estimate, standard error, and confidence interval of 
# the annual changes of cumulative relative survival for interval=5.
aapc(fit2, type="RelChgHaz")
}
