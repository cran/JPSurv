\name{aapc.multiints}
\alias{aapc.multiints}
\title{
Trend summary measures for multiple selected intervals
}
\description{
A function that returns the trend summary measures including the 
 Absolute Change in Survival (ACS), the Relative Change in Survival (RCS) 
and the Relative Change in Hazard (RCH) for each joinpoint segment given multiple selected intervals.
It can also return the Weighted Average Abosulte Change in Survival for the user-specified ACS year range.}
\usage{
aapc.multiints(fit, type = "AbsChgSur", int.select = NULL, 
               ACS.range = NULL, ACS.out = NULL)
}
\arguments{
  \item{fit}{
Joinpoint object with predicted values and joinpoint model selections.
}
  \item{type}{
Type of trend summary measure.  
Supported measures are: 
RelChgHaz - Hazard of cancer death, 
AbsChgSur - Cumulative Relative Survival Average Annual Absolute Change.
RelChgSur - Cumulative Relative Survival Average Annual Relative Change. The default is AbsChgSur.
}
  \item{int.select}{
The single or multiple interval values selected for the trend measures calculation. The default is NULL.
}
  \item{ACS.range}{
It is needed for type="AbsChgSur" only. The minimum and maximum of the year range specified for the Average Absolute Change in Survival trend measure. 
The default is NULL. If ACS.range=NULL, the average absolute change in survival for joinpoint segments will be returned. 
Otherwise, the weighted average absolute change in survival for the desired range of years along with the the average absolute change in survival 
for joinpoint segments will be produced. 
Note that, the ACS.range should be within the entire range of data. 
}
  \item{ACS.out}{
It is used to specify which ACS trend results will be returned when type="AbsChgSur". ACS.out can be defined as "user" or "both". 
If ACS.range=NULL, ACS.out should always be NULL and the ACS trend results between joinpoints will be returned. If ACS.range is the year range, the ACS trend measure results between the selected years will be returned when
ACS.out="user" and the ACS trend results for both between joinpoints and user-selected years will be returned when ACS.out="both". The default is ACS.out=NULL.
}
}
\value{ A list of data frames corresponding to each interval from \code{int.select}. 
Each data frame contains the columns \code{start.year}, \code{end.year}, \code{estimate},
      \code{std.error}, \code{lowCI}, \code{upCI} for the 
    estimates, standard errors, and confidence interval of the trend summary measure.
}
\author{
Fanni Zhang <zhangf@imsweb.com>
}
\seealso{
\code{\link{aapc}}
}
\examples{
data("fit1", package="JPSurv")

aapc.multiints(fit1, type="RelChgHaz",int.select=c(1,3,5)) 
}

