\name{Plot.surv.int.multiyears}
\alias{Plot.surv.int.multiyears}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Figure - Cumulative Survival by Interval
}
\description{
A function that returns a plot for Cumulative Survival by Interval supporting multiple selected years.
}
\usage{
Plot.surv.int.multiyears(plotdata, fit, nJP, yearvar, 
                         obscumvar = "Relative_Survival_Cum", 
                         predcumvar = "Predicted_Survival_Cum", 
                         interval = "Interval", year.select = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plotdata}{
Either of the following data set would work.
1) the graph data returned by function download.data with downloadtype="graph" and all intervals needed for int.col.
2) the full data returned by function download.data with downloadtype="full". 
}
  \item{fit}{
The joinpoint object containing the model output.
}
  \item{nJP}{
The number of joinpoints in the model.
}
  \item{yearvar}{
The variable name for year of diagnosis used in argument 'year' of the function joinpoint.
}
  \item{obscumvar}{
The variable name for observed relative cumulative survival. The default is "Relative_Survival_Cum" for relative survival data. 
For cause-specific data, it needs to be changed accordingly.
}
  \item{predcumvar}{
The variable name for predicted cumulative survival. The default is "Predicted_Survival_Cum".
}
  \item{interval}{
The variable name for year since diagnosis. The default is 'Interval'.
}
  \item{year.select}{
The year values selected for the plot. The default is NULL.
}
}

\author{
Fanni Zhang <zhangf@imsweb.com>
}
\value{An object of class \code{ggplot} containing the plot.}
\seealso{
 \code{\link{download.data}},\code{\link{Plot.dying.year.annotate}}, 
 \code{\link{Plot.surv.year.annotate}}.
}
\examples{
data("breast.example", package="JPSurv")
data("fit3", package="JPSurv")

yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
interval<-"Interval"
nJP<-2
data.full<-download.data(breast.example,fit3,nJP,yearvar,"full",
                         subsetStr3,interval="Interval")
plot<-Plot.surv.int.multiyears(data.full,fit3,nJP,yearvar,obscumvar,
                    predcumvar,interval,year.select=c(1985,1990))
}
