\name{Plot.surv.year.annotate}
\alias{Plot.surv.year.annotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Figure - Average Change in Cumulative Survival by Diagnosis Year
}
\description{
A function that returns a plot (or a list including the Absolute Change in Survival trend measures and the plot) for average change in cumulative survival by diagnosis year using ggplot.
The annotation feature is available for nJP<=3 and the number of multiple intervals selected <=3.
}
\usage{
Plot.surv.year.annotate(plotdata, fit, nJP, yearvar, 
                        obscumvar = "Relative_Survival_Cum", 
                        predcumvar = "Predicted_Survival_Cum", 
                        interval = "Interval", annotation = 0, 
                        trend = 0, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plotdata}{
The graph data returned by function download.data with downloadtype="graph".
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
  \item{annotation}{
The indicator for the annotation feature. The default is 0 (no annotation on the plot).Two plots with and without annotations will be returned in a list when annotation=1.
}
  \item{trend}{
The indicator for returning the Absolute Change in Survival trend measure tables. The default is 0 (no trend tables returned).
}
  \item{title}{
The title used for the plot. If title=NULL, the default title "Relative Survival by Diagonosis Year" or "Cause-Specific Survival by Diagnosis Year" will be used
according to the input data type.
}
}

\value{
  The returned object is a list and depends on the values for \code{trend} and \code{annotation}. 
  The list will contain some or all of the names \code{plot_anno}, \code{plot_no_anno}, and \code{trends},
  where \code{plot_anno} is an object of class \code{ggplot} for a plot with annotations,
  \code{plot_no_anno} is an object of class \code{ggplot} for a plot without annotations, and
  \code{trends} is a list of data frames for the change in hazard trend measures with columns
  \code{start.year}, \code{end.year}, \code{estimate}, \code{std.error}, 
  \code{lowCI}, \code{upCI}, and \code{interval}.
}

\author{
Fanni Zhang <zhangf@imsweb.com>
}

\seealso{
  \code{\link{download.data}}, \code{\link{Plot.dying.year.annotate}}, 
  \code{\link{Plot.surv.int.multiyears}},\code{\link{aapc.multiints}}.
}
\examples{
data("breast.example", package="JPSurv")
data("fit3", package="JPSurv")

yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
interval<-"Interval"
nJP<-2
data.graph<-download.data(breast.example,fit3,nJP,yearvar,"graph",
                          subsetStr3,interval,int.select=c(1,5))

out<-Plot.surv.year.annotate(data.graph,fit3,nJP,yearvar,obscumvar,
                              predcumvar,interval,annotation=0,trend=0)

}
