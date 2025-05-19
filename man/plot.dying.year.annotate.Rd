\name{Plot.dying.year.annotate}
\alias{Plot.dying.year.annotate}
\title{
Figure - Percent Change in the Annual Probability of Dying of Cancer by Diagnosis Year
}
\description{
A function that returns a plot (or a list including the Relative Change in Hazard trend measures and the plot) for percent change in the annual probability of dying of cancer by diagnosis year using ggplot. 
The annotation feature is available for nJP<=3 and the number of multiple intervals selected <=3. 
}
\usage{
Plot.dying.year.annotate(plotdata, fit, nJP, yearvar, 
                         obsintvar = "Relative_Survival_Interval", 
                         predintvar = "Predicted_ProbDeath_Int", 
                         interval = "Interval", 
                         annotation = 0, topanno = 1, trend = 0,
                         title = NULL)
}
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
  \item{obsintvar}{
The variable name for observed interval survival. The default is "Relative_Survival_Interval" for relative survival data. 
For cause-specific data, it needs to be changed accordingly.
}
  \item{predintvar}{
The variable name for predicted interval survival. The default is "Predicted_ProbDeath_Int".
}
  \item{interval}{
The variable name for year since diagnosis. The default is 'Interval'.
}
  \item{annotation}{
The indicator for the annotation feature. The default is 0 (no annotation on the plot).Two plots with and without annotations will be returned in a list when annotation=1.
}
  \item{topanno}{
The indicator for showing the top curve annotation. The default is 1 (annotation for the top curve).
}
  \item{trend}{
The indicator for returning the Relative Change in Hazard trend measures . The default is 0 (no trend tables returned).
}
  \item{title}{
The title used for the plot. If title=NULL, the default title "Annual Probability of Dying of Cancer by Diagnosis Year" will be used.
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
 \code{\link{download.data}}, \code{\link{Plot.surv.year.annotate}},
 \code{\link{Plot.surv.int.multiyears}},\code{\link{aapc}}.
}
\examples{
data("breast.example", package="JPSurv")
data("fit1", package="JPSurv")

yearvar<-"Year_of_diagnosis_1975"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
interval<-"Interval"
nJP<-3
data.graph<-download.data(breast.example,fit1,nJP,yearvar,"graph",
                          subsetStr1,interval,int.select=c(1,3,5))
out.anno<-Plot.dying.year.annotate(data.graph,fit1,nJP,yearvar,obsintvar,
                     predintvar,interval,annotation=1,topanno=1,trend=1)
out<-Plot.dying.year.annotate(data.graph,fit1,nJP,yearvar,obsintvar,
                     predintvar,interval,annotation=0,topanno=0,trend=1)
trend.rch<-out.anno$trends
plot.rch.anno<-out.anno$plot_anno
plot.rch<-out$plot_no_anno
}
