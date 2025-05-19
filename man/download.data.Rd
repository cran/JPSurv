\name{download.data}
\alias{download.data}
\title{
Combine inputs and outputs
}
\description{
A function that returns a merged data set including the selected cohort input data
and the accompanying data for plot.surv.year/plot.dying.year for the download feature in the JPSurv web app.
}
\usage{
download.data(input, fit, nJP, yearvar, downloadtype, subset = NULL, 
              interval = "Interval", int.select = NULL)
}
\arguments{
  \item{input}{
The input dataset read in by function joinpoint.seerdata.
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
  \item{downloadtype}{
Either "graph" for graph data or "full" for full data.
}
  \item{subset}{
An optional string specifying a subset of observations used in the fitting process.
}
  \item{interval}{
The variable name for year since diagnosis. The default is 'Interval'.
}
  \item{int.select}{
The interval values selected for the plot if the downloadtype="graph". The default is NULL.
}
}
\value{
A data frame containing the selected cohort input data and the below predicted columns.
  \item{Predicted_Survival_Int}{The predicted interval survival.}
  \item{Predicted_ProbDeath_Int}{The predicted probability of Dying of Cancer.}
  \item{Predicted_Survival_Cum}{The predicted cumulative survival.}
  \item{Predicted_Survival_Int_SE}{The standard error of the predicted interval survival.}
  \item{Predicted_Survival_Cum_SE}{The standard error of the predicted probability of Dying of Cancer.}
  \item{Predicted_ProbDeath_Int_SE}{The standard error of the predicted cumulative survival.}
}
\references{
The JPSurv web app https://analysistools.nci.nih.gov/jpsurv/.
}
\author{
Fanni Zhang <zhangf@imsweb.com>
}
\seealso{
\code{\link{Plot.dying.year.annotate}}, \code{\link{Plot.surv.year.annotate}}, \code{\link{Plot.surv.int.multiyears}}.
}
\examples{
data("breast.example", package="JPSurv")
data("fit1", package="JPSurv")

nJP<-2                 
data.graph135<-download.data(breast.example,fit1,nJP,"Year_of_diagnosis_1975","graph",
                             subsetStr1,interval="Interval",int.select=c(1,3,5)) 
data.full<-download.data(breast.example,fit1,nJP,"Year_of_diagnosis_1975","full",subsetStr1)          
}
