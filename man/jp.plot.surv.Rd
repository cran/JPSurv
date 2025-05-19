\name{jp.plot.surv}
\alias{jp.plot.surv}
\title{Plot of conditional survival versus interval}
\description{Plot of conditional survival versus interval}
\usage{jp.plot.surv(fit.uncond, fit.cond, start.interval=NULL, end.interval=NULL,
                    year.col="Year", interval.col="Interval", 
                    relSurvInt.col="Relative_Survival_Interval", addToYear=0,
                    ylim=NULL, yearsToPlot=NULL, legend.pos="bottom")
}
\arguments{
\item{fit.uncond}{Object returned from \code{\link{joinpoint}}.}
\item{fit.cond}{Object returned from \code{\link{joinpoint.conditional}}.}
\item{start.interval}{NULL or the starting interval. If NULL, then the smallest
                      starting interval will be chosen. }
\item{end.interval}{NULL or the ending interval. If NULL, then the largest
                    interval will be chosen.}
\item{year.col}{The name of the year at diagnosis column in \code{fit.cond}.
                The default is "Year".}
\item{interval.col}{The name of the interval column in \code{fit.cond}.
                The default is "Interval".}
\item{relSurvInt.col}{The name of the relative survival interval column in \code{fit.cond}.
                The default is "Relative_Survival_Interval".}
\item{addToYear}{Integer to add to the year at diagnosis column to give the correct
                 years at diagnosis for displaying. The default is 0.}
\item{ylim}{NULL or the y-axis limits of the plot. If NULL, then the limits
            will be determined from the data.}
\item{yearsToPlot}{NULL or the years at diagnosis to display in the plot. If NULL,
                   then no more than five years will be selected for the plot.}
\item{legend.pos}{Character string to give the legend position in the plot.
                  The default is "bottom".}
}

\details{
  The lines on the plot use the \code{pred_cum} column in \code{fit.cond},
  while the points are from the \code{relSurvInt.col} column.
}

\value{
  NULL
}
\examples{
#Load the provided SEER 18 breast cancer example data.
data("breast.example", package="JPSurv")
 
# Subset of observations to use
subset <- "Age_groups == '00-49' & Breast_stage == 'Localized'"

# Fit the unconditional survival join point model
fit <- joinpoint(breast.example, subset,
                 year="Year_of_diagnosis_1975",
                 model.form=NULL, maxnum.jp=0)

# Compute conditional survival S(10 | 5) = P(T>10 | T>5)
ret <- joinpoint.conditional(fit, 5, 10) 

jp.plot.surv(fit, ret, year.col="Year_of_diagnosis_1975")
}

\seealso{ \code{\link{joinpoint}}, \code{\link{joinpoint.conditional}} }

%\references{}
