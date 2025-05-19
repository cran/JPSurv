\name{joinpoint.cond}
\alias{joinpoint.cond}
\title{Fitting a join point conditional relative survival model}
\description{Fitting a joinpoint conditional relative survival model}
\usage{joinpoint.cond(data, subset, start.interval, end.interval=NULL,
	year="Year", interval="Interval",
	number.event="Died", number.alive="Alive_at_Start", 
        number.loss="Lost_to_Followup",
	expected.rate="Expected_Survival_Interval",
	model.form=NULL, maxnum.jp=0, proj.year.num=5,
	op=list(), delLastIntvl=FALSE, add.data.cols="_ALL_")
}
\arguments{
\item{data}{Data frame containing all variables in the model.}
\item{subset}{A logical vector of length \code{nrow(data)}, a character string or 
              NULL to include a particular subset
              of \code{data} in the analysis. 
              See \bold{Details} and \bold{Examples}.}
\item{start.interval}{A positive integer giving the number of intervals to condition on. }
\item{end.interval}{A positive integer > \code{start.interval} giving the end number of intervals. }
\item{year}{Column name of \code{data} giving the year or year code. This
            column must be numeric.
            The default is 'Year'.}
\item{interval}{Column name of \code{data} giving the time interval elapsed 
                from the starting time to the event time. This column must be numeric.
                The default is 'Interval'.}
\item{number.event}{Column name of \code{data} giving the number of events or deaths. 
                    This column must be numeric.
                    The default is 'Died'.}
\item{number.alive}{Column name of \code{data} giving the number of subjects alive. 
                    This column must be numeric.
                    The default is 'Alive_at_Start'.}
\item{number.loss}{Column name of \code{data} giving the number of subjects lost to followup. 
                   This column must be numeric.
                   The default is 'Lost_to_Followup'.}
\item{expected.rate}{Column name of \code{data} giving the interval survival. 
                     This column must be numeric.
                     The default is 'Expected_Survival_Interval'.}
\item{model.form}{an object of class "formula": a symbolic description of covariates. Example: ~-1+age+as.factor(stage)}
\item{maxnum.jp}{The maximum number of join points allowed.  The default is zero, which is equivalent to a proportional hazard relative survival model.}
\item{proj.year.num}{The number of projection years for use in the prediction step. Default value is 5 years, with a valid range of 0 to 30 years.}
\item{op}{List of more options. Details --- 
  \itemize{
    \item{\bold{numbetwn}:} integer value, number of skipped obs between joinpoints exclusive
                (not count for the joinpoints). Default is 2.
    \item{\bold{numfromstart}:} integer value, number of skipped obs from the first obs to 
                   joinpoints exclusive (not count for the joinpoint). Default is 3.
    \item{\bold{numtoend}:} integer value, number of skipped obs from the first obs to joinpoints 
                      exclusive (not count for the joinpoint). Default is 4.
  }
}

\item{delLastIntvl}{an logical value indicating whether or not deleting records of last intervals of all years. 
                    The default is FALSE.}
\item{add.data.cols}{Character vector of column names in \code{data} to add onto the
                     returned data frames of results. Use "_ALL_" to add all columns and
                     use NULL to not add any columns.
                     The default is "_ALL_".}
}

\details{
  The data to be included in the analysis must contain unique 
   \code{year}, \code{interval} pairs. If not, then an error will be thrown.
  The \code{subset} option can be used to ensure that there are unique pairs
  of \code{year} and \code{interval}. 

  This function will set up the data based on the value of \code{start} by
  removing all rows with \code{interval} less than or equal to \code{start}
  and then call the \code{\link{joinpoint}} function.
}

\value{
An object of class "joinpoint" will be returned with attributes:
\item{coefficients}{a named vector of coefficients and standard errors}
\item{jp}{the estimates of the join points}
\item{converged}{convergence status}
\item{predicted}{the fitted relative survival rates}
\item{fullpredicted}{the full output matrix, with all year/interval combinations and projections}
\item{xbeta}{the linear predictor}
\item{ll}{log likelihood}
\item{aic}{AIC}
\item{bic}{BIC}
\item{FitList}{a list that contains fitting results for the number of joinpoints = 0,1,...,numJPoints respectively.}
}
\examples{
#Load the provided SEER 18 breast cancer example data.
data("breast.example", package="JPSurv")
 
# Subset of observations to use
subset <- "Age_groups == '00-49' & Breast_stage == 'Localized'"

# Fit the conditional survival join point model with starting
#   interval 5 
fit <- joinpoint.cond(breast.example, subset, 5, 
                 year="Year_of_diagnosis_1975",
                 model.form=NULL, maxnum.jp=0)
}

\seealso{ \code{\link{joinpoint}}, \code{\link{joinpoint.conditional}} }

\references{
Yu, B., Huang, L., Tiwari, R. C., Feuer, E. J. and Johnson, K. A. (2009), 
Modeling population-based cancer survival trends by using join point models for grouped survival data.  
Journal of the Royal Statistical Society: Series A, 172, 405-425.
}
