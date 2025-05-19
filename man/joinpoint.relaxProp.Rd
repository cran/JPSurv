\name{joinpoint.relaxProp}
\alias{joinpoint.relaxProp}
\title{ Relaxing the proportionality assumption}
\description{Fitting a joinpoint survival model by relaxing the proportionality assumption}
\usage{joinpoint.relaxProp(data, subset, max.cutpoint=5,
       year="Year", interval="Interval", number.event="Died", 
       number.alive="Alive_at_Start", number.loss="Lost_to_Followup",
       expected.rate="Expected_Survival_Interval", 
       observed.rate="Observed_Survival_Interval", 
       model.form=NULL, maxnum.jp=0, 
       proj.year.num=5, op=list(), delLastIntvl=FALSE, add.data.cols=NULL)
}
\arguments{
\item{data}{Data frame containing all variables in the model.}
\item{subset}{A logical vector of length \code{nrow(data)}, a character string or 
              NULL to include a particular subset
              of \code{data} in the analysis. 
              See \bold{Details} and \bold{Examples}.}
\item{max.cutpoint}{A positive integer or NULL giving the number of cutpoints to consider.
                    If NULL, then it will be set to the number of intervals minus one.
                    The default is 5. }
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
\item{expected.rate}{Column name of \code{data} giving the expected interval survival. 
                     This column must be numeric.
                     The default is 'Expected_Survival_Interval'.}
\item{observed.rate}{Column name of \code{data} giving the observed interval survival. 
                     This column must be numeric.
                     The default is 'Observed_Survival_Interval'.}
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
                     The default is NULL.}
}

\details{
  This function finds the optimal clustering of intervals (1, ..., I),
  where I is the number of intervals, 
  such that there are at most two ordered clusters of the form
  (1, ..., j) and (j+1, ..., I). For each ordered cluster,
  a model is fit and the BIC is computed.
  The algorithm is as follows: \cr
  1. Fit the (unconditional) joinpoint survival model on intervals (1, ..., I) 
      and compute the BIC and call it BIC-0. \cr
  2. For each cutpoint j, j = 1, ..., \code{max.cutpoint}, fit the 
      (unconditional) joinpoint survival model on intervals (1, ..., j) and fit the
      conditional joinpoint survival model on intervals (j+1, ..., I). Compute
      the BIC and label it  BIC-j. \cr
  3. The optimal clustering is the one with minimum BIC = min(BIC-0, BIC-1, ...)
  
}

\value{
A list of class "jp.relaxProp" with the following objects:
\item{fit.info}{A data frame containing fitting information from the joinpoint and
               conditional joinpoint models at each step of the algorithm. The data
               frame contains the joinpoints and number of joinpoints.}
\item{predicted}{The fitted relative survival rates}
\item{fullpredicted}{The full output matrix, with all year/interval combinations and projections}
\item{fit.uncond}{The fitted model from the (unconditional) joinpoint model
                  corresponding to the best fit.}
\item{fit.cond}{The fitted model from the conditional joinpoint model
                  corresponding to the best fit.}
\item{all.results}{A list containing all the results at each cutpoint. 
                   Each element of \code{all.results} is a list containing
                   \code{cutpoint}, \code{fit.cond}, and \code{fit.uncond}.}
}
\examples{
#Load the provided SEER 18 breast cancer example data.
data("breast.example", package="JPSurv")
 
# Subset of observations to use
subset <- "Age_groups == '00-49' & Breast_stage == 'Localized'"
fit    <- joinpoint.relaxProp(breast.example, subset, max.cutpoint=2,
                              year="Year_of_diagnosis_1975")
}

\seealso{ \code{\link{joinpoint}}, \code{\link{joinpoint.cond}} }

\references{
Yu, B., Huang, L., Tiwari, R. C., Feuer, E. J. and Johnson, K. A. (2009), 
Modeling population-based cancer survival trends by using join point models for grouped survival data.  
Journal of the Royal Statistical Society: Series A, 172, 405-425.
}
