\name{joinpoint.conditional}
\alias{joinpoint.conditional}
\title{Fitting a join point conditional relative survival model from the unconditional model}
\description{Fitting a joinpoint conditional relative survival model from the unconditional model}
\usage{joinpoint.conditional(fit.uncond, start.intervals, end.intervals, njp=NULL)
}
\arguments{
\item{fit.uncond}{Object returned from \code{\link{joinpoint}}.}
\item{start.intervals}{Vector of integers giving the intervals to condition on. }
\item{end.intervals}{Vector of integers giving the end number of intervals.
                     This vector must have the same length and order as \code{start.intervals}
                     with \code{end.intervals[i] > start.intervals[i]}. }
\item{njp}{NULL or the number of joinpoints corresponding to one of the fitted models
           in \code{fit.uncond}. If \code{njp = k}, then the conditional probabilities
           will be based on the model with \code{k} joinpoints.
           If NULL, then the model corresponding to the best fit returned by
           \code{\link{joinpoint}} will be used.}
}

\details{
  This function computes the conditional survival
  \deqn{ P(T > t_{j+k} | T > t_{j}) = \frac{P(T > t_{j+k})}{P(T > t_{j})}, \hspace{0.1in} k = 1, \ldots , m }{
        P(T > t_{j+k} | T > t_{j}) = P(T > t_{j+k})/P(T > t_{j}),  k = 1, ..., m}
}

\value{
  A data frame similar to the \code{fullpredicted} data frame returned from 
  \code{\link{joinpoint}} except that it will only the contain rows corresponding
  to the \code{start.intervals} and \code{end.intervals} that were specified.
  The data frame will also contain the additional column "Start.interval",
  and will be grouped by the \code{start.intervals}.
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

}

\seealso{ \code{\link{joinpoint}}, \code{\link{joinpoint.seerdata}} }

\references{
Yu, B., Huang, L., Tiwari, R. C., Feuer, E. J. and Johnson, K. A. (2009), 
Modeling population-based cancer survival trends by using join point models for grouped survival data.  
Journal of the Royal Statistical Society: Series A, 172, 405-425.
}
