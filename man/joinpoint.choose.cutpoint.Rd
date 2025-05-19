\name{joinpoint.choose.cutpoint}
\alias{joinpoint.choose.cutpoint}
\title{Choose cutpoint from joinpoint.relaxProp results}
\description{Choose the cutpoint from the joinpoint.relaxProp returned object}
\usage{joinpoint.choose.cutpoint(obj, cutpoint)
}
\arguments{
\item{obj}{Returned object from \code{\link{joinpoint.relaxProp}}.}
\item{cutpoint}{The cutpoint to choose. It should correspond to a
                set of results in \code{obj$all.results}. If not,
                then an error will be thrown. }
}

%\details{}

\value{
A list of class with the following objects:
\item{predicted}{The fitted relative survival rates}
\item{fullpredicted}{The full output matrix, with all year/interval combinations and projections}
\item{fit.uncond}{The fitted model from the (unconditional) joinpoint model
                  corresponding to the cutpoint.}
\item{fit.cond}{The fitted model from the conditional joinpoint model
                  corresponding to the cutpoint.}
}
%\examples{}

\seealso{ \code{\link{joinpoint.relaxProp}} }

\references{
Yu, B., Huang, L., Tiwari, R. C., Feuer, E. J. and Johnson, K. A. (2009), 
Modeling population-based cancer survival trends by using join point models for grouped survival data.  
Journal of the Royal Statistical Society: Series A, 172, 405-425.
}
