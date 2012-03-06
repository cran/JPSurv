\name{predict.joinpoint}
\alias{predict.joinpoint}
\title{Predict method for join point models}
\description{Predicted values based on join point object}
\usage{\method{predict}{joinpoint}(object, ...)}
\arguments{
\item{object}{Object of class "joinpoint"}
\item{...}{Two more arguments: years, intervals may be used to specify which years and intervals to predict.}
}
\value{A data frame which contains the predicted interval survivals and cumulative survivals}
\examples{
data(prostate);
fit = joinpoint(~Year, data=prostate, numJPoints = 0);
pred1 = predict(fit);
pred2 = predict(fit, years = 1978, intervals = 5);
pred3 = predict(fit, years = 1976:1980, intervals = c(2, 3, 5));
}
