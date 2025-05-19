\name{input.valid}
\alias{input.valid}
\title{
Input data validation
}
\description{
A function that validates the selected cohort data and returns an indicator along with a warning message if the selected cohort is not valid.
}
\usage{
input.valid(input, subset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
Input dataset imported by function joinpoint.seerdata.
}
  \item{subset}{
An optional string specifying a subset of observations used in the fitting process.
}
}
\value{
If the cohort selection is available, value 1 will be returned; 
otherwise, value 0 will be returned and a warning message printed.
}

