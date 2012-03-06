\name{read.seerstat}
\alias{read.seerstat}
\title{Reading a SEER*Stat dataset}
\description{Read the SEER*Stat source files}
\usage{read.seerstat(fileName, keep.missing=TRUE)}
\arguments{
\item{fileName}{The input SEER*Stat file name, no missing values allowed.}
\item{keep.missing}{Whether to delete the missing values.  }
}
\value{The SEER*Stat matrix.}
\examples{
#X = read.seerstat("allstages.bysiteage.sex.year.yearly");
}
