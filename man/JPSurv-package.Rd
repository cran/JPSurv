\name{JPSurv-package}
\alias{JPSurv-package}
\alias{JPSurv}
\docType{package}
\title{
\packageTitle{JPSurv}
}
\description{
This package is to analyze trends in survival with respect to year of diagnosis using joinpoint model.
}
\details{
Survival data includes two temporal dimensions that are important to consider: the calendar year of diagnosis 
and the time since diagnosis. The JPSurv model is an extension of Cox proportional hazard model and of Hakulinen and Tenkanen
in the case of relative survival, and fits a proportional hazard joinpoint model to survival data on the log hazard scale. 
Joinpoint models consist of linear segments connected through joinpoints. The probability (hazard) of cancer death is specified as 
the product of a baseline hazard (on time since diagnosis) and a multiplicative factor describing the effect of year of diagnosis 
and possibly other covariates. The effect of year of diagnosis is modeled as joined linear segments on the log scale. 
The number and location of joinpoints are estimated from data and represent the times at which trends changed. 
This model implies that the probability of cancer death as a function of time since diagnosis is proportional for individuals 
diagnosed in different calendar years. The software uses discrete-time survival data, i.e. survival data grouped by years 
since diagnosis in the life table format. The package accommodates both relative survival and cause-specific survival.
}

\references{
Yu B, Huang L, Tiwari RC, Feuer EJ, Johnson KA. 
Modelling population-based cancer survival trends by using join point models for grouped survival data. 
Journal of the Royal Statistical Society Series a-Statistics in Society. 2009;172:405-25.

Hakulinen T, Tenkanen L. Regression Analysis of Relative Survival Rates. Applied Statistics. 1987;36(3):309-17.

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
}
