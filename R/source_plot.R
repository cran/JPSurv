jp.plot.cond.year <- function(fit.cond, start.interval=NULL, end.interval=NULL,
                              year.col="Year", interval.col="Interval", 
                              relSurvInt.col="Relative_Survival_Interval", addToYear=0,
                              ylim=NULL) {

  # Check arguments
  check_dataframe(fit.cond, nm="fit.cond")
  if (!is.null(start.interval)) check_integer(start.interval, "start.interval", min=1)
  if (!is.null(end.interval)) check_integer(end.interval, "end.interval", min=2)
  if (!is.null(addToYear)) check_integer(addToYear, "addToYear", min=0)
  check_dataVar(fit.cond, year.col, "year.col", num=0, allow.miss=1)
  check_dataVar(fit.cond, interval.col, "interval.col", num=0, allow.miss=1)
  check_dataVar(fit.cond, relSurvInt.col, "relSurvInt.col", num=0, allow.miss=1)
  check_numeric(addToYear, "addToYear")
  if (length(ylim)) check_numvec(ylim, "ylim", len=2)

  if (is.null(start.interval)) {
    start.interval <- min(fit.cond[, "Start.interval", drop=TRUE], na.rm=TRUE)
  }  
  tmp      <- fit.cond[, "Start.interval", drop=TRUE] %in% start.interval
  fit.cond <- fit.cond[tmp, , drop=FALSE]
  if (!nrow(fit.cond)) stop("ERROR: start.interval is not valid")
  if (is.null(end.interval)) {
    end.interval <- max(fit.cond[, interval.col, drop=TRUE], na.rm=TRUE)
  }     
  tmp             <- fit.cond[, interval.col, drop=TRUE] <= end.interval
  tmp[is.na(tmp)] <- FALSE
  fit.cond        <- fit.cond[tmp, , drop=FALSE]
  if (!nrow(fit.cond)) stop("ERROR: end.interval is not valid")

  intervals         <- as.numeric(fit.cond[, interval.col, drop=TRUE])
  years             <- addToYear + as.numeric(fit.cond[, year.col, drop=TRUE])
  yvec              <- as.numeric(fit.cond[, "pred_cum", drop=TRUE])
  tmp0              <- intervals == end.interval
  tmp0[is.na(tmp0)] <- FALSE
  main              <- paste0("Conditional Survival at ", end.interval, 
                              " Years Given Survival to ", start.interval, " Year")
  RelSurv           <- as.numeric(fit.cond[, relSurvInt.col, drop=TRUE])

  uyears            <- unique(years)
  obs_cond          <- rep(NA, length(uyears))
  names(obs_cond)   <- as.character(uyears)
  tmp               <- intervals %in% end.interval
  nms               <- as.character(years[tmp])
  obs_cond[nms]     <- RelSurv[tmp]
  nms0              <- nms
  a                 <- start.interval + 1
  b                 <- end.interval - 1
  if (a <= b) {
    for (int in b:a) {
      tmp           <- intervals %in% int
      nms           <- as.character(years[tmp])
      if (length(nms)) obs_cond[nms] <- obs_cond[nms]*RelSurv[tmp]
      tmp <- !(nms0 %in% nms)
      if (any(tmp)) obs_cond[nms0[tmp]] <- NA
    }
  }
  if (is.null(ylim)) {
    tmp   <- c(yvec[tmp0], obs_cond)
    ylim  <- c(min(tmp, na.rm=TRUE), max(tmp, na.rm=TRUE)) 
  }

  plot(years[tmp0], yvec[tmp0], lty=1, lwd=2, type="l", ylim=ylim, 
       xlab="Year at diagnosis", ylab="Relative Survival", 
       main=main, col="black")

  tmp <- is.finite(obs_cond)
  if (any(tmp)) points(uyears[tmp], obs_cond[tmp], pch=20, col="blue")


  NULL
}

jp.plot.death <- function(fit.uncond, fit.cond, start.interval=NULL, end.interval=NULL,
                              year.col="Year", interval.col="Interval", 
                              relSurvInt.col="Relative_Survival_Interval", addToYear=0,
                              ylim=NULL) {

  # Check arguments
  check_fit.uncond(fit.uncond)
  check_dataframe(fit.cond, nm="fit.cond")
  if (!is.null(start.interval)) check_integer(start.interval, "start.interval", min=1)
  if (!is.null(end.interval)) check_integer(end.interval, "end.interval", min=2)
  if (!is.null(addToYear)) check_integer(addToYear, "addToYear", min=0)
  check_dataVar(fit.cond, year.col, "year.col", num=0, allow.miss=1)
  check_dataVar(fit.cond, interval.col, "interval.col", num=0, allow.miss=1)
  check_dataVar(fit.cond, relSurvInt.col, "relSurvInt.col", num=0, allow.miss=1)
  check_numeric(addToYear, "addToYear")
  if (length(ylim)) check_numvec(ylim, "ylim", len=2)

  if (is.null(start.interval)) {
    start.interval <- min(fit.cond[, "Start.interval", drop=TRUE], na.rm=TRUE)
  }  
  tmp      <- fit.cond[, "Start.interval", drop=TRUE] %in% start.interval
  fit.cond <- fit.cond[tmp, , drop=FALSE]
  if (!nrow(fit.cond)) stop("ERROR: start.interval is not valid")
  if (is.null(end.interval)) {
    end.interval <- max(fit.cond[, interval.col, drop=TRUE], na.rm=TRUE)
  }     
  tmp      <- fit.cond[, interval.col, drop=TRUE] %in% end.interval
  fit.cond <- fit.cond[tmp, , drop=FALSE]
  if (!nrow(fit.cond)) stop("ERROR: end.interval is not valid")

  obs_prob  <- 1 - as.numeric(fit.cond[, relSurvInt.col, drop=TRUE])
  obs_years <- addToYear + as.numeric(fit.cond[, year.col, drop=TRUE])
  fullpred  <- fit.uncond$fullpredicted
  tmp       <- fullpred[, interval.col, drop=TRUE] %in% end.interval
  years     <- addToYear + as.numeric(fullpred[tmp, year.col, drop=TRUE])
  yvec      <- 1 - as.numeric(fullpred[tmp, "pred_int", drop=TRUE])
  main      <- "Annual Probability of Dying of Cancer"
           
  if (is.null(ylim)) {
    tmp   <- c(yvec, obs_prob)
    ylim  <- c(min(tmp, na.rm=TRUE), max(tmp, na.rm=TRUE)) 
  }

  plot(years, yvec, lty=1, lwd=2, type="l", ylim=ylim, 
       xlab="Year at diagnosis", ylab="Annual Probability of Cancer Death", 
       main=main, col="black")

  tmp <- is.finite(obs_prob)
  if (any(tmp)) points(obs_years[tmp], obs_prob[tmp], pch=20, col="blue")

  NULL
}

jp.plot.surv <- function(fit.uncond, fit.cond, start.interval=NULL, end.interval=NULL,
                              year.col="Year", interval.col="Interval", 
                              relSurvInt.col="Relative_Survival_Interval", addToYear=0,
                              ylim=NULL, yearsToPlot=NULL, legend.pos="bottom") {

  # Check arguments
  check_fit.uncond(fit.uncond)
  check_dataframe(fit.cond, nm="fit.cond")
  if (!is.null(start.interval)) check_integer(start.interval, "start.interval", min=1)
  if (!is.null(end.interval)) check_integer(end.interval, "end.interval", min=2)
  if (!is.null(addToYear)) check_integer(addToYear, "addToYear", min=0)
  check_dataVar(fit.cond, year.col, "year.col", num=0, allow.miss=1)
  check_dataVar(fit.cond, interval.col, "interval.col", num=0, allow.miss=1)
  check_dataVar(fit.cond, relSurvInt.col, "relSurvInt.col", num=0, allow.miss=1)
  check_numeric(addToYear, "addToYear")
  if (length(ylim)) check_numvec(ylim, "ylim", len=2)
  if (length(yearsToPlot)) check_intvec(yearsToPlot, "yearsToPlot")
  check_legend.pos(legend.pos)

  if (is.null(start.interval)) {
    start.interval <- min(fit.cond[, "Start.interval", drop=TRUE], na.rm=TRUE)
  }  
  tmp      <- fit.cond[, "Start.interval", drop=TRUE] %in% start.interval
  fit.cond <- fit.cond[tmp, , drop=FALSE]
  if (!nrow(fit.cond)) stop("ERROR: start.interval is not valid")
  if (is.null(end.interval)) {
    end.interval <- max(fit.cond[, interval.col, drop=TRUE], na.rm=TRUE)
  }     
  tmp             <- fit.cond[, interval.col, drop=TRUE] <= end.interval
  tmp[is.na(tmp)] <- FALSE
  fit.cond <- fit.cond[tmp, , drop=FALSE]
  if (!nrow(fit.cond)) stop("ERROR: end.interval is not valid")

  ints     <- start.interval:end.interval
  nints    <- length(ints)  
  fullpred <- fit.uncond$fullpredicted
  tmp      <- as.numeric(fullpred[, interval.col, drop=TRUE]) %in% ints
  fullpred <- fullpred[tmp, , drop=FALSE]
  if (!nrow(fullpred)) stop("ERROR with start.interval and/or end.interval")

  # Sort by interval and year
  tmp      <- order(as.numeric(fullpred[, interval.col, drop=TRUE]))
  fullpred <- fullpred[tmp, , drop=FALSE]
  tmp      <- order(as.numeric(fullpred[, year.col, drop=TRUE]))
  fullpred <- fullpred[tmp, , drop=FALSE]

  # Get all years for start.interval
  intervals         <- as.numeric(fullpred[, interval.col, drop=TRUE])
  years             <- as.numeric(fullpred[, year.col, drop=TRUE])
  pred_cum          <- as.numeric(fullpred[, "pred_cum", drop=TRUE])
  tmp               <- intervals %in% start.interval
  years1            <- years[tmp]
  nyears            <- length(years1)
  pred_cum1         <- pred_cum[tmp]
  surv_int_mat      <- matrix(data=NA, nrow=nyears, ncol=nints)
  surv_int_mat[, 1] <- 1
  for (j in 2:nints) {
    tmp0 <- intervals %in% ints[j]
    yrs  <- years[tmp0]
    tmp  <- all.equal(years1, yrs)
    if (!is.logical(tmp) || !tmp) stop("ERROR with fit.uncond$fullpredicted years and/or intervals")
    surv_int_mat[, j] <- pred_cum[tmp0]/pred_cum1
  }

  # Get the rows to plot
  if (!length(yearsToPlot)) {
    if (nyears < 5) {
      yearsToPlot <- years1
    } else {
      yearsToPlot <- round(quantile(years1))
    }
  }

  # Check yearsToPlot
  tmp <- (yearsToPlot %in% years1) | (yearsToPlot %in% (addToYear + years1))
  yearsToPlot <- yearsToPlot[tmp]
  if (!length(yearsToPlot)) stop("ERROR with yearsToPlot")

  tmp <- (years1 %in% yearsToPlot) | ((addToYear + years1) %in% yearsToPlot)
  if (!any(tmp)) stop("ERROR with yearsToPlot")
  rows   <- (1:nyears)[tmp]
  nrows  <- length(rows)
  colors <- rainbow(nrows)

  # Sort by interval and year
  tmp       <- order(as.numeric(fit.cond[, interval.col, drop=TRUE]))
  fit.cond  <- fit.cond[tmp, , drop=FALSE]
  tmp       <- order(as.numeric(fit.cond[, year.col, drop=TRUE]))
  fit.cond  <- fit.cond[tmp, , drop=FALSE]
  intervals <- as.numeric(fit.cond[, interval.col, drop=TRUE])
  years     <- as.numeric(fit.cond[, year.col, drop=TRUE])
  relsurv   <- as.numeric(fit.cond[, relSurvInt.col, drop=TRUE])
  obs_mat      <- matrix(data=NA, nrow=nyears, ncol=nints)
  obs_mat[, 1] <- 1
  if (nints > 2) {
    for (j in 2:nints) {
      tmp0 <- intervals %in% ints[j]
      yrs  <- years[tmp0]
      tmp  <- all.equal(years1, yrs)
      if (!is.logical(tmp) || !tmp) stop("ERROR with fit.cond years and/or intervals")
      obs_mat[, j] <- obs_mat[, j-1]*relsurv[tmp0]
    }
  }
  if (!length(ylim)) {
    tmp  <- c(as.vector(surv_int_mat[rows, ]), as.vector(obs_mat[rows, ]))
    ylim <- c(min(tmp, na.rm=TRUE), max(tmp, na.rm=TRUE)) 
  }
  
  main <- paste0("Conditional Survival Given Survival to ", start.interval, " Year")
  plot(ints, surv_int_mat[rows[1], ], lty=1, lwd=2, type="l", ylim=ylim, xlab="Interval", 
       ylab="Relative Survival", main=main, col=colors[1], xaxt="n")
  axis(1, at=seq(ints[1], ints[nints], 1))
  if (nrows > 1) {
    for (j in 2:nrows) {
      lines(ints, surv_int_mat[rows[j], ], lty=1, lwd=2, col=colors[j])
    }
  }
  abline(v=start.interval,col="gray",lty=2,lwd=3)
  for (j in 1:nrows) points(ints, obs_mat[rows[j], ], pch=20, col=colors[j])
  if (all(yearsToPlot < 1900)) yearsToPlot <- yearsToPlot + addToYear
  legend(legend.pos, legend=as.character(yearsToPlot), pch=20, lty=1, lwd=2, col=colors)
  #legend(x=start.interval+0.1, y=ylim[1], legend=as.character(yearsToPlot), pch=20, lty=1, lwd=2, col=colors)

  NULL
}

