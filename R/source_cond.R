joinpoint.conditional <- function(fit.uncond, start.intervals, end.intervals, njp=NULL) {

  check_fit.uncond(fit.uncond) 
  check_intervals(start.intervals, end.intervals)
  check_njp(njp, fit.uncond)

  ret <- jp.cond.main(fit.uncond, start.intervals, end.intervals, njp)

  ret
}

jp.cond.main <- function(fit, start.ints, end.ints, njp) {

  # First save vectors
  interval     <- fit$interval
  year         <- fit$year
  number.alive <- fit$number.alive
  number.event <- fit$number.event
  number.loss  <- fit$number.loss

  # Get the correct fit
  if (length(njp)) fit <- fit$FitList[[njp+1]]
  
  # Put vectors in fitted object
  fit$interval     <- interval
  fit$year         <- year
  fit$number.alive <- number.alive
  fit$number.event <- number.event
  fit$number.loss  <- number.loss


  int.col    <- fit$interval
  year.col   <- fit$year
  df         <- fit$fullpredicted
  intervals  <- getDataVec(df, int.col)
  uints      <- unique(intervals)
  tmp        <- (start.ints %in% uints) & (end.ints %in% uints)
  start.ints <- start.ints[tmp]
  end.ints   <- end.ints[tmp]
  if (!length(start.ints)) stop("ERROR: intervals are not valid")

  # Get unique intervals
  tmp        <- unique(cbind(start.ints, end.ints))
  start.ints <- tmp[, 1, drop=TRUE]
  end.ints   <- tmp[, 2, drop=TRUE]

  # Data should already be ordered
  df <- orderDataByIntYear(df, year.col, int.col)

  # Set missing alive, dead, loss to previous
  df <- jp.cond.setMiss(df, fit)

  # New column in output for starting interval
  fit$interval.start <- "Start.interval"

  # Compute conditional estimates and standard errors for each start/end in all years
  intervals <- getDataVec(df, int.col)
  years     <- getDataVec(df, year.col)
  uyears    <- unique(years)
  nyears    <- length(uyears)
  ret       <- NULL
  for (i in 1:nyears) {
    tmp <- years %in% uyears[i]
    x   <- jp.cond.getEstSE.1yr(df[tmp, , drop=FALSE], start.ints, end.ints, fit)
    if (!length(ret)) {
      ret <- x
    } else {
      ret <- rbind(ret, x)
    }
  }

  ret <- jp.cond.setReturn(ret, start.ints, end.ints, fit)
  ret
}

jp.cond.setReturn <- function(ret, startvec, endvec, obj) {

  n          <- nrow(ret)
  st.int.col <- obj$interval.start
  avec       <- ret[, st.int.col, drop=TRUE]
  runvec     <- rep(Inf, n)
  ustart     <- sort(unique(startvec))
  for (i in 1:length(ustart)) {
    tmp <- avec %in% ustart[i]
    if (any(tmp)) runvec[tmp] <- i
  }
  ord <- order(runvec)
  ret <- ret[ord, , drop=FALSE]
  ret
}

jp.cond.getEstSE.1yr <- function(x, startvec, endvec, obj) {

  ustartvec   <- unique(startvec)
  nstart      <- length(ustartvec)
  int.col     <- obj$interval
  st.int.col  <- obj$interval.start
  intervals   <- getDataVec(x, int.col)
  rownames(x) <- as.character(intervals)
  x0          <- x
  ret         <- NULL

  # Compute conditional probs for all possible columns in x
  #vv  <- c("pred_int", "pred_cum", 
  #         "Expected_Survival_Interval", "Expected_Survival_Cum",
  #         "Observed_Survival_Interval", "Observed_Survival_Cum",   
  #         "Relative_Survival_Interval", "Relative_Survival_Cum")
  #tmp <- vv %in% colnames(x)
  #vv  <- vv[tmp]

  estv1       <- "pred_int"
  estv2       <- "pred_cum"
  #estv3       <- "Expected_Survival_Interval"

  for (i in 1:nstart) {
    x    <- x0
    t0   <- ustartvec[i]
    tmp  <- startvec %in% t0  
    m    <- sum(tmp)
    if (!m) next
    end2 <- endvec[tmp]
    tmp  <- end2 %in% intervals
    end2 <- end2[tmp]
    m    <- length(end2)
    if (!m) next
    t1   <- max(end2)

    # Compute estimates
    tmp <- (intervals >= t0) & (intervals <= t1)
    if (!any(tmp)) stop("INTERNAL CODING ERROR 1")
    x          <- x[tmp, , drop=FALSE]
    ct0        <- as.character(t0)
    beta_int0  <- x[ct0, estv1, drop=TRUE]
    beta_cum0  <- x[ct0, estv2, drop=TRUE]
    #beta_exp0  <- x[ct0, estv3, drop=TRUE]
    x[, estv1] <- jp.checkProbVec(x[, estv1]/beta_int0)
    x[, estv2] <- jp.checkProbVec(x[, estv2]/beta_cum0)
    #x[, estv3] <- jp.checkProbVec(x[, estv3]/beta_exp0)

    # Compute standard errors
    x <- cond_se(x, obj)

    # Subset by the intervals we want
    tmp <- getDataVec(x, int.col) > t0
    tmp[is.na(tmp)] <- FALSE
    x   <- x[tmp, , drop=FALSE]
    if (!nrow(x)) stop("INTERNAL CODING ERROR 2")
    
    # Order by interval, year is constant
    tmp <- order(getDataVec(x, int.col))
    x   <- x[tmp, , drop=FALSE]

    # Add starting interval
    x[, st.int.col] <- t0
    if (!length(ret)) {
      ret <- x
    } else {
      ret <- rbind(ret, x)
    }
  }

  ret
}

jp.checkProbVec <- function(vec) {

  tmp <- vec > 1
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) vec[tmp] <- 1
  vec
}

jp.cond.setMiss <- function(df, obj) {

  # df must be sorted by year and interval

  cola      <- obj$number.alive
  cold      <- obj$number.event
  coll      <- obj$number.loss
  alive     <- getDataVec(df, cola)
  died      <- getDataVec(df, cold)
  lost      <- getDataVec(df, coll)
  intervals <- getDataVec(df, obj$interval)
  years     <- getDataVec(df, obj$year)

  # See if there is missing data
  tmp1 <- !is.finite(alive) 
  tmp2 <- !is.finite(died) 
  tmp3 <- !is.finite(lost)
  tmp  <- tmp1 | tmp2 | tmp3
  if (!any(tmp)) return(df)

  # Start walking through the data
  uyears    <- unique(years)
  nr        <- nrow(df)
  vec       <- 1:nr
  imp.rows  <- vec[tmp]
  if (imp.rows[1] < 2) return(df) # Need data for at least the first row
  for (row in imp.rows) {
    flag1   <- tmp1[row]
    flag2   <- tmp2[row]
    flag3   <- tmp3[row]
    year    <- years[row]
    int     <- intervals[row]
    row0    <- row - 1
    year0   <- years[row0]
    int0    <- intervals[row0]
    row.imp <- 0
    # First check if the previous row (row0) is the same year as current row
    if (year0 == year) {
      # Use this row for imputed values. Since data is sorted, it will be previous interval.
      row.imp <- row0
    } else {
      # Use row from year0
      tmp0 <- years == year0
      # Get the ideal interval from previous year
      use.int <- year - year0 + int
      tmp     <- tmp0 & (intervals <= use.int)
      if (!any(tmp)) next
      rows    <- vec[tmp]
      row.imp <- max(rows) 
    }
    if (row.imp) {
      if (flag1) alive[row] <- alive[row.imp]
      if (flag2) died[row]  <- died[row.imp]
      if (flag3) lost[row]  <- lost[row.imp]
    }
  }
  df[, cola] <- alive
  df[, cold] <- died
  df[, coll] <- lost

  df
}

jp.set.miss.to.last <- function(vec) {

  # vec is ordered
  tmp <- !is.finite(vec)
  if (!any(tmp) || all(tmp)) return(vec)

  # Check for consecutive indices
  inds <- 1:length(vec)
  ok   <- inds[!tmp]
  if (all(ok == 1:length(ok))) {
    ind      <- max(ok)
    vec[tmp] <- vec[ind]
  } else {
    ivec <- inds[tmp]
    ivec <- ivec[ivec > 1]
    if (!length(ivec)) return(vec) 
    for (i in ivec) vec[i] <- vec[i-1]
  }
  vec

}

joinpoint.cond <- function(data, subset, start.interval, end.interval=NULL, 
                      year="Year", interval="Interval", number.event="Died", 
                      number.alive="Alive_at_Start", number.loss="Lost_to_Followup", 
                      expected.rate="Expected_Survival_Interval", model.form=NULL,
                      maxnum.jp=0, proj.year.num=5, op=list(), delLastIntvl=FALSE,
                      add.data.cols="_ALL_") {

  # Check for errors
  observedrelsurv <- NULL
  check_dataframe(data) 
  data.cols <- colnames(data) 
  check_subset(subset, nrow(data))
  check_integer(start.interval, "start.interval", min=1) 
  check_dataVar(data, year, "year")
  check_dataVar(data, interval, "interval")
  check_dataVar(data, number.event, "number.event")
  check_dataVar(data, number.alive, "number.alive")
  check_dataVar(data, number.loss, "number.loss")
  check_dataVar(data, expected.rate, "expected.rate", allow.miss=1)
  check_formula(model.form, data.cols)
  check_integer(maxnum.jp, "maxnum.jp", valid=0:10) 
  check_integer(proj.year.num, "proj.year.num", valid=0:30)
  check_logical(delLastIntvl, "delLastIntvl") 
  op <- check_op(op)
  end.interval <- check_end.interval(end.interval, start.interval, data, interval)
  add.data.cols <- check_add.data.cols(add.data.cols)

  # Get the correct subset of data
  data <- getDataForJoinpoint(data, NULL, subset, start.interval, year, op) 

  objlist <- list(year=year, interval=interval, number.event=number.event, 
                  number.alive=number.alive, number.loss=number.loss, expected.rate=expected.rate, 
                  observedrelsurv=observedrelsurv, model.form=model.form, maxnum.jp=maxnum.jp,
                  proj.year.num=proj.year.num, start.interval=start.interval, end.interval=end.interval, 
                  op=op, delLastIntvl=delLastIntvl, add.data.cols=add.data.cols)
  checkDataObjects(data, objlist) 

  ret <- joinpoint.cond_main(data, objlist)

  ret
}

joinpoint.cond_main <- function(data, objlist) {

  DEBUG <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: joinpoint.cond_main\n")
 
  data <- applyStartToData(data, objlist)

  # Data is now set up for original joinpoint function
  ret <- joinpoint(data, subset=NULL, year=objlist$year, interval=objlist$interval,
             number.event=objlist$number.event, number.alive=objlist$number.alive, 
             number.loss=objlist$number.loss, expected.rate=objlist$expected.rate, 
             observedrelsurv=objlist[["observedrelsurv", exact=TRUE]],
             model.form=objlist[["model.form", exact=TRUE]], maxnum.jp=objlist$maxnum.jp, 
             proj.year.num=objlist$proj.year.num, op=objlist$op, delLastIntvl=objlist$delLastIntvl, 
             add.data.cols=objlist[["add.data.cols", exact=TRUE]])

  # Update full and predicted matrices for best fit and for each fit in FitList
  ret     <- computeCondSurv.fit(ret, objlist) 
  FitList <- ret[["FitList", exact=TRUE]]
  n       <- length(FitList)
  if (n) {
    for (i in 1:n) {
      fit          <- FitList[[i]]
      FitList[[i]] <- computeCondSurv.fit(fit, objlist) 
    }
    ret$FitList <- FitList
  }

  # Add interval number to ret to be able to shift intervals later. See email
  #  on 2024-08-28 and 2024-10-01 related to error with aapc.multiints.
  # For relax.prop function, start.interval is set to int in objlist
  int                <- objlist$start.interval
  ret$shift.interval <- int
  lst <- ret[["FitList", exact=TRUE]]
  n   <- length(lst)
  if (n && is.list(lst)) {
    for (i in 1:n) lst[[i]]$shift.interval <- int
    ret$FitList <- lst
  } 

  if (DEBUG) cat("End: jointpoint.cond_main\n")

  ret
}

updateIntervalVec <- function(ret, objlist) {

  if (!is.list(ret)) return(ret)
  add <- objlist[["start.interval", exact=TRUE]]
  if (!length(add)) stop("ERROR 1")
  vec <- ret[["Interval", exact=TRUE]]
  if (is.numeric(vec)) ret$Interval <- vec + add
  FitList <- ret[["FitList", exact=TRUE]]
  if (!is.list(FitList)) return(ret)
  n <- length(FitList)
  if (!n) return(ret)
  for (i in 1:n) {
    tmp <- FitList[[i]]
    if (!is.list(tmp)) next
    vec <- tmp[["Interval", exact=TRUE]]
    if (is.numeric(vec)) {
      tmp$Interval <- vec + add
      FitList[[i]] <- tmp  
    }
  }
  ret$FitList <- FitList
  
  ret
}

computeCondSurv.fit <- function(cox_fit, objlist) {

  # Previously, we removed rows of data with interval <= start, so now
  #  interval 1 is from start - <start+1 years
  
  DEBUG    <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: computeCondSurv.fit\n")
  int.var  <- objlist$interval
  year.var <- objlist$year
  start    <- objlist$start.interval
  end      <- objlist$end.interval
  pred     <- cox_fit$predicted
  full     <- cox_fit$fullpredicted  

  # Add start to interval
  pred[, int.var] <- pred[, int.var] + start
  full[, int.var] <- full[, int.var] + start
  start.cond      <- start + 1  

  # Subset data by start to end
  tmp  <- pred[, int.var] %in% (start):(end)
  pred <- pred[tmp, , drop=FALSE]
  if (!nrow(pred)) stop("ERROR 1")
  tmp  <- full[, int.var] %in% (start):(end)
  full <- full[tmp, , drop=FALSE]
  if (!nrow(full)) stop("ERROR 2")

  # Get row ids for predicted and predictedfull
  pred.ids <- paste(pred[, year.var], pred[, int.var], sep=":")
  full.ids <- paste(full[, year.var], full[, int.var], sep=":")
 
  # Subset
  full            <- as.data.frame(full, stringsAsFactors=FALSE, check.names=FALSE)
  ints            <- full[, int.var]
  full[, "label"] <- paste0("P(T>", ints, "|T>", ints-1, ")") 

  tmp                   <- full.ids %in% pred.ids
  pred                  <- full[tmp, , drop=FALSE]
  cox_fit$predicted     <- pred
  cox_fit$fullpredicted <- full

  if (DEBUG) cat("End: computeCondSurv.fit\n")
  cox_fit

}

cond_se <- function(df, objlist) {

  # Data should already be ordered by year and interval
  
  alive    <- df[, objlist$number.alive, drop=TRUE]
  died     <- df[, objlist$number.event, drop=TRUE]
  lost     <- df[, objlist$number.loss, drop=TRUE]
  years    <- df[, objlist$year, drop=TRUE]
  pred_int <- df[, "pred_int", drop=TRUE]
  pred_cum <- df[, "pred_cum", drop=TRUE]
  uyears   <- unique(years)
  nyears   <- length(uyears)
  se.int   <- rep(NA, nrow(df))
  se.cum   <- se.int
  for (i in 1:nyears) {
    tmp         <- years %in% uyears[i]
    d2          <- died[tmp]
    a2          <- alive[tmp]
    l2          <- lost[tmp]
    se.int[tmp] <- cond_se_1year.int(pred_int[tmp], d2, a2, l2)
    se.cum[tmp] <- cond_se_1year(pred_cum[tmp], d2, a2, l2)
  }
  df[, "pred_int_se"] <- se.int
  df[, "pred_cum_se"] <- se.cum
  df
}

cond_se_1year <- function(cond.est, died, alive, lost) {

  effss    <- alive - lost/2
  mhat     <- died/effss
  oneMmhat <- 1 - mhat
  surv     <- cumprod(oneMmhat)
  probDie  <- 1 - surv
  vec      <- probDie/(effss - died)
  csumvec  <- cumsum(vec)
  se       <- cond.est*sqrt(csumvec)

  se
}

cond_se_1year.int <- function(cond.est, died, alive, lost) {

  effss    <- alive - lost/2
  vec      <- died/(effss*(effss - died))
  csumvec  <- cumsum(vec)
  se       <- cond.est*sqrt(csumvec)

  se
}

computeCondSurv.fit0 <- function(cox_fit, objlist) {

  # Want to compute, for example,  P(T>10 | T > 5) = S(10)/S(5)
  # Previously, we removed rows of data with interval <= start, so now
  #  interval 1 is from start - <start+1 years
  
  DEBUG    <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: computeCondSurv.fit\n")
  int.var  <- objlist$interval
  year.var <- objlist$year
  start    <- objlist$start.interval
  end      <- objlist$end.interval
  pred     <- cox_fit$predicted
  full     <- cox_fit$fullpredicted  
  Interval <- cox_fit[["Interval", exact=TRUE]]
  iZ0      <- cox_fit[["iZ0", exact=TRUE]]

  # Add start to interval
  pred[, int.var] <- pred[, int.var] + start
  full[, int.var] <- full[, int.var] + start
  start.cond      <- start + 1  

  # Subset data by start+1 to end
  tmp  <- pred[, int.var] %in% (start+1):(end)
  pred <- pred[tmp, , drop=FALSE]
  if (!nrow(pred)) stop("ERROR 1")
  tmp  <- full[, int.var] %in% (start+1):(end)
  full <- full[tmp, , drop=FALSE]
  if (!nrow(full)) stop("ERROR 2")

  # Get row ids for predicted and predictedfull
  pred.ids <- paste(pred[, year.var], pred[, int.var], sep=":")
  full.ids <- paste(full[, year.var], full[, int.var], sep=":")
 
  # pred and full are ordered by year and interval
  # loop over each row of full and compute covar matrix for each year
  nr          <- nrow(full)
  year0       <- -1
  int0        <- -1
  cond.int    <- rep(NA, nr)
  cond.int.se <- rep(NA, nr)
  cond.cum    <- rep(NA, nr)
  cond.cum.se <- rep(NA, nr)

  for (i in 1:nr) {
    condSurv.int    <- NA
    condSurv.int.se <- NA
    condSurv.cum    <- NA
    condSurv.cum.se <- NA
    year            <- full[i, year.var]
    int             <- full[i, int.var]
    S_predInt       <- full[i, "pred_int"]
    S_predCum       <- full[i, "pred_cum"]
    if (year != year0) {
      # Compute covar matrix for S(t_j)
      tmp     <- calcCovPred(cox_fit, year, Interval, iZ0)
      cov.int <- tmp$cov.int
      cov.cum <- tmp$cov.cum

      # Initialize
      S0_predInt <- NA
      S0_predCum <- NA
    }
    if (int == start.cond) {
      # save survival prob P(T > start)
      S0_predInt <- S_predInt
      S0_predCum <- S_predCum
    } else {
      condSurv.int   <- S_predInt/S0_predInt
      condSurv.cum   <- S_predCum/S0_predCum

      # Get subset of covariance matrix
      cols <- c(start.cond, int) - start
      if (!all(cols %in% 1:ncol(cov.int))) stop("ERROR with column numbers")

      # Compute variance by delta method
      condSurv.int.se <- sqrt(calcCondSurvVariance(S_predInt, S0_predInt, cov.int[cols, cols]))
      condSurv.cum.se <- sqrt(calcCondSurvVariance(S_predCum, S0_predCum, cov.cum[cols, cols]))
    }
    cond.int[i]    <- condSurv.int
    cond.int.se[i] <- condSurv.int.se
    cond.cum[i]    <- condSurv.cum
    cond.cum.se[i] <- condSurv.cum.se
    year0          <- year
    int0           <- int
  }
  # cond surv should be less than 1
  tmp <- cond.int > 1
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) cond.int[tmp] <- 1
  tmp <- cond.cum > 1
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) cond.cum[tmp] <- 1

  full[, "pred_int"]    <- cond.int
  full[, "pred_int_se"] <- cond.int.se
  full[, "pred_cum"]    <- cond.cum
  full[, "pred_cum_se"] <- cond.cum.se
  
  # Subset
  tmp             <- full[, int.var] != start.cond
  full            <- full[tmp, , drop=FALSE]
  full.ids        <- full.ids[tmp]
  full            <- as.data.frame(full, stringsAsFactors=FALSE, check.names=FALSE)
  ints            <- full[, int.var]
  full[, "label"] <- paste0("P(T>", ints-1, "|T>", start, ")") 
  full[, int.var] <- full[, int.var] - 1
  tmp             <- full.ids %in% pred.ids
  pred            <- full[tmp, , drop=FALSE]

  cox_fit$predicted      <- pred
  cox_fit$fullpredicted  <- full

  if (DEBUG) cat("End: computeCondSurv.fit\n")
  cox_fit

}

getOption <- function(op, nm, retIfMiss=NULL) {
  
  if (!is.list(op)) return(retIfMiss)
  if (nm %in% names(op)) {
    ret <- op[[nm, exact=TRUE]]
  } else {
    ret <- retIfMiss 
  }
  ret
}

calcCondSurvVariance <- function(S_jplusk, S_j, cov) {

  # Multivariate delta method to compute variance of S(t_{j+k} | x)/S(t_j | x)
  gradient      <- c(1/S_j, -S_jplusk/(S_j*S_j))
  dim(gradient) <- c(2, 1)
  ret           <- t(gradient) %*% cov %*% gradient
  dim(ret)      <- NULL
  ret

}

calcJacobPred <- function(cox_fit, year, Interval, iZ0) {

  ret <- Predict(cox_fit, year, 0, Interval, cox_fit$jp, years=year, intervals=NULL,
                 beta_input=NULL, Z0=iZ0, gamma_input=NULL, ret.jacob=TRUE)
  if (any(is.na(ret$jacob.int))) stop("ERROR 1 with Jacobian")
  if (any(is.na(ret$jacob.cum))) stop("ERROR 2 with Jacobian")

  ret
}

calcCovPred <- function(cox_fit, year, Interval, iZ0) {

  # Compute covariance matrix for predicted survival probs for each interval
  # Use deleta method: J %*% VCOV * t(J)
  tmp       <- calcJacobPred(cox_fit, year, Interval, iZ0)
  jacob.int <- tmp$jacob.int
  jacob.cum <- tmp$jacob.cum

  cov       <- cox_fit$covariance
  ret.int   <- jacob.int %*% cov %*% t(jacob.int)
  ret.cum   <- jacob.cum %*% cov %*% t(jacob.cum)

  list(cov.int=ret.int, cov.cum=ret.cum)
}

checkDataObjects <- function(data, objlist) {

  DEBUG   <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: checkDataObjects\n")
  start   <- objlist[["start.interval", exact=TRUE]]
  fup.var <- objlist$interval
  vec     <- getDataVec(data, fup.var)
  maxint  <- max(vec, na.rm=TRUE)
  if (length(start) && (start >= maxint)) stop(paste0("ERROR: the value of start.interval cannot exceed ", maxint-1))
  
  # Check for repeated year of diagnosis and follow up year
  yr.var <- objlist$year
  tmp    <- duplicated(data[, c(yr.var, fup.var), drop=FALSE])
  if (any(tmp)) {
    #print(data[tmp, c(yr.var, fup.var), drop=FALSE])
    stop("ERROR: data has repeated year/interval pairs")
  }

  # Check that for each year, intervals are from 1, 2, ....
  ints <- getDataVec(data, fup.var)
  yrs  <- getDataVec(data, yr.var)
  uyrs <- unique(yrs)
  nyrs <- length(uyrs)
  for (i in 1:nyrs) {
    tmp    <- yrs %in% uyrs[i]
    uint   <- unique(ints[tmp]) 
    maxint <- max(uint, na.rm=TRUE)
    if (!all(uint %in% 1:maxint)) {
      stop(paste0("ERROR: incorrect intervals for year = ", uyrs[i]))
    }
  }

  if (DEBUG) cat("End: checkDataObjects\n")

  NULL

}

orderDataByIntYear <- function(data, yr.var, int.var) {

  data <- data[order(getDataVec(data, int.var)), , drop=FALSE]
  data <- data[order(getDataVec(data, yr.var)), , drop=FALSE]
  data

}

applyStartToData <- function(data, objlist) {

  start    <- objlist$start.interval
  fup.var  <- objlist$interval
  yr.var   <- objlist$year
  rate.var <- objlist$expected.rate
  vec      <- getDataVec(data, fup.var)
  maxint   <- max(vec, na.rm=TRUE)
  if (start >= maxint) stop(paste0("ERROR: the value of start.interval cannot exceed ", maxint-1))

  vec             <- vec - start
  data[, fup.var] <- vec 
  tmp             <- vec >= 1
  if (any(is.na(tmp))) stop("ERROR 1")
  data <- data[tmp, , drop=FALSE]
  nr   <- nrow(data)
  if (nr < 1) stop("ERROR: all rows of data removed after applying start.interval")
  
  # Cumulative survival, compute year by year
  # Order by year and interval
  data     <- orderDataByIntYear(data, yr.var, fup.var)
  
  #years   <- getDataVec(data, yr.var)
  #uyears  <- sort(unique(years))
  #nyears  <- length(uyears)
  #cumvec  <- rep(NA, nrow(data))
  #for (i in 1:nyears) {
  #  tmp          <- years == uyears[i]
  #  vec          <- getDataVec(data[tmp, , drop=FALSE], rate.var)
  #  cumvec[tmp]  <- cumprod(vec)
  #}

  data
}

subsetDataObj <- function(data, objlist) {

  nms  <- c("year", "interval", "number.event", "number.alive", "number.loss",
            "expected.rate", "observedrelsurv")
  nr0  <- nrow(data)
  prnt <- objlist$op$print
  keep <- rep(TRUE, nr0)
  vars <- NULL
  for (v in nms) {
    var  <- objlist[[v, exact=TRUE]]
    if (is.null(var)) next
    vec  <- as.numeric(getDataVec(data, var))
    tmp  <- is.finite(vec)  
    keep <- keep & tmp
    if (prnt) {
      m <- sum(!tmp)
      if (m) message(paste0(m, " rows removed due to missing values in ", var))
    }
    data[, var] <- vec
    vars        <- c(vars, var) 
  } 
  if (!all(keep)) data <- data[keep, , drop=FALSE]
  if (nrow(data) < 1) stop("ERROR: no rows left in data after removing missing values")

  form    <- objlist[["model.form", exact=TRUE]]
  allvars <- NULL
  if (!is.null(form)) allvars <- all.vars(form)
  allvars <- unique(c(allvars, vars))
  data    <- data[, allvars, drop=FALSE]
  data

}

cond_modifyData <- function(data, obj.list) {

  start   <- obj.list$start.year.num
  nms     <- c("expected.rate", "observedrelsurv")
  int.var <- obj.list$interval

  vec <- getDataVec(data, int.var)
  tmp <- vec <= start
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) {
    data <- data[tmp, , drop=FALSE]
  }
  data[, int.var] <- vec - start
  
  ord      <- order(data[, int.var])
  data     <- data[ord, , drop=FALSE]
  ord      <- order(data[, obj.list$year])
  data     <- data[ord, , drop=FALSE]
  for (nm in nms) {
    var <- obj.list[[nm, exact=TRUE]]
    if (!length(var)) next
    vec <- getDataVec(data, int.var)
    data[, var] <- cumprod(vec)
  }
  max.year <- obj.list$max.year
  tmp      <- getDataVec(data, obj.list$year) <= max.year
  tmp[is.na(tmp)] <- FALSE
  data     <- data[tmp, , drop=FALSE]

  list(data=data, obj.list=obj.list)
}

getDataVec <- function(data, col) {
  unlist(data[, col, drop=TRUE])
}

isString <- function(x) {
  (length(x) == 1) && is.character(x)
}

isNumber <- function(x) {
  (length(x) == 1) && is.numeric(x)
}

getDicFile <- function(data) {
  ret <- paste0(data, ".dic")
  if (!file.exists(ret)) ret <- paste0(ret, ".gz")
  if (!file.exists(ret)) stop("ERROR")
  ret
}

getDataFile <- function(data) {
  ret <- paste0(data, ".txt")
  if (!file.exists(ret)) ret <- paste0(ret, ".gz")
  if (!file.exists(ret)) stop("ERROR")
  ret
}

getDataCols <- function(data) {
  if (is.data.frame(data)) {
    ret <- colnames(data)
  } else {  
    dicf <- getDicFile(data)
    ret  <- dic.getColnames(dicf) 
  }
  ret
}

getDataNrows <- function(data) {
  if (is.data.frame(data)) {
    ret <- nrow(data)
  } else {  
    datf <- getDataFile(data)
    ret  <- length(scan(datf, what="char", sep="\n", quiet=TRUE)) 
  }
  ret
}

subsetDataForJoinpoint <- function(data, subset, DEBUG) {

  if (is.null(subset)) return(data)

  if (is.character(subset)){
    data <- try(subset(data, eval(parse(text=subset))), silent=TRUE)
    if ("try-error" %in% class(data)) {
      stop(paste0("ERROR applying subset=", subset, " to data"))
    }
    if (nrow(data) < 1) {
      stop(paste0("ERROR: no rows of data remain after applying subset=", subset))
    }
    if (DEBUG) cat(paste0("After applying subset, nrow(data) = ", nrow(data), "\n"))
  } else if (is.logical(subset)) {
    data <- data[subset, , drop=FALSE]
  }
  
  data
}

getDataForJoinpoint <- function(data, dic.file, subset, start, year.var, op) {

  DEBUG <- op$DEBUG
  if (DEBUG) {
    cat("Begin: getDataForJoinpoint\n")
    cat(paste0("nrow(data) = ", nrow(data), "\n"))
  }

  # Apply subset
  data <- subsetDataForJoinpoint(data, subset, DEBUG) 
  if (start <= 0) return(data)

  # Get information about the data
  #tmp      <- cond.getFileInfo(dic.file)  
  #fup.year <- tmp$fup.year
  #nint     <- tmp$nint
  #diag.mat <- tmp$diag.year.mat

  # Check for option to get minimum number of follow up years
  nfup.ints <- op[["nFupInts", exact=TRUE]]
  if (is.null(nfup.ints)) nfup.ints <- 0

  if (nfup.ints) {
    stop("ERROR: this option is not implemented yet")
    fup.year <- NULL
    diag.mat <- NULL
    # Determine the last year for diagnosis
    max.year <- fup.year - nfup.ints

    if (DEBUG) cat(paste0("Cut year of diagnosis at ", max.year, "\n"))

    # first column of diag.mat is code, second column is year
    tmp <- diag.mat[, 2] <= max.year
    if (all(tmp)) return(data) # no need to subset
    max.code <- max(diag.mat[tmp, 1])
    if (DEBUG) cat(paste0("max.year=", max.year, ", max.code=", max.code, "\n"))
    if (!is.finite(max.code)) stop("ERROR in dictionary file: year")

    # Subset the data by max.year (or max.code). Try the code first
    vec <- getDataVec(data, year.var)
    tmp <- vec <= max.code
    if (!any(tmp)) tmp <- vec <= max.year
    if (!all(tmp)) data <- data[tmp, , drop=FALSE]
    if (!nrow(data)) {
      stop(paste0("ERROR: all rows of data have been removed.\nCheck the data and the value for start.interval."))
    } 
  }
  if (DEBUG) cat("End: getDataForJoinpoint\n")

  data
}

#####################################################################
####################### Dictionary file #############################
#####################################################################

cond.getFileInfo <- function(dic.file) {

  # dic.file can also be a list for the relaxed method
  if (!isString(dic.file)) return(dic.file)

  x         <- dic.scanfile(dic.file)
  fup.year  <- dic.getFollowUpYear(x) 
  nint      <- dic.getNumberOfInts(x)
  delim     <- dic.getDelimiter(x)
  diagYears <- dic.getDiagYears(x)

  list(fup.year=fup.year, nint=nint, delim=delim, diag.year.mat=diagYears)

}

cond.modify.dic.file <- function(dic.file, out, start.time=5) {

  x        <- dic.scanfile(dic.file)
  fup.year <- dic.getFollowUpYear(x) 
  nint     <- dic.getNumberOfInts(x)
  delim    <- dic.getDelimiter(x)

  ret      <- x
  ret      <- dic.changeNumOfInt(x, nint, nint - start.time) 
  ret      <- dic.changeYearDiag(ret, fup.year - start.time)
  ret      <- dic.changeIntervals(ret, nint - start.time)
  if (!is.null(out)) write(ret, file=out, ncolumns=1)

  list(fup.year=fup.year, nint=nint, delim=delim)
}

dic.getColnames <- function(dic.file) {

  x   <- dic.scanfile(dic.file)
  x   <- dic.getChunk(x, chunk="[Life Page Variables]") 
  x   <- parseDelimVec(x, "=", ncol=2, numeric=0)
  ret <- x[, 2]
  ret

}

dic.changeIntervals <- function(x, max.int) {

  rng <- dic.getChunkRange(x, chunk="[Format=Interval]") 
  vec <- rng[1]:rng[2]
  mat <- parseDelimVec(x[vec], "=", ncol=2, numeric=0)
  tmp <- as.numeric(mat[, 1]) > max.int
  tmp[is.na(tmp)] <- FALSE
  if (!any(tmp)) return(x)
  rem <- vec[tmp]
  x   <- x[-rem]
  x
}

dic.changeYearDiag <- function(x, max.year) {

  rng <- dic.getChunkRange(x, chunk="[Format=Year of diagnosis") 
  vec <- rng[1]:rng[2]
  mat <- parseDelimVec(x[vec], "=", ncol=2, numeric=1)
  tmp <- mat[, 2] > max.year
  tmp[is.na(tmp)] <- FALSE
  if (!any(tmp)) return(x)
  rem <- vec[tmp]
  x   <- x[-rem]
  x
  
}

dic.changeNumOfInt <- function(x, from, to) {

  str1 <- paste0("NumberOfIntervals=", from)
  str2 <- paste0("NumberOfIntervals=", to)
  ii   <- grep(str1, x, fixed=TRUE)
  len  <- length(ii)
  if (!len) stop("ERROR 1")
  if (len > 1) stop("ERROR 2")
  x[ii] <- str2 
  x

}

dic.getNumberOfInts <- function(x) {

  x2  <- dic.getChunk(x, chunk="[Session Options]") 
  ii  <- dic.getLineNumber(x2, "NumberOfIntervals")
  ret <- dic.getLineValue(x2[ii], sep="=")
  ret <- as.numeric(ret)
  if (!is.finite(ret)) stop("ERROR 1")
  ret

}

dic.getDelimiter <- function(x) {

  x2  <- dic.getChunk(x, chunk="[Export Options]") 
  ii  <- dic.getLineNumber(x2, "Field delimiter")
  ret <- dic.getLineValue(x2[ii], sep="=")
  ret <- tolower(ret)
  if (ret == "tab") {
    ret <- "\t"
  } else if (ret == "comma") {
    ret <- ","
  }

  ret

}

dic.changeLineValue <- function(x, chunk, line, val.from, val.to) {

  nx    <- length(x)
  tmp   <- dic.getChunkRange(x, chunk=chunk)
  i1    <- tmp[1]
  i2    <- tmp[2]
  lines <- dic.getLineNumbers(x, line)
  tmp   <- (lines >= i1) & (lines <= i2)
  lines <- lines[tmp]
  if (length(lines) != 1) stop("ERROR 1") 
  x[lines] <- gsub(val.from, val.to, x[lines], fixed=TRUE)
  x
}

dic.getIntervals <- function(x) {

  x2  <- dic.getChunk(x, chunk="[Format=Interval]") 
  x2  <- parseDelimVec(x2, "=", ncol=2, numeric=0)
  x2  <- x2[, 2]
  x2  <- gsub("yr", "", x2, fixed=TRUE)
  x2  <- gsub("<", "", x2, fixed=TRUE)
  x2  <- removeWhiteSpace(x2)
  tmp <- grepl("-", x2[1], fixed=TRUE)
  if (!tmp) x2[1] <- paste0("0-", x2[1])
  mat <- parseDelimVec(x2, "0", ncol=2, numeric=1)
  colnames(mat) <- c("start", "end")
  mat
}

dic.getMaxInterval <- function(x) {

  mat <- dic.getIntervals(x)
  ret <- max(mat[, 2], na.rm=TRUE)
  ret
}

dic.getDiagYears <- function(x) {

  x2  <- try(dic.getChunk(x, chunk="[Format=Year of"), silent=TRUE)
  if ("try-error" %in% class(x2)) {
    x2 <- dic.getChunk(x, chunk="[Format=Diagnosis Year]")
  } 
  x2  <- parseDelimVec(x2, "=", ncol=2, numeric=1)
  x2

}

dic.getMaxDiagYear <- function(x) {

  mat <- dic.getDiagYears(x)
  ret <- max(mat[, 2], na.rm=TRUE)
  ret
}

dic.getFollowUpYear <- function(x) {

  x2 <- dic.getChunk(x, chunk="[Session Options]") 
  ii <- try(dic.getLineNumber(x2, "LostToFollowupDate"), silent=FALSE)
  if ("try-error" %in% class(ii)) ii <- try(dic.getLineNumber(x2, "StudyCutoffDate"), silent=FALSE)
  if ("try-error" %in% class(ii)) {
    stop("ERROR: dictionary file is missing LostToFollowupDate and StudyCutoffDate")
  }
  str <- dic.getLineValue(x2[ii], sep="=")
  vec <- getVecFromStr(str, delimiter="/")
  len <- length(vec)
  if (len > 1) vec <- vec[len]
  ret <- as.numeric(vec)
  if (!is.finite(ret) || !(ret %in% 1900:2050)) {
    #print(vec)
    stop("ERROR 2")
  }
  ret

}

dic.getLineValue <- function(line, sep="=") {

  vec <- getVecFromStr(line, delimiter=sep)
  if (length(vec) != 2) {
    #print(line)
    stop(paste0("ERROR parsing line")) 
  }
  ret <- vec[2]
  ret
}

dic.getLineNumber <- function(x, line) {

  xl <- tolower(x)
  ll <- tolower(line)
  i  <- grep(ll, xl, fixed=TRUE)
  if (length(i) != 1) stop("ERROR 1")
  if (!is.finite(i)) stop("ERROR 2")
  i
}

dic.getLineNumbers <- function(x, line) {

  xl  <- tolower(x)
  ll  <- tolower(line)
  ret <- grep(ll, xl, fixed=TRUE)
  if (!length(ret)) stop("ERROR 1")
  ret
}


dic.getChunkRange <- function(x, chunk="[Session Options]") {

  # Chunks should begin with "["
  xl <- tolower(x)
  cl <- tolower(chunk)
  nx <- length(x)

  i1 <- dic.getLineNumber(xl, cl) + 1
  i2 <- grep("[", xl, fixed=TRUE)
  m  <- length(i2)
  if (!m) {
    stop("ERROR 3")
  } else if (m > 1) {
    i2  <- i2 - 1
    tmp <- i2 > i1
    i2  <- i2[tmp]
    len <- length(i2)
    if (len >= 1) {
      i2 <- i2[1]
    } else if (!len) {
      i2 <- nx
    } 
  } else {
    stop("ERROR 4")
  }
  if (i2 <= i1) stop("ERROR 5")

  vec <- 1:nx
  tmp <- (vec %in% i1:i2) & (nchar(x) > 0)
  vec <- vec[tmp]
  len <- length(vec)
  if (len < 2) stop("ERROR 6")
  i1  <- vec[1]
  i2  <- vec[len]
  ret <- c(i1, i2)

  ret
}

dic.getChunk <- function(x, chunk="[Session Options]") {

  range <- dic.getChunkRange(x, chunk=chunk)
  ret   <- x[range[1]:range[2]]  
  ret

}

dic.scanfile <- function(f) {

  x <- scan(f, what="character", quiet=TRUE, sep="\n", blank.lines.skip=FALSE)
  x <- removeWhiteSpace(x)
  x
}

#####################################################################
####################### Utility #####################################
#####################################################################
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} 

removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

}

getVecFromStr <- function(string, delimiter="|") {

  # string       String to break apart. No default
  # delimiter    Delimiter used in string. The default is "|". 

  strsplit(string, delimiter, fixed=TRUE)[[1]]

} 

# Function to break up character vector
parseDelimVec0 <- function(vec, sep, ncol, numeric=0) {

  mat <- unlist(strsplit(vec, sep, fixed=TRUE))
  if (length(mat) != length(vec)*ncol) {
    stop("ERROR: check ncol or if some elements of the vector are missing delimiters")
  }
  if (numeric) mat <- as.numeric(mat)
  mat <- matrix(mat, byrow=TRUE, ncol=ncol)
  return(mat)

  mat   

} # END: parseDelimVec0

parseDelimVec <- function(vec, sep, numeric=0, ncol=0) {

  # Determine if sep is white space
  lensep <- nchar(removeWhiteSpace(sep))
  if (ncol) {
    # Try efficient way
    ret <- try(parseDelimVec0(vec, sep, ncol, numeric=numeric), silent=TRUE)
    if (!("try-error" %in% class(ret))) return(ret)
  }

  if (!lensep) {
    add     <- "[?{"
    addFlag <- 1
  } else {
    add     <- " "
    addFlag <- 0
  }
  add1Flag <- 0

  N   <- length(vec)
  if (lensep) vec <- removeWhiteSpace(vec)
  len <- nchar(vec)

  if (lensep) {
    miss  <- !grepl(sep, vec, fixed=TRUE)
  } else {
    if (length(unique(len)) > 1) stop("ERROR: all elements of vec must have the same number of chars when nchar(sep)=0")
    miss  <- vec %in% ""
  }
  nmiss <- sum(miss)
  if (nmiss == N) stop("ERROR: delimiter not found in vector")

  # Is some elements contain data but not delimiter, then throw error
  temp <- miss & (len > 0)
  if (any(temp)) {
    #print(vec[temp][1])
    stop("ERROR: some vector elements are not missing but have no delimiter")
  }

  # Check str to add
  if (addFlag) {
    temp  <- grepl(add, vec, fixed=TRUE)
    if (any(temp)) stop("ERROR: this function cannot be used with your data") 
  }

  if (lensep) {
    temp     <- substr(vec, 1, lensep) == sep
    if (any(temp)) {
      vec[temp] <- paste(add, vec[temp], sep="")
      add1Flag  <- 1
    }
  }

  # Get the number of columns
  str  <- vec[!miss][1]
  temp <- getVecFromStr(str, delimiter=sep)
  ncol <- length(temp)
  l0   <- nchar(str)
  if (substr(str, l0-lensep+1, l0) == sep) ncol <- ncol + 1
  if (ncol < 2) stop("ERROR: in function")

  # Add to the end
  if (lensep) vec <- paste(vec, add, sep="")

  # For the empty ones
  if (nmiss) {
    str <- rep(paste(add, sep, add, sep=""), ncol-1)
    if (length(str) > 1) str <- paste(str, collapse="", sep="")
    vec[miss] <- str
  }

  mat <- unlist(strsplit(vec, sep, fixed=TRUE))
  if (length(mat) != N*ncol) {
    #print(paste("length(mat)=", length(mat), ", N=", N, ", ncol=", ncol, sep=""))
    stop("Check that each non-missing element of vec has the correct number of delimiters")
  }
  mat <- matrix(mat, byrow=TRUE, ncol=ncol)

  if (nmiss) mat[miss, ] <- ""

  if (add1Flag) mat[, 1] <- gsub(add, "", mat[, 1], fixed=TRUE)
  mat[, ncol] <- gsub(add, "", mat[, ncol], fixed=TRUE)
  mat <- removeWhiteSpace(mat)

  if (numeric) mat <- matrix(as.numeric(mat), byrow=FALSE, ncol=ncol)

  mat   

} 
