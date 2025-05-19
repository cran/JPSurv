joinpoint.relaxProp <- function(data, subset, max.cutpoint=5, 
                      year="Year", interval="Interval", number.event="Died", 
                      number.alive="Alive_at_Start", number.loss="Lost_to_Followup", 
                      expected.rate="Expected_Survival_Interval", 
                      observed.rate="Observed_Survival_Interval", model.form=NULL,
                      maxnum.jp=0, proj.year.num=5, op=list(), delLastIntvl=FALSE,
                      add.data.cols=NULL) {

  # Check for errors
  observedrelsurv <- NULL
  check_dataframe(data) 
  data.cols <- colnames(data) 
  check_subset(subset, nrow(data))
  check_dataVar(data, year, "year")
  check_dataVar(data, interval, "interval")
  check_dataVar(data, number.event, "number.event")
  check_dataVar(data, number.alive, "number.alive")
  check_dataVar(data, number.loss, "number.loss")
  check_dataVar(data, expected.rate, "expected.rate", allow.miss=1)
  check_dataVar(data, observed.rate, "observed.rate", allow.miss=1)
  check_formula(model.form, data.cols)
  check_integer(maxnum.jp, "maxnum.jp", valid=0:10) 
  check_integer(proj.year.num, "proj.year.num", valid=0:30)
  check_logical(delLastIntvl, "delLastIntvl") 
  op <- check_op(op)
  check_max.cutpoint(max.cutpoint, data, interval)
  add.data.cols <- check_add.data.cols(add.data.cols)
  add.data.cols <- unique(c(add.data.cols, observed.rate))

  objlist <- list(year=year, interval=interval, number.event=number.event, 
                  number.alive=number.alive, number.loss=number.loss, expected.rate=expected.rate, 
                  observedrelsurv=observedrelsurv, model.form=model.form, maxnum.jp=maxnum.jp,
                  proj.year.num=proj.year.num, max.cutpoint=max.cutpoint, op=op, 
                  delLastIntvl=delLastIntvl, subset=subset, add.data.cols=add.data.cols,
                  observed.rate=observed.rate)
  ret <- jp_relax_main(data, objlist)
  
  ret
}

jp_relax_main <- function(data, objlist) {

  DEBUG <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: jp_relax_main\n")

  # Set up data. getDataForJP_relax must be called before checkDataObjects
  tmp     <- getDataForJP_relax(data, objlist)
  data    <- tmp$data
  objlist <- tmp$objlist
  rm(tmp)
  gc()
  checkDataObjects(data, objlist) 

  # Data frame of information on model fits
  I      <- objlist$end.interval  # maximum interval, not the max cutpoint
  infoDF <- jp_relax_initRes(objlist$max.cutpoint + 1)

  # Fit unconditional model using all data (1 to end.interval) 
  fit0 <- joinpoint(data, subset=NULL, year=objlist$year, interval=objlist$interval,
             number.event=objlist$number.event, number.alive=objlist$number.alive, 
             number.loss=objlist$number.loss, expected.rate=objlist$expected.rate, 
             model.form=objlist[["model.form", exact=TRUE]], maxnum.jp=objlist$maxnum.jp, 
             proj.year.num=objlist$proj.year.num, op=objlist$op, delLastIntvl=objlist$delLastIntvl,
             add.data.cols=objlist[["add.data.cols", exact=TRUE]])
  bic0            <- fit0$bic
  fit0$bic.uncond <- bic0
  fit0$aic.uncond <- fit0$aic
  fit0$fit.uncond <- fit0
  BIC0            <- jp_relax_compute_BIC_from_Fit(fit0, NULL, 0, objlist) 
  fit0$BIC        <- BIC0

  if (DEBUG) cat(paste0("BIC0=", BIC0, "\n"))
  infoDF <- jp_relax_saveRes(infoDF, fit0, 0, I) 
  fit0$bic.uncond <- NULL
  fit0$aic.uncond <- NULL
  fit0$fit.uncond <- NULL

  # Fit unconditional and conditional models over 1,...,j and j+1,...,I respectively
  # Best fit over clustered fits (min BIC) is returned
  tmp    <- jp_relax_minBIC(data, objlist, infoDF) 
  ret    <- tmp$best
  infoDF <- tmp$info
  all    <- tmp$all
  rm(tmp)
  gc()
  if (!is.null(ret)) {
    if (BIC0 < ret$BIC) {
      ret <- list(fit.uncond=fit0, bic.uncond=bic0, bic=bic0, interval=0)
    }
  }
  if (DEBUG) cat(paste0("Minimum BIC occurs at interval ", ret$interval, "\n"))

  ret <- jp_relax_setReturn(ret, infoDF, objlist, all, fit0)

  if (DEBUG) cat("End: jp_relax_main\n")
  ret
}

jp_relax_compute_BIC0 <- function(intRelSurv, intRelSurvHat, nJP, nInt) {

  N <- length(intRelSurv)
  if (N != length(intRelSurvHat)) stop("ERROR 1")
  tmp <- intRelSurv - intRelSurvHat
  rss <- sum(tmp*tmp, na.rm=TRUE)
  p0  <- 2*nJP + 1 + nInt
  ret <- log(rss/N) + p0*log(N)/N

  ret
}

jp_relax_compute_BIC <- function(intRelSurv1, intRelSurvHat1, intRelSurv2, intRelSurvHat2,
                                 nJP1, nJP2, nInt) {

  n1 <- length(intRelSurv1)
  if (n1 != length(intRelSurvHat1)) stop("ERROR 1")
  n2 <- length(intRelSurv2)
  if (n2 != length(intRelSurvHat2)) stop("ERROR 2")

  tmp  <- intRelSurv1 - intRelSurvHat1
  rss1 <- sum(tmp*tmp, na.rm=TRUE)
  tmp  <- intRelSurv2 - intRelSurvHat2
  rss2 <- sum(tmp*tmp, na.rm=TRUE)
  rss  <- rss1 + rss2
  N    <- n1 + n2
  p    <- 1 + 2*(nJP1 + nJP2) + 2 + nInt
  ret  <- log(rss/N) + p*log(N)/N
  ret
}

jp_relax_compute_BIC_from_Fit <- function(fit.uncond, fit.cond, int, objlist) {

  ret <- NA
  if ("try-error" %in% class(fit.uncond)) return(ret)
  if (int && ("try-error" %in% class(fit.cond))) return(ret)
  ivar    <- objlist$interval
  obsv    <- objlist$observed.rate
  hatv    <- "pred_int"
  x       <- fit.uncond$fullpredicted
  njp1    <- length(fit.uncond[["jp", exact=TRUE]])
  irs1    <- x[, obsv, drop=TRUE]
  irs1hat <- x[, hatv, drop=TRUE] 
  tmp     <- is.finite(irs1) & is.finite(irs1hat)
  nInt1   <- length(unique(x[tmp, ivar, drop=TRUE]))
  irs1    <- irs1[tmp]
  irs1hat <- irs1hat[tmp]

  if (int) {
    x       <- fit.cond$fullpredicted
    njp2    <- length(fit.cond[["jp", exact=TRUE]])
    irs2    <- x[, obsv, drop=TRUE]
    irs2hat <- x[, hatv, drop=TRUE] 
    tmp     <- is.finite(irs2) & is.finite(irs2hat)
    nInt2   <- length(unique(x[tmp, ivar, drop=TRUE]))
    irs2    <- irs2[tmp]
    irs2hat <- irs2hat[tmp]
    ret <- jp_relax_compute_BIC(irs1, irs1hat, irs2, irs2hat, njp1, njp2, nInt1+nInt2)
  } else {
    ret <- jp_relax_compute_BIC0(irs1, irs1hat, njp1, nInt1)
  }
  ret
}

jp_relax_initRes <- function(n) {

  nvec <- rep(NA, n)
  cvec <- rep("", n)
  ret  <- data.frame(Cutpoint=nvec, Intervals.1=cvec, Intervals.2=cvec,
                     AIC=nvec, AIC.1=nvec, AIC.2=nvec, 
                     BIC=nvec, 
                     N.JP.1=nvec, N.JP.2=nvec,
                     JP.1=cvec, JP.2=cvec, 
                     "BestFit(Min BIC)"=cvec, Message=cvec,
                     stringsAsFactors=FALSE, check.names=FALSE)
  ret 
}

jp_relax_saveRes <- function(all, res, interval, end.interval) {

  row                     <- interval + 1
  all[row, "Cutpoint"]    <- interval
  if (interval == 0) {
    all[row, "Intervals.1"] <- paste0("1 - ", end.interval)
  } else {
    all[row, "Intervals.1"] <- paste0("1 - ", interval)
    all[row, "Intervals.2"] <- paste0((interval+1), " - ", end.interval)
  }

  if (!length(res)) return(all)
  if ("try-error" %in% class(res)) {
    all[row, "Message"] <- getErrorMsgFromTryError(res)
    return(all)
  }

  x   <- res[["BIC", exact=TRUE]];        if (!is.null(x)) all[row, "BIC"]   <- x
  #x   <- res[["bic.uncond", exact=TRUE]]; if (!is.null(x)) all[row, "BIC.1"] <- x
  #x   <- res[["bic.cond", exact=TRUE]];   if (!is.null(x)) all[row, "BIC.2"] <- x
  x   <- res[["aic", exact=TRUE]];        if (!is.null(x)) all[row, "AIC"]   <- x
  x   <- res[["aic.uncond", exact=TRUE]]; if (!is.null(x)) all[row, "AIC.1"] <- x
  x   <- res[["aic.cond", exact=TRUE]];   if (!is.null(x)) all[row, "AIC.2"] <- x

  fit <- res[["fit.uncond", exact=TRUE]]
  jp  <- fit[["jp", exact=TRUE]]
  njp <- length(jp)
  all[row, "N.JP.1"] <- njp
  if (njp) all[row, "JP.1"] <- paste0(jp, collapse=" ")
  fit <- res[["fit.cond", exact=TRUE]]
  jp  <- fit[["jp", exact=TRUE]]
  njp <- length(jp)
  all[row, "N.JP.2"] <- njp
  if (njp) all[row, "JP.2"] <- paste0(jp, collapse=" ")

  all
}

jp_relax_minBIC <- function(data, objlist, infoDF) {

  DEBUG  <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: jp_relax_minBIC\n")

  maxcut <- objlist$max.cutpoint  
  I      <- objlist$end.interval
  minbic <- Inf
  ret    <- NULL
  all    <- list()
  for (i in 1:maxcut) {
    if (DEBUG) cat(paste0("  interval = ", i, "\n"))
    tmp      <- try(jp_relax_BIC(data, i, objlist), silent=TRUE)
    all[[i]] <- tmp
    #tmp <- jp_relax_BIC(data, i, objlist)
    if (!("try-error" %in% class(tmp))) {
      bic <- tmp$BIC
      if (bic < minbic) {
        # save results
        minbic <- bic
        ret    <- tmp
      } 
    } else {
      #if (DEBUG) print(tmp)
    }
    infoDF <- jp_relax_saveRes(infoDF, tmp, i, I)
  }
  if (DEBUG) cat("End: jp_relax_minBIC\n")

  list(best=ret, info=infoDF, all=all)
}

jp_relax_BIC <- function(data, int, objlist) {

  DEBUG  <- objlist$op$DEBUG
  if (DEBUG) cat("Begin: jp_relax_BIC\n")

  # At this point, intervals in data are from 1, ..., end.interval

  # Fit jointpoint for intervals 1, ..., int
  tmp  <- data[, objlist$interval] %in% 1:int
  fit1 <- joinpoint(data[tmp, , drop=FALSE], subset=NULL, year=objlist$year, interval=objlist$interval,
             number.event=objlist$number.event, number.alive=objlist$number.alive, 
             number.loss=objlist$number.loss, expected.rate=objlist$expected.rate, 
             model.form=objlist[["model.form", exact=TRUE]], maxnum.jp=objlist$maxnum.jp, 
             proj.year.num=objlist$proj.year.num, op=objlist$op, delLastIntvl=objlist$delLastIntvl,
             add.data.cols=objlist[["add.data.cols", exact=TRUE]])
  bic1 <- fit1$bic
  aic1 <- fit1$aic

  # Fit conditional joinpoint to int+1, ..., N=end.interval
  # The data must be set up first
  objlist$start.interval <- int
  data <- getDataForJoinpoint(data, objlist$dic.list, NULL, int,
                              objlist$year, objlist$op) 
  # start is applied in joinpoint.cond_main
  fit2 <- joinpoint.cond_main(data, objlist) # int is the start.interval in objlist
  bic2 <- fit2$bic
  aic2 <- fit2$aic
  bic  <- bic1 + bic2
  aic  <- aic1 + aic2

  BIC  <- jp_relax_compute_BIC_from_Fit(fit1, fit2, int, objlist)

  if (DEBUG) {
    cat(paste0("BIC=", bic, ", BIC1=", bic1, ", BIC2=", bic2, "\n"))
    cat(paste0("AIC=", aic, ", AIC1=", aic1, ", AIC2=", aic2, "\n"))
    cat(paste0("BIC (new) =", BIC, "\n"))
  }

  if (DEBUG) cat("End: jp_relax_BIC\n")
  list(aic=aic, bic=bic, fit.uncond=fit1, bic.uncond=bic1, aic.uncond=aic1,
       fit.cond=fit2, bic.cond=bic2, aic.cond=aic2, interval=int, BIC=BIC)
}

getDataForJP_relax <- function(data, objlist) {

  DEBUG <- objlist$op$DEBUG
  if (DEBUG) {
    cat("Begin: getDataForJP_relax\n")
    cat(paste0("nrow(data) = ", nrow(data), "\n"))
  }

  # Apply subset (only needs to be applied once)
  data <- subsetDataForJoinpoint(data, objlist[["subset", exact=TRUE]], DEBUG) 
  objlist$subset <- NULL

  # Get information about the data from dictionary file
  #objlist$dic.list <- cond.getFileInfo(objlist$dic.file)  
 
  # Remove intervals from data if needed
  vec    <- getDataVec(data, objlist$interval)
  maxint <- max(vec, na.rm=TRUE)
  objlist$end.interval <- maxint

  if (DEBUG) cat("End: getDataForJP_relax \n")

  list(data=data, objlist=objlist)
}

jp_relax_setLabel <- function(x, int.var) {

  if (length(x)) {
    x[, "label"] <- paste0("P(T>", x[, int.var, drop=TRUE], ")")
  }
  x
}

jp_relax_setReturn <- function(ret, infoDF, objlist, all.res, fit0) {

  int.var    <- objlist$interval
  year.var   <- objlist$year
  fit.uncond <- ret[["fit.uncond", exact=TRUE]]
  fit.cond   <- ret[["fit.cond", exact=TRUE]]
  pred       <- NULL
  full       <- NULL

  if (!is.null(fit.uncond)) {
    pred <- jp_relax_setLabel(fit.uncond$predicted, int.var)
    full <- jp_relax_setLabel(fit.uncond$fullpredicted, int.var)
  }
  if (!is.null(fit.cond)) {
    pred <- rbind(pred, fit.cond$predicted)
    full <- rbind(full, fit.cond$fullpredicted)
  }


  if (length(pred)) pred <- orderDataByIntYear(pred, year.var, int.var)
  if (length(full)) full <- orderDataByIntYear(full, year.var, int.var)

  tmp <- infoDF[, "Cutpoint"] == ret$interval
  infoDF[tmp, "BestFit(Min BIC)"] <- "*"

  all      <- list()
  all[[1]] <- list(cutpoint=0, fit.uncond=fit0)
  for (i in 1:length(all.res)) {
    tmp <- all.res[[i]]
    if (!is.list(tmp)) next
    all[[length(all)+1]] <- list(cutpoint=tmp$interval, fit.uncond=tmp$fit.uncond,
                                 fit.cond=tmp$fit.cond)
  }

  ret <- list(fit.info=infoDF, predicted=pred, fullpredicted=full,
          fit.uncond=fit.uncond, fit.cond=fit.cond, all.results=all)
  class(ret) <- "jp.relaxProp"
  ret
}

getErrorMsgFromTryError <- function(obj) {

  ret <- ""
  msg <- attr(obj, "condition")
  tmp <- msg[["message", exact=TRUE]]
  if (isString(tmp)) ret <- tmp
  ret
} 

joinpoint.choose.cutpoint <- function(obj, cutpoint) {

  check_jp.relaxPropObj(obj)
  check_integer(cutpoint, "cutpoint", valid=NULL, min=0, max=length(obj$all.results)-1) 
  jp.choose.cutpoint.main(obj, cutpoint)
}

jp.choose.cutpoint.main <- function(obj, cutpoint) {

  flag <- 0
  all  <- obj[["all.results", exact=TRUE]]
  nall <- length(all)
  for (i in 1:nall) {
    lst <- all[[i]]
    if (lst$cutpoint == cutpoint) {
      flag <- 1
      break
    }
  }
  if (!flag) stop(paste0("ERROR: cutpoint ", cutpoint, " not found in all.results"))

  fit.uncond <- lst[["fit.uncond", exact=TRUE]] # Should not be NULL
  fit.cond   <- lst[["fit.cond", exact=TRUE]]
  int.var    <- fit.uncond$interval
  year.var   <- fit.uncond$year
  pred       <- NULL
  full       <- NULL
  if (!is.null(fit.uncond)) {
    pred <- jp_relax_setLabel(fit.uncond$predicted, int.var)
    full <- jp_relax_setLabel(fit.uncond$fullpredicted, int.var)
  }
  if (!is.null(fit.cond)) {
    pred <- rbind(pred, fit.cond$predicted)
    full <- rbind(full, fit.cond$fullpredicted)
  }
  if (length(pred)) pred <- orderDataByIntYear(pred, year.var, int.var)
  if (length(full)) full <- orderDataByIntYear(full, year.var, int.var)

  list(predicted=pred, fullpredicted=full, 
       fit.uncond=fit.uncond, fit.cond=fit.cond)
}

