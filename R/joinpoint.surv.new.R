aapc <- function(fit, type="AbsChgSur", interval=5) {
  Z0<-NULL
  ACS.range<-NULL
#calculate for each segment: AAPC measure estimate and SE
#output columns: start.year end.year     estimate    std.error PredInterval
#For SE, calculate the derivative of aapc(bete,gamma) by approximation
#Z0: only 1 row
	stopifnot(type %in% c("RelChgHaz","AbsChgSur","RelChgSur")); #Removed "HAZ_AC(CS)", "HAZ_APC(CS)",; 
	if(!is.null(Z0)) Z0 = t(data.frame(as.vector(Z0)));
	if(!is.null(ACS.range)){
	  if(length(ACS.range)!=2){
	    stop("ACS.range should be defined using the minimum and maximum of the year range for Average Absolute Change in Survival trend measure.")
	  }
	  if(ACS.range[1]>=ACS.range[2]){
	    stop("The maximum should be defined greater than the minimum in the ACS.range.")
	  }
	  if(ACS.range[1]<min(fit$apc$start.year) | ACS.range[2]>max(fit$apc$end.year)){
	    stop("ACS.range should be within the entire range.")
	  }
	}
	fup = interval;
	getAapc = function(fit, beta, gamma, Z0) {
		epsilon = 1e-5;
		getAapc1 = function(years) {
			#pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			#pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
                     pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+epsilon,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
                     pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years-epsilon,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

			#deriv = (log(pred1$pred_cum) - log(pred2$pred_cum)) / epsilon / 2;
			deriv = ((pred1$pred_cum) - (pred2$pred_cum)) / epsilon / 2;
			return(deriv);
		}
		nSeg = dim(fit$apc)[1];
		newApc = rep(NA, nSeg);
		for (i in 1:nSeg) {
			t0 = fit$apc$start.year[i];
			t1 = fit$apc$end.year[i]-1;
			newApc[i] = mean(getAapc1(t0:t1));
		}
		return(newApc);
	}

	getAapc_log = function(fit, beta, gamma, Z0) {
		epsilon = 1e-5;
		getAapc1 = function(years) {
			#pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			#pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
                     pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+epsilon,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
                     pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years-epsilon,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

			deriv = (log(pred1$pred_cum) - log(pred2$pred_cum)) / epsilon / 2;
			#deriv = ((pred1$pred_cum) - (pred2$pred_cum)) / epsilon / 2;
			return(deriv);
		}
		nSeg = dim(fit$apc)[1];
		newApc = rep(NA, nSeg);
		for (i in 1:nSeg) {
			t0 = fit$apc$start.year[i];
			t1 = fit$apc$end.year[i]-1;
			newApc[i] = mean(getAapc1(t0:t1));
		}
		return(newApc);
	}

	getAapc_hr = function(fit, beta, gamma, Z0) {

    nJP = length(as.vector(fit$jp));
		result = rep(NA, nJP + 1);
		result[1] = beta[nJP + 1];
		if (nJP > 0) for (i in 1:nJP) result[i + 1] = result[i] + beta[i]
		return(exp(result) - 1);
	}

	getAVEAAPC = function(fit, beta, gamma, Z0) {
		getAapc1 = function(years) {
			#pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			#pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
                     pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+1,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)


			AAPC = ((pred2$pred_cum) - (pred1$pred_cum));
			return(AAPC);
		}
		nSeg = dim(fit$apc)[1];
		newApc = rep(NA, nSeg);
		for (i in 1:nSeg) {
			t0 = fit$apc$start.year[i];
			t1 = fit$apc$end.year[i]-1;
			newApc[i] = mean(getAapc1(t0:t1));
		}
		return(newApc);
	}
	getACSweight<-function(x_l1,x_l2,tau_k1,tau_k2){
	  if(tau_k1>=x_l2 | x_l1>=tau_k2){
	    weight<-0
	  }else{
	    I_start<-max(x_l1,tau_k1)
	    I_end<-min(x_l2,tau_k2)
	    weight<-I_end-I_start
	  }
	  return(weight)
	}
	getAVEAAPC_range = function(fit, beta, gamma, Z0) {
	  getAapc1 = function(years) {
	    #pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
	    #pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
           pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years,
                            intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
           pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+1,
                            intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

	    AAPC = ((pred2$pred_cum) - (pred1$pred_cum))
	    return(AAPC)
	  }
	  nSeg = dim(fit$apc)[1]
	  ACS.seg = rep(NA, nSeg)
	  w.seg<-rep(NA, nSeg)
	  
	  for (i in 1:nSeg) {
	    tau_k1 = fit$apc$start.year[i]
	    tau_k2 = fit$apc$end.year[i]
	    ACS.seg[i] = mean(getAapc1(tau_k1:(tau_k2-1)))
	    w.seg[i]=getACSweight(ACS.range[1],ACS.range[2],tau_k1,tau_k2)
	  }
	  
	  wACS<-sum(w.seg*ACS.seg)/sum(w.seg)
	  return(wACS)
	}

	getAVEARPC = function(fit, beta, gamma, Z0) {
		getAapc1 = function(years) {
			#pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			#pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
                     pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
                     pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+1,
                                      intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

			ARPC=0;
			zero.pos<-which(pred1$pred_cum==0)
			nonzero.pos<-which(pred1$pred_cum!=0)
			if(length(zero.pos)==0){
			  ARPC = ((pred2$pred_cum) - (pred1$pred_cum))/(pred1$pred_cum);
			}else{
			  ARPC[nonzero.pos]<-((pred2$pred_cum[nonzero.pos]) - (pred1$pred_cum[nonzero.pos]))/(pred1$pred_cum[nonzero.pos]);
      }
			return(ARPC);
		}
		nSeg = dim(fit$apc)[1];
		newApc = rep(NA, nSeg);
		for (i in 1:nSeg) {
			t0 = fit$apc$start.year[i];
			t1 = fit$apc$end.year[i]-1;
			newApc[i] = mean(getAapc1(t0:t1));
		}
		return(newApc);
	}

	result = NULL;

	#if (type == "HAZ_AC(CS)") result = .getSeAapc(fit, getAapc,Z0)
	#else if (type == "HAZ_APC(CS)") result = .getSeAapc(fit, getAapc_log,Z0)
	if (type == "RelChgHaz") result = .getSeAapc(fit, getAapc_hr,Z0,ACS.range=NULL)
	else if (type == "AbsChgSur" & is.null(ACS.range)) result = .getSeAapc(fit, getAVEAAPC,Z0,ACS.range=NULL)
	else if (type == "AbsChgSur" & !is.null(ACS.range)){
	  result = .getSeAapc(fit, getAVEAAPC,Z0,ACS.range=NULL)
	  result.range = .getSeAapc(fit, getAVEAAPC_range,Z0,ACS.range)
	  result<-list(result,result.range)
	} else if (type == "RelChgSur") result = .getSeAapc(fit, getAVEARPC,Z0,ACS.range=NULL) 
	return(result);

} # END: aapc

aapc.multiints <- function(fit, type="AbsChgSur",int.select=NULL,ACS.range=NULL,ACS.out=NULL){

  int0           <- int.select
  int.flag       <- 0
  shift.interval <- 0

  # 2024-09-03 Check for shift.interval in fit to shift int.select
  if (length(int.select)) {
    shift.interval  <- fit[["shift.interval", exact=TRUE]]
    if (length(shift.interval) == 1) {
      int.select <- int.select - shift.interval
      int.flag   <- 1
    }
    tmp <- !(int.select %in% fit$Interval)
    if (any(tmp)) {
      err <- paste0(sort(int0[tmp]), collapse=", ") 
      msg <- paste0("int.select = ", err, " is/are not valid")
      stop(msg) 
    }  
  }

  Z0<-NULL
  #calculate for each segment: AAPC measure estimate and SE
  #output columns: 4 start.year end.year     estimate    std.error PredInterval
  #For SE, calculate the derivative of aapc(bete,gamma) by approximation
  #Z0: only 1 row
  stopifnot(type %in% c("RelChgHaz","AbsChgSur","RelChgSur")); #Removed "HAZ_AC(CS)", "HAZ_APC(CS)",; 
  stopifnot(ACS.out %in% c(NULL,"user","both"))
  if(!is.null(Z0)) Z0 = t(data.frame(as.vector(Z0)));
  if(!is.null(ACS.range)){
    if(type!="AbsChgSur"){
      stop("ACS.range needs to be defined for trend measure type='AbsChgSur' only.")
    }
    if(length(ACS.range)!=2){
      stop("ACS.range should be defined using the minimum and maximum of the year range for Average Absolute Change in Survival trend measure.")
    }
    if(ACS.range[1]>=ACS.range[2]){
      stop("The maximum should be defined greater than the minimum in the ACS.range.")
    }
    if(ACS.range[1]<min(fit$apc$start.year) | ACS.range[2]>max(fit$apc$end.year)){
      stop("ACS.range should be within the entire range.")
    }
    if(is.null(ACS.out)){
      stop("ACS.out should be either 'user' or 'both' when the year range ACS.range is defined.")
    }
  }else if(is.null(ACS.range)){
    if(!is.null(ACS.out)){
      stop("ACS.out should be NULL when ACS.range is NULL.")
    }
  }

  getAapc = function(fit, beta, gamma, Z0) {
    epsilon = 1e-5;
    getAapc1 = function(years) {
      #pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      #pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+epsilon,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
      pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years-epsilon,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

      #deriv = (log(pred1$pred_cum) - log(pred2$pred_cum)) / epsilon / 2;
      deriv = ((pred1$pred_cum) - (pred2$pred_cum)) / epsilon / 2;
      return(deriv);
    }
    nSeg = dim(fit$apc)[1];
    newApc = rep(NA, nSeg);
    for (i in 1:nSeg) {
      t0 = fit$apc$start.year[i];
      t1 = fit$apc$end.year[i]-1;
      newApc[i] = mean(getAapc1(t0:t1));
    }
    return(newApc);
  }
  
  getAapc_log = function(fit, beta, gamma, Z0) {
    epsilon = 1e-5;
    getAapc1 = function(years) {
      #pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      #pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+epsilon,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
      pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years-epsilon,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

      deriv = (log(pred1$pred_cum) - log(pred2$pred_cum)) / epsilon / 2;
      #deriv = ((pred1$pred_cum) - (pred2$pred_cum)) / epsilon / 2;
      return(deriv);
    }
    nSeg = dim(fit$apc)[1];
    newApc = rep(NA, nSeg);
    for (i in 1:nSeg) {
      t0 = fit$apc$start.year[i];
      t1 = fit$apc$end.year[i]-1;
      newApc[i] = mean(getAapc1(t0:t1));
    }
    return(newApc);
  }
  
  getAapc_hr = function(fit, beta, gamma, Z0) {
 
    nJP = length(as.vector(fit$jp));
    result = rep(NA, nJP + 1);
    result[1] = beta[nJP + 1];
    if (nJP > 0) for (i in 1:nJP) result[i + 1] = result[i] + beta[i]
    return(exp(result) - 1);
  }
  
  getAVEAAPC = function(fit, beta, gamma, Z0) {

    getAapc1 = function(years) {

      #pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      #pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);    
      pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
      pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+1,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
  
      AAPC = ((pred2$pred_cum) - (pred1$pred_cum));
      return(AAPC);
    }
    nSeg = dim(fit$apc)[1];
    newApc = rep(NA, nSeg);
    for (i in 1:nSeg) {
      t0 = fit$apc$start.year[i];
      t1 = fit$apc$end.year[i]-1;
      newApc[i] = mean(getAapc1(t0:t1));
    }
    return(newApc);
  }
  # get ACS weights for each jp segments [tau_k1,tau_k2] over the user-specified range [x_l1,x_l2]
  getACSweight<-function(x_l1,x_l2,tau_k1,tau_k2){
    if(tau_k1>=x_l2 | x_l1>=tau_k2){
      weight<-0
    }else{
      I_start<-max(x_l1,tau_k1)
      I_end<-min(x_l2,tau_k2)
      weight<-I_end-I_start
    }
    return(weight)
  }
  ### the Absolute change in survival between 2 calendar years the users select
  getAVEAAPC_range = function(fit, beta, gamma, Z0) {
    getAapc1 = function(years) {
      #pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
      #pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
      pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma) 
      pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+1,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma) 

      AAPC = ((pred2$pred_cum) - (pred1$pred_cum))
      return(AAPC)
    }
    nSeg = dim(fit$apc)[1]
    ACS.seg = rep(NA, nSeg)
    w.seg<-rep(NA, nSeg)
    
    for (i in 1:nSeg) {
      tau_k1 = fit$apc$start.year[i]
      tau_k2 = fit$apc$end.year[i]
      ACS.seg[i] = mean(getAapc1(tau_k1:(tau_k2-1)))
      w.seg[i]=getACSweight(ACS.range[1],ACS.range[2],tau_k1,tau_k2)
    }
    
    wACS<-sum(w.seg*ACS.seg)/sum(w.seg)
    return(wACS)
  }
  
  getAVEARPC = function(fit, beta, gamma, Z0) {
    getAapc1 = function(years) {
      #pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      #pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred1 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)
      pred2 <- Predict(fit, fit$Year, fit$proj.year.num, fit$Interval, fit$jp, years=years+1,
                       intervals=fup, beta_input=beta, Z0=Z0, gamma_input=gamma)

      ARPC=rep(0,length(pred1$pred_cum))
      zero.pos<-which(pred1$pred_cum==0)
      nonzero.pos<-which(pred1$pred_cum!=0)
      if(length(zero.pos)==0){
        ARPC = ((pred2$pred_cum) - (pred1$pred_cum))/(pred1$pred_cum);
      }else{
        ARPC[nonzero.pos]<-((pred2$pred_cum[nonzero.pos]) - (pred1$pred_cum[nonzero.pos]))/(pred1$pred_cum[nonzero.pos]);
      }
      return(ARPC);
    }
    nSeg = dim(fit$apc)[1];
    newApc = rep(NA, nSeg);
    for (i in 1:nSeg) {
      t0 = fit$apc$start.year[i];
      t1 = fit$apc$end.year[i]-1;
      newApc[i] = mean(getAapc1(t0:t1));
    }
    return(newApc);
  }
  
  aapc.results<-list()
  if(type == "AbsChgSur" & !is.null(ACS.range)){
    aapc.results.range<-list()
  }

  for(i in 1:length(int.select)){
    fup = int.select[i];

    result = NULL
    result.range=NULL
    #if (type == "HAZ_AC(CS)") result = .getSeAapc(fit, getAapc,Z0)
    #else if (type == "HAZ_APC(CS)") result = .getSeAapc(fit, getAapc_log,Z0)

    if (type == "RelChgHaz"){
      result = .getSeAapc(fit, getAapc_hr,Z0,ACS.range=NULL)
    } else if (type == "AbsChgSur" & is.null(ACS.range)){
      result = .getSeAapc(fit, getAVEAAPC,Z0,ACS.range=NULL)
    } else if (type == "AbsChgSur" & !is.null(ACS.range)){
      if (ACS.out=="user"){
        result.range = .getSeAapc(fit, getAVEAAPC_range,Z0,ACS.range)
      }else if (ACS.out=="both"){
        result = .getSeAapc(fit, getAVEAAPC,Z0,ACS.range=NULL)
        result.range = .getSeAapc(fit, getAVEAAPC_range,Z0,ACS.range)
      }
    } else if (type == "RelChgSur") {
      result = .getSeAapc(fit, getAVEARPC,Z0,ACS.range=NULL) 
    }

    if(!is.null(result)){
      result[,"interval"]<-fup
      aapc.results[[i]]<-result
    }
    if(!is.null(result.range)){
      result.range[,"interval"]<-fup
      aapc.results.range[[i]]<-result.range
    }
  }
  if(type == "AbsChgSur" & !is.null(ACS.range)){
    if(ACS.out=="both"){
      aapc.results<-list(aapc.results,aapc.results.range)
      names(aapc.results)<-c("ACS.jp","ACS.user")
    }else if(ACS.out=="user"){
      aapc.results <- aapc.results.range
    }
  }
  if (int.flag) {
    tmp <- try(aapc.update.interval(aapc.results, shift.interval), silent=TRUE)
    if (!inherits(tmp, "try-error")) aapc.results <- tmp
  }

  return(aapc.results)

} # END: aapc.multiints

aapc.update.interval <- function(retlist, shift.interval) {

  # retlist could be  list of data frames or a list of sublists of data frames
  n <- length(retlist)
  if (!n || !is.list(retlist)) return(retlist)
  for (i in 1:n) {
    x <- retlist[[i]]
    if (is.data.frame(x)) {
      if ("interval" %in% colnames(x)) {
        x[, "interval"] <- x[, "interval", drop=TRUE] + shift.interval
      }
    } else if (is.list(x)) {
      x <- aapc.update.interval(x, shift.interval)
    }
    retlist[[i]] <- x
  } 
  retlist
}

# compute the predicted cumulative survival rates.
Predict <- function(cox_fit, Year, proj.year.num, Interval, jpoints, years=NULL, intervals=NULL,
                    beta_input=NULL, Z0=NULL, gamma_input=NULL, ret.jacob=FALSE) {

  # ret.jacob : TRUE of FALSE to return Jacobian matrices for a single year
  
  if (ret.jacob && length(years) != 1) stop("ERROR: length(years) must be 1 when ret.jacob is TRUE")

  nJP <- length(jpoints)
  if (is.null(years)) years = min(Year):(max(Year) + proj.year.num)
  if (is.null(intervals)) intervals = (1:max(Interval))
  if (is.null(Z0)) Z0 = cox_fit$Z0
  if (is.null(beta_input)) beta_input = cox_fit$coefficients[1:cox_fit$p, 1]
  if (cox_fit$numparaCovar > 0) {
    if (is.null(gamma_input)) gamma_input = cox_fit$coefficients[cox_fit$p+1:cox_fit$numparaCovar, 1]
  } else {
    gamma_input=NULL
  }

  maxInterval = max(intervals)
  stopifnot(maxInterval <= max(Interval))
		
  beta2 = list()
  #cox_fit$beta = cox_fit$coefficients[, 1]
  #intercept = cox_fit$beta[nJP + 1]
  intercept = 0
  if (nJP > 0) beta2$Seg = beta_input[1:nJP]
  beta2$Year = beta_input[nJP + 1]
  #beta2$Interval = c(0, beta_input[(nJP+2):length(beta_input)])
  #beta2$Interval = beta2$Interval + intercept
  beta2$Interval = beta_input[(nJP+2):cox_fit$p]
  if (cox_fit$numparaCovar > 0) {
    beta2$betagamma = c(beta_input,gamma_input)
  } else {
    beta2$betagamma = c(beta_input)
  }
			
  nrowcovars=1
  if (!is.null(Z0)) {
    nrowcovars = nrow(Z0)
  } else {
    Z0=0
  }
  result = data.frame(matrix(0, length(years) * maxInterval * nrowcovars, 4+2+length(cox_fit$covarnames)));
  names(result) = c("Year", "Interval", "pred_int", "pred_cum","pred_int_se","pred_cum_se",cox_fit$covarnames);

  if (ret.jacob) {
    jacobmat.int           <- matrix(data=NA, nrow=maxInterval, ncol=length(beta2$betagamma))
    colnames(jacobmat.int) <- names(beta2$betagamma)
    jacobmat.cum           <- jacobmat.int 
  }

  for (ZIndex in 1:nrowcovars){
    surv_int   = matrix(NA, length(years), maxInterval)
    surv_cum   = surv_int
    surv_xbeta = surv_int
    z0         = numeric(0)
    zgamma     = 0
    if (!is.null(gamma_input) && cox_fit$numparaCovar>0) {
      zgamma=Z0[ZIndex,,drop=FALSE] %*% gamma_input
      z0=Z0[ZIndex,,drop=FALSE]
    }
  
    for (i in 1:length(years)) {
      for (j in 1:maxInterval) {
        xbeta2 = beta2$Interval[j] + years[i] * beta2$Year
        if (nJP > 0) for (k in 1:nJP) {
                       xseg = years[i] - jpoints[k]
                       if (xseg < 0) xseg = 0
                       xbeta2 = xbeta2 + beta2$Seg[k] * xseg
                     }
        surv_xbeta[i, j] = xbeta2+zgamma
        surv_int[i, j] = exp(-exp(xbeta2+zgamma))
        if (j == 1) surv_cum[i, j] = surv_int[i, j]
  	 else surv_cum[i, j] = surv_cum[i, j-1] * surv_int[i, j]
  	 idx = (ZIndex-1)*length(years)*maxInterval + (i - 1) * maxInterval + j
        result[idx, 1] = years[i]
        result[idx, 2] = j
        result[idx, 3] = surv_int[i, j]
        result[idx, 4] = surv_cum[i, j]

        DerivPredInt = .calcDerivPredInt(cox_fit,jpoints,beta2$betagamma,yearvalue=years[i],intervalvalue=j,z0=z0)
  	result[idx, 5] = sqrt(t(DerivPredInt)%*% cox_fit$covariance %*%DerivPredInt)
        DerivPredCum = .calcDerivPredCum(cox_fit,jpoints,beta2$betagamma,yearvalue=years[i],intervalvalue=j,z0=z0)
        result[idx, 6] = sqrt(t(DerivPredCum)%*% cox_fit$covariance %*%DerivPredCum)
        if(cox_fit$numparaCovar>0){
          result[idx, 7:ncol(result)] = Z0[ZIndex,,drop=TRUE] 
        }

        if (ret.jacob) {
          jacobmat.int[j, ] <- DerivPredInt
          jacobmat.cum[j, ] <- DerivPredCum
        }
      }
    }
  }

  if (ret.jacob) return(list(jacob.int=jacobmat.int, jacob.cum=jacobmat.cum))

  #result1 = subset(result, Interval %in% intervals)
  requested = match(result$Interval, intervals,0L)

  result1 = result[requested>0,]

  return(result1)

} # END: Predict

GetApc <- function(cox_fit, jpoints, Year) {

  nJP    <- length(jpoints)
  result = matrix(NA, nJP+1, 3)
  result = data.frame(result)
  dimnames(result)[[2]] = c("start.year", "end.year", "estimate")
  #dimnames(result)[[1]] = paste("Seg", 1:(nJP+1), sep = " ")
  rownames(result) = NULL
  result$estimate[1] = cox_fit$coefficients[nJP+1, 1]
  result$start.year[1] = min(Year)
  result$end.year[nJP + 1] = max(Year)
  if (nJP > 0) {
    for (i in 1:nJP) {
      result$estimate[i + 1] = result$estimate[i] + cox_fit$coefficients[i, 1]
      result$end.year[i] = jpoints[i]
      result$start.year[i + 1] = jpoints[i]
    }
  }
  return(result)

} # END: GetApc

# Indv_bic() computes the BIC for the given join points.
Indv_bic <- function(jpoints,Year,nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                     proj.year.num) {

  lineSeg = NULL
  nJP = length(jpoints)
  if (nJP > 0) {
    for (i in 1:nJP) {
      seg = Year - jpoints[i]
      seg[seg < 0] = 0
      lineSeg = cbind(lineSeg, seg)
    }
    colnames(lineSeg) = paste("jp", 1:nJP, sep = "_")
  }
		
  X1 = cbind(lineSeg, X)
  cox_fit = .CoxFit(X1, nAlive, nDied, nLost, ExpSurv,Z)
		
  # run checks on proj.year.num
  if (is.null(proj.year.num)) { 
    proj.year.num = 5
  } else{
    if (proj.year.num<0 | proj.year.num>30) proj.year.num = 5
  }
		
  # compute the predicted cumulative survival rates.
  #res.predict <- Predict(cox_fit, Year, proj.year.num, Interval, jpoints, years=NULL, intervals=NULL,
  #                       beta_input=NULL, Z0=NULL, gamma_input=NULL)
  cox_fit$Predict <- Predict		
  cox_fit$apc     <- GetApc(cox_fit, jpoints, Year)
  
  return(cox_fit)

} # END: Indv_bic

.Handle.op=function(x,default=3){
    if(is.null(x)){
      res=default;
    }else{
      res=x;
    }
    return(res);  
}

.Get.integerlist=function(n,multiplier){
    if(n>0){
      res=(1:n)*multiplier;
    }else{
      res=integer(0);
    }
    return(res);
}

.GetJP <- function(nJP, op, Year, nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                     proj.year.num) {

  if(nJP>0){
    #op=list(), #options 
    #             numbetwn: number of skipped obs between joinpoints exclusive (not count for the joinpoints);
    #             numfromstart: number of skipped obs from the first obs to joinpoints exclusive (not count for the joinpoint);
    #             numtoend: number of skipped obs from the first obs to joinpoints exclusive (not count for the joinpoint);
    op$numbetwn=.Handle.op(op$numbetwn,2)
    op$numfromstart=.Handle.op(op$numfromstart,3)
    op$numtoend=.Handle.op(op$numtoend,5)

    intervalSize = op$numbetwn+1;
    #what implemented before:  so called intervalSize fixed as 3
    #   initial joinpoints: numfromstart=3, numbetwn=2, numtoend=5
    #   find next: numbetwn=2, numtoend=3
    #   exit: numtoend=5

    stopifnot(max(Year) >= min(Year) + op$numfromstart + (nJP-1) * intervalSize + op$numtoend)# here it is consistent with what's claimed for op 
    jpoints = min(Year) + op$numfromstart + c(0,.Get.integerlist(nJP-1, intervalSize))# here it is consistent with what's claimed for op 
    endFlag = FALSE;
    #numtoend.temp=op$numtoend-1;#It is not consistent with what's claimed for op, in future replace op$numtoend for numtoend.temp below
    numtoend.temp=op$numtoend;# here it is consistent with what's claimed for op 

    nextjp = function(oldjp) {
      oldjp[nJP] = oldjp[nJP] + 1;
      if (nJP >= 2) for (k in nJP:2) {
		        if (oldjp[k] > max(Year) - intervalSize * (nJP - k) - numtoend.temp) {
                        oldjp[k - 1] = oldjp[k - 1] + 1
                        oldjp[k:nJP] = oldjp[k - 1] + intervalSize * (1:(nJP - k + 1))
                      }
                    }
      if (oldjp[1] > max(Year) - intervalSize * (nJP-1)-op$numtoend) endFlag <<- TRUE# here it is consistent with what's claimed for op 
      return(oldjp)
    }

    result = list()
    result$jp = NA
    result$bic = Inf
    nIter = 0
    jp.debug = FALSE  # setting the debug option
    while (!endFlag) {
      new_bic <- Indv_bic(jpoints,Year,nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                          proj.year.num)
      if (new_bic$bic < result$bic) {
        result = new_bic
        result$jp = jpoints
      }
      jpoints = nextjp(jpoints)
      nIter = nIter + 1
      if (jp.debug & nIter > 3) break
    }
  } else {
    result <- Indv_bic(NULL, Year,nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                       proj.year.num)
    result$jp = NULL
  }

  # Add objects to return list
  result$Year          <- Year
  result$proj.year.num <- proj.year.num
  result$Interval      <- Interval

  return(result)

} # END: .GetJP

.GetBestJP = function(nJP, op, Year, nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                      proj.year.num) {

  FitList=list()
  bestFit <- .GetJP(0, op, Year, nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                     proj.year.num)
  FitList[[length(FitList)+1]]=bestFit
  if (nJP > 0) {
    for (k in 1:nJP) {
      message(paste0("Computing estimates when number of join points is ", k))
      newFit = .GetJP(k, op, Year, nAlive, nDied, nLost, ExpSurv,X,Z,Interval,
                      proj.year.num)
      FitList[[length(FitList)+1]]=newFit
      if (newFit$bic < bestFit$bic) bestFit = newFit
    }
  }
  return(list(bestFit=bestFit, FitList=FitList))

} # END: .GetBestJP

#
# R package JPSurv is for fitting the joinpoint survival model
#
# Key components of codes: 
# 1) accessory functions 
# 2) function joinpoint.surv(...) used for fitting the model
#
# joinpoint arguments:
# data,                  # data set names
# subset=NULL,               # subset option
# na.action = na.fail,
# year="Year",         # year variable
# interval="Interval",             # interval survival times
# number.event="Died",                # number of death
# number.alive="Alive_at_Start",     # number of alive (alive in the beginning of interval)
# number.loss="Lost_to_Followup",    # number of people lost to follow-up (if this variable does not exist, it means 0)
# expected.rate="Expected_Survival_Interval", # expected interval survival rates: Expected_Survival_Interval in seer data
# observedrelsurv = NULL, # Relative_Survival_Cum in seer data
# model.form = NULL,     # =~-1+age+as.factor(stage)
# maxnum.jp = 0,
# proj.year.num = NULL,   #Added 8/24/2015 DM - Number of projection years, range 0-30, default 5 set below in code;
# op=list(), #options 
#             numbetwn: number of skipped obs between joinpoints exclusive (not count for the joinpoints);
#             numfromstart: number of skipped obs from the first obs to joinpoints exclusive (not count for the joinpoint);
#             numtoend: number of skipped obs from the first obs to joinpoints exclusive (not count for the joinpoint);
# delLastIntvl=F
joinpoint <- function(data, subset=NULL, na.action = na.fail, year="Year", interval="Interval", 
                      number.event="Died", number.alive="Alive_at_Start",number.loss="Lost_to_Followup", 
                      expected.rate="Expected_Survival_Interval", observedrelsurv = NULL, model.form = NULL,
                      maxnum.jp = 0, proj.year.num = 5, # Edited 04/07/2020 FZ;Added 8/24/2015 DM - Number of projection years, range 0-30, default 5 set below in code
                      op=list(), delLastIntvl=FALSE, add.data.cols="_ALL_"){
# About input data:
#
# 1) year,interval, number.event, number.alive, number.loss, expected.rate, observedrelsurv: 
#     could be a vector of numeric or a character string giving a column name of the argument 'data'
# 2) interval should be integers, not a string like "1-<2 yr"
# 3) model.form is used for defining covariates and it is a object of "formula". 
#    The associated data are contained in the argument 'data'.
# 4) Rows of year, interval, number.event, number.alive, number.loss, expected.rate, or observedrelsurv if given as numeric vectors
#    should match, in the same order, rows of 'data' without subsetting, i.e. subset will be applied on those
#    arguments too.
#input data structure:
#rows: records
#columns:
#   must contain: 
#      for survival info: 
#             (alive, died (number.event), lost, expected, [observedsurvival])
#                  the meaning of expected = ExpectedSurvAtEnd/ExpectedSurvAtStart = ExpectedSurvInterval
#             [] means optional
#             observedsurvival is observed Relative_Survival_Cum
#
#   optional:
#      covariates in Z specified by model.form: use -1+ to remove the intercept

#major steps:
#1) handle arguments
#2) prepare fitting inputs
#  survMatrix: columns
#    "year","interval",
#    "number.event","number.alive", "number.loss", "expected.rate",
#    "observedrelsurv"
#             the meaning of expected = ExpectedSurvAtEnd/ExpectedSurvAtStart = ExpectedSurvInterval
#
#
#   covarMatrix of Z from model.form, on the right of ~
#   
#3) fit the model and return results
#      method of fitting: maximum likelihood esitimation by iteratively reweighted Least square method

  .merge.formula <- function(form1, form2, ...){
  
  	# get character strings of the names for the responses 
  	# (i.e. left hand sides, lhs)
  	lhs1 <- as.character(as.list(form1))[1];
  	lhs2 <- as.character(as.list(form2))[1];
  	if(lhs1 != lhs2) stop('both formulas must have the same response')
  
  	# get character strings of the right hand sides
  	rhs1 <- as.character(as.list(form1))[2];
  	rhs2 <- as.character(as.list(form2))[2];
  
  	# create the merged rhs and lhs in character string form
  	opr="";
  	if(substr(rhs2,1,1)!="-") opr="+";
  	rhs <- paste(rhs1,opr ,rhs2, sep="")
  	lhs <- as.character("~")
  	formula.str=paste(lhs,rhs,sep="");
  	# put the two sides together with the amazing 
  	# reformulate function
  	out <- as.formula(formula.str);
  
  	# set the environment of the formula (i.e. where should
  	# R look for variables when data aren't specified?)
  	environment(out) <- parent.frame()
  
  	return(out)
  }

  # Check expected.rate column 2023-05-26
  if (!length(expected.rate)) expected.rate <- "Expected_Survival_Interval"
  if (!(expected.rate %in% colnames(data))) data[, expected.rate] <- 1

  # Argument add.data.cols added on 2023-05-18
  add.data.cols <- check_add.data.cols(add.data.cols)
  
  # Changed on 2022-10-04
  # Set observedrelsurv to NULL since it is not used in the code. If non-NULL, then
  #  rows of data could mistakenly be removed due to missing values with observedrelsurv.
  observedrelsurv <- NULL

#1) handle arguments
  mfcall <- match.call(expand.dots=FALSE);
  #mfcall.datainfo.index <- match(c("data", "subset", "na.action"), names(mfcall), 0L);
  mfcall.datainfo.index <- match(c("data", "eval(parse(text=subset))", "na.action"), names(mfcall), 0L)  ### FZ 07/16/2019
  mfcall.datainfo <- mfcall[c(1L, mfcall.datainfo.index)];
  mfcall.datainfo$drop.unused.levels <- TRUE;
  #inputdata=mfcall.datainfo$data;
  
# extract data of survinfo: year, interval, number.event, number.alive, number.loss, expected.rate, observedrelsurv
  #crntdata=data;
  if(is.null(subset)){
    crntdata<-data
  }else{
    crntdata<-subset(data,eval(parse(text=subset))) ### FZ 07/17/2019
  }
  if(dim(crntdata)[1]==0){
    stop.str<-paste("No data available for the cohort selection '",subset,"' in the input data file.",sep="")
    stop(stop.str)
  }
  classtypes = c(
    class(year), class(interval), 
    class(number.event), class(number.alive), class(number.loss), class(expected.rate),
    class(observedrelsurv)
  );
  allvarnames=c(
    "year","interval",
    "number.event","number.alive", "number.loss", "expected.rate",
    "observedrelsurv"
  );
  classtypes.index <- match(classtypes, c("character"), 0L);

  # R CMD check warning BW 12/03/2020
  #eval(parse(text=paste("varfromargs=c(",paste0(allvarnames[classtypes.index==1],collapse=", "),")",sep="")),environment()); #those are given in args and assumed to be the colnames of data
  varfromargs <- eval(parse(text=paste("c(",paste0(allvarnames[classtypes.index==1],collapse=", "),")",sep="")),environment()); #those are given in args and assumed to be the colnames of data
  
  if(length(varfromargs)>0) allvarnames[classtypes.index==1]=varfromargs;
  varnotfromargs=allvarnames[classtypes.index==0 & classtypes != "NULL"];
  varnotfromargs.index <- match(varnotfromargs, colnames(crntdata), 0L);
  if(length(varnotfromargs)>0)for(Index in 1:length(varnotfromargs)){
    if(varnotfromargs.index[Index]>0){
      exprstr=paste("crntdata[,",as.character(varnotfromargs.index[Index]),"]=",varnotfromargs[Index],sep="");
      eval(parse(text=exprstr));
    }else{
      exprstr=paste("crntdata=cbind(crntdata,",varnotfromargs[Index],"=",varnotfromargs[Index],");",sep="");
      eval(parse(text=exprstr));
    }
  }
  varnames=allvarnames[classtypes != "NULL"];
  crntdata=.covert.numeric(crntdata,varnames);
  crntformula=paste("~",paste0(varnames,collapse="+"),sep="");

  if(delLastIntvl){
    crntdata=deleteLastInterval(data=crntdata,byvarnames=c(),yearcol=varnames[1],intervalcol=varnames[2]);
  }
  mfcall.datainfo$data=as.name("crntdata");
  mfcall.datainfo$formula = as.formula(crntformula);
  mfcall.datainfo[[1]] <- as.name("model.frame")
  m.survinfo <- eval(mfcall.datainfo, list(sys.frame(sys.parent()),environment()));#
  mterms.survinfo <- attr(m.survinfo,"terms")
  survMatrix = model.frame(mterms.survinfo, m.survinfo);
  NumsurvVars=dim(survMatrix)[2];

# extract data of covarmatrix: covariates (model.form)
  covarMatrix=NULL;
  m.covars=NULL;
  covar.names=NULL;
  referencelevel=NULL;
 
  if(!is.null(model.form)){
    #mfcall.datainfo$data=inputdata;
    mfcall.datainfo$data=as.name("crntdata");
    mfcall.datainfo$formula = .merge.formula(model.form,as.formula(crntformula));
    mfcall.datainfo[[1]] <- as.name("model.frame")
    #m.covars <- eval(mfcall.datainfo, sys.frame(sys.parent()));
    m.covars <- eval(mfcall.datainfo, list(sys.frame(sys.parent()),environment()));
    mterms.covars <- attr(m.covars,"terms")
    covarMatrix = model.matrix(mterms.covars,m.covars);
    covar.names=as.character(attr(mterms.covars,"variables"));
    covar.names=covar.names[2:(length(covar.names)-NumsurvVars)];

    .get.reference=function(varnames, data){
      res=character(length(varnames));
      cols.data=colnames(data);
      if(length(cols.data)>0) for(index.cols.data in 1:length(cols.data)){
        index.in.varnames = .search.substr(cols.data[index.cols.data],varnames,fixed=TRUE);
        if(!is.null(index.in.varnames) && (length(index.in.varnames)>0)){
          if(length(index.in.varnames)>1) stop("colnames repeat in covariate list!");
          res[index.in.varnames]=levels(factor(as.character(data[,index.cols.data])))[1];
        }
      }
      return(res);
    }
    referencelevel=.get.reference(covar.names, crntdata);

    ColStartSurv=dim(covarMatrix)[2]-NumsurvVars+1;

    survMatrix=covarMatrix[,c(ColStartSurv:(dim(covarMatrix)[2])),drop=FALSE];

  }

  nAlive = survMatrix[, 4];
  nDied = survMatrix[, 3];
  nLost = survMatrix[, 5];
  ExpSurv = survMatrix[, 6];
  Interval = survMatrix[, 2];
  Year = survMatrix[, 1];
  RelSurvCum = NULL;
  if (!is.null(observedrelsurv)) RelSurvCum = survMatrix[, 7];
  if(length(unique(Interval))==1){                   ## FZ 07/26/2019 fix the issue when interval=1
    Interval_<-Interval
  }else{
    Interval_ = as.factor(Interval);
  }
  X = model.matrix(~-1+Year+Interval_);
  if(!is.null(covarMatrix) && (2<=(ColStartSurv-1))){
    Z = covarMatrix[,c(2:(ColStartSurv-1)),drop=FALSE];
  }else{
    Z = NULL;
  }

  TempStartTime <- Sys.time();
  fitres = .GetBestJP(maxnum.jp, op, Year, nAlive, nDied, nLost, ExpSurv,X,Z,
                      Interval, proj.year.num);

  cox_fit=fitres$bestFit;
 	
  #parameters: alpha, joinpoints, beta, gamma
  #infomat: for alpha, beta, gamma
  #   infodimname
  #   para: alpha, beta, gamma
  #   se: alpha, beta, gamma
 	
  #prediction of survival: require (alpha, joinpoints, beta, gamma; year, interval (range only[1,J]), Z)
  #   se of survival

  TempEndTime <- Sys.time();
  TimeFittingmodel = difftime(TempEndTime,TempStartTime, units = c("secs"));


  covar.col.names=colnames(cox_fit$Z0);
  year.name.input=allvarnames[1];
  interval.input=allvarnames[2];
  input.merge=c(covar.col.names,year.name.input,interval.input);
  output.merge=c(covar.col.names,"Year","Interval");
  iZ0=NULL;
  if(!is.null(covar.col.names)){
    iZ0 = diag(length(covar.col.names)+1)[,1:length(covar.col.names)];
  }

  first.y = min(survMatrix[,1]);
  last.y = max(survMatrix[,1]);
  intnum=(last.y-first.y)+1;
  maxintnum=intnum-1;
  empty.triangle=NULL;
  for(y in first.y:last.y){
    currenty = cbind(y,1:intnum);
    empty.triangle = rbind(empty.triangle[,1:2],currenty);
    intnum=intnum-1;
  }
  colnames(empty.triangle)[1:2]=c("Year","Interval");
  #pred = cox_fit$Predict(Z0=iZ0);

  pred <- Predict(cox_fit, Year, proj.year.num, Interval, cox_fit$jp, years=NULL, intervals=NULL,
                    beta_input=NULL, Z0=iZ0, gamma_input=NULL)
  pred_tri = merge(empty.triangle,pred,by=c("Year","Interval"));

  cox_fit$predicted = merge(survMatrix, pred_tri, by.x=input.merge, by.y=output.merge, all.y=TRUE);
  cox_fit$predicted = cox_fit$predicted[order(cox_fit$predicted[,1],cox_fit$predicted[,2]),];
  cox_fit$fullpredicted = merge(survMatrix, pred, by.x=input.merge, by.y=output.merge, all.y=TRUE);
  cox_fit$fullpredicted = cox_fit$fullpredicted[order(cox_fit$fullpredicted[,1],cox_fit$fullpredicted[,2]),];

  BestNumJP=length(cox_fit$jp);
  FitList=fitres$FitList;
  for(Index in 0:maxnum.jp){
    if(BestNumJP == Index){
      FitList[[Index+1]]=cox_fit;
    }else{
      #pred = FitList[[Index+1]]$Predict(Z0=iZ0);
      pred <- Predict(FitList[[Index+1]], Year, proj.year.num, Interval, FitList[[Index+1]]$jp, 
                      years=NULL, intervals=NULL, beta_input=NULL, Z0=iZ0, gamma_input=NULL)

      pred_tri = merge(empty.triangle,pred,by=c("Year","Interval"));
      FitList[[Index+1]]$predicted = merge(survMatrix, pred_tri, by.x=input.merge, by.y=output.merge, all.y=TRUE);
      FitList[[Index+1]]$predicted = FitList[[Index+1]]$predicted[order(FitList[[Index+1]]$predicted[,1],FitList[[Index+1]]$predicted[,2]),];
      FitList[[Index+1]]$fullpredicted = merge(survMatrix, pred, by.x=input.merge, by.y=output.merge, all.y=TRUE);
      FitList[[Index+1]]$fullpredicted = FitList[[Index+1]]$fullpredicted[order(FitList[[Index+1]]$fullpredicted[,1],FitList[[Index+1]]$fullpredicted[,2]),];
      FitList[[Index+1]]$jp = rep(NULL,Index);
      for(j in 1:Index){
        if(Index > 0){
          FitList[[Index+1]]$jp[j] = FitList[[Index+1]]$apc[j,2];#Getting the joinpoints from the APC output for nonselected models. 
        }
      }
      #FitList[[Index+1]]$plot.surv = plot.surv; #Added to test plotting of FitList object.
    }	
  }

  cox_fit$FitList=FitList;
  #cox_fit$plot.surv = plot.surv;
  cox_fit$TimeFittingmodel=TimeFittingmodel;
  cox_fit$covar.names=covar.names;
  cox_fit$Z=Z;
  cox_fit$referencelevel=referencelevel;
  #cox_fit$survMatrix=survMatrix;
  #cox_fit$covarMatrix=covarMatrix;

  # For conditional survival
  cox_fit$interval     <- interval
  cox_fit$year         <- year
  cox_fit$number.alive <- number.alive
  cox_fit$number.event <- number.event
  cox_fit$number.loss  <- number.loss
  
  # Add columns to returned data frames
  if (length(add.data.cols)) {
    cox_fit$predicted     <- jp_addDataCols(cox_fit$predicted,     add.data.cols, crntdata, year, interval)
    cox_fit$fullpredicted <- jp_addDataCols(cox_fit$fullpredicted, add.data.cols, crntdata, year, interval)

    # Add to each fit in FitList
    FitList <- cox_fit$FitList
    for (i in 1:length(FitList)) {
      tmp               <- FitList[[i]]
      tmp$predicted     <- jp_addDataCols(tmp$predicted,     add.data.cols, crntdata, year, interval)
      tmp$fullpredicted <- jp_addDataCols(tmp$fullpredicted, add.data.cols, crntdata, year, interval)
      FitList[[i]]      <- tmp
    }
    cox_fit$FitList <- FitList
  }

  class(cox_fit) <- "joinpoint";
  return(cox_fit);

} # END: joinpoint

jp_addDataCols <- function(df, add, crntdata, year.var, int.var) {

  # Get all column names
  if (any(add %in% "_ALL_")) add <- colnames(crntdata)
  
  tmp  <- !(add %in% colnames(df)) & (add %in% colnames(crntdata))
  add  <- add[tmp]
  nadd <- length(add)
  if (!nadd) {
    warning("No columns from add.data.cols added to results")
    return(df)
  }

  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors=FALSE, check.names=FALSE)
  df.ids <- paste0(trimws(df[, year.var, drop=TRUE]), ":", trimws(df[, int.var, drop=TRUE]))
  if (any(duplicated(df.ids))) {
    warning("Non-unique year/interval pairs in results data")
    return(df)
  }
  cr.ids <- paste0(trimws(crntdata[, year.var, drop=TRUE]), ":", trimws(crntdata[, int.var, drop=TRUE]))
  if (any(duplicated(cr.ids))) {
    warning("Non-unique year/interval pairs in input data")
    return(df)
  }
  rows <- match(df.ids, cr.ids)
  tmp  <- !is.na(rows)
  rows <- rows[tmp]
  if (!length(rows)) {
    warning("No matching rows between input and results data")
    return(df)
  }
  for (v in add) {
    df[, v]    <- NA
    df[tmp, v] <- crntdata[rows, v, drop=TRUE]
  }
  df
}

