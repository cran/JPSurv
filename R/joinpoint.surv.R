### Below are the R code changes by FZ
### - aapc.multiints was added for multiple selected intervals. 08/15/2019
### - .getSeAapc and aapc.multiints were updated for user-specified year range. 09/09/2019
### - Functions plot.surv() and .comment() was removed. 04/01/2020
### - Functions called in deleteLastInterval were copied from CoxFit.R. 04/06/2020
### - The warning/error message was added if no data available for the selected cohort. 04/06/2020
### - aapc.multiints was modified to return ACS trend for user-specified range only. 05/19/2020
################################################################################################
# .CoxFit() is the core function to fit the proportional hazard relative survival models.
# Inputs:
#     X1:     the design matrix; 	X1 = joinpointMatrix + model.matrix(~-1+Year+Interval_);
#             para: delta, beta_year, alpha
#     nAlive: number of patients at risk
#     nDied:  number of people died
#     nLost:  number of peple lost to follow-up
#     expSurv:expected survival
#     Z: design matrix for covariates
#        para: gamma

.CoxFit = function(X1, nAlive, nDied, nLost, expSurv,Z=NULL) {
	stopifnot(sum(expSurv <= 0) == 0);
	N = nAlive - nLost / 2;
	Y = N - nDied;
	p = dim(X1)[2];
	numparaCovar=0;
	if(!is.null(Z)) numparaCovar = dim(Z)[2];

	min_value = 1e-6;
	min_value2 = 1e-15;
	
	# RelativeSurvivalLinkFun is the link function for the relative survival model.
	RelativeSurvivalLinkFun = list();
	RelativeSurvivalLinkFun$LogLik = function(beta1,gamma1) {
		xbeta1 = X1 %*% beta1;
		zgamma1 = 0;
		if(numparaCovar>0) zgamma1 = Z %*% gamma1;
		prob = exp(-exp(xbeta1+zgamma1)) * expSurv; # is p_j(x)*E_j(x)
		nlogl = sum(Y * (log(prob + min_value2) - log((Y+min_value2)/N)) 
				+ (N - Y) * (log(1 - prob + min_value2) - log((N-Y+min_value2)/N))); # including the calculation of BIC
		return(nlogl);
	}
	
	RelativeSurvivalLinkFun$UW = function(beta,gamma) {
		xbeta = X1 %*% beta;
		zgamma = 0;
		if(numparaCovar>0) zgamma = Z %*% gamma;
		prob = exp(-exp(xbeta+zgamma)) * expSurv;
		index = (prob > expSurv - min_value);
		prob[index] = (expSurv - min_value)[index];
		prob[prob < min_value] = min_value;
		grad_ = 1 / prob / log(prob / expSurv);
		U = (Y / N - prob) * grad_;
		W = prob * (1 - prob) / N * grad_ * grad_;
		
		result = list();
		result$U = U;
		result$W = W;
		return(result);
	}

	RelativeSurvivalLinkFun$f = function(xbeta,zgamma) {
		result = exp(-exp(xbeta+zgamma));
		return(result);
	}
	
	# GetEstimate() is the function to compute the estimates using the iteratively reweighted least square algorithm.
	GetEstimate = function(linkfun) {
		old_beta = rep(0, p);
		old_gamma = numeric(0);
		if(numparaCovar>0) old_gamma = rep(0, numparaCovar);
		old_ll = linkfun$LogLik(old_beta,old_gamma);
		converged = FALSE;
		lm_fit = NULL;
		fitNACols = rep(FALSE,p+numparaCovar);
		for (i in 1:100) {
			UW = linkfun$UW(old_beta,old_gamma);
			X = X1; #local variable
			if(!is.null(Z)){
			 lm_fit = lm(UW$U ~ -1 + X + Z, weights = 1 / UW$W);
			}else{
			 lm_fit = lm(UW$U ~ -1 + X, weights = 1 / UW$W);
			}

			deltagamma = lm_fit$coefficients;
			NACols = is.na(deltagamma);
			deltagamma[NACols] = 0;
			scale = 1.0;
			new_beta = NA;
			new_gamma = NA;
			new_ll = NA;
			for (j in 1:20) {
				new_beta = old_beta + deltagamma[1:p] * scale;
				new_gamma = old_gamma;
				if(numparaCovar>0) new_gamma = old_gamma + deltagamma[p+1:numparaCovar] * scale;
				new_ll = linkfun$LogLik(new_beta,new_gamma);
				if (new_ll > old_ll) break
				scale = scale / 2;
			}
			if (new_ll < old_ll + 1e-4) {
				if (new_ll > old_ll) {
					converged = TRUE;
				}
				fitNACols = NACols;
				break;
			}
			old_beta = new_beta;
			old_gamma = new_gamma;
			old_ll = new_ll;
		}

		if (i == 100) converged = FALSE;
		
		result = list();
		result$X1names=colnames(X1);
		result$X1=X1; ### FZ 04/06/2020
		result$covarnames=character(0);
		if(!is.null(Z)) {
      result$Z0=Z[1,,drop=FALSE];#t(apply(Z,2,min));
  		result$covarnames=colnames(Z);
    }
    varmatrixnames = c(result$X1names,result$covarnames);
		#result$beta = new_beta;
		result$xbeta = X1 %*% new_beta;
		if(!is.null(Z)) result$zgamma = Z %*% new_gamma;
		#pred = linkfun$f(result$xbeta);
		#result$pred = pred;
		result$ll = new_ll;
		Estimates = c(new_beta,new_gamma);
		summ_fit = summary(lm_fit);
		Std.Error = summ_fit$coefficients[, "Std. Error"];
		Std.Error[fitNACols]=0;
		result$coefficients = cbind(Estimates, Std.Error);
		.fillNA=function(mdata,NAs,varmatrixnames){
		  NAs = as.vector(NAs);
		  res = mdata;
		  if(length(NAs)>0)for(Index in 1:length(NAs)){
		    if(NAs[Index]){
		      if(Index == 1){
		        if(nrow(res)>=1){
		          res = rbind(matrix(0,1,ncol(res)),res[Index:nrow(res),]);
		        }else{
		          res = matrix(0,1,ncol(res));
		        }
		      }else if (Index == length(NAs) && Index >1){
		        res = rbind(res[1:(Index-1),],matrix(0,1,ncol(res)));
          }else{
            if(nrow(res)>=Index){
		          res = rbind(res[1:(Index-1),],matrix(0,1,ncol(res)),res[Index:nrow(res),]);
		        }else{
		          res = rbind(res[1:(Index-1),],matrix(0,1,ncol(res)));
		        }
          }
		    }
		  }
		  if(length(NAs)>0)for(Index in 1:length(NAs)){
		    if(NAs[Index]){
		      if(Index == 1){
		        if(ncol(res)>=1){
		          res = cbind(matrix(0, nrow(res),1),res[,Index:ncol(res)]);
		        }else{
		          res = matrix(0, nrow(res),1);
		        }
		      }else if (Index == length(NAs) && Index >1){
		        res = cbind(res[,1:(Index-1)],matrix(0, nrow(res),1));
          }else{
            if(ncol(res)>=Index){
		          res = cbind(res[,1:(Index-1)],matrix(0, nrow(res),1),res[,Index:ncol(res)]);
		        }else{
		          res = cbind(res[,1:(Index-1)],matrix(0, nrow(res),1));
		        }
          }
		    }
		  }
		  colnames(res)[NAs] = varmatrixnames[NAs];
		  rownames(res) = colnames(res);
		  return(res);
		}
		result$covariance = .fillNA(vcov(lm_fit),fitNACols,varmatrixnames);
		result$aic = 2 * (p+numparaCovar) - 2 * new_ll;
		result$bic = (p+numparaCovar) * log(sum(N)) -2 * new_ll;
		result$converged = converged;
		result$p=p;
		result$numparaCovar=numparaCovar;
		
		return(result);
	}

	fit = GetEstimate(RelativeSurvivalLinkFun);

	return(fit);
}
.getSeAapc = function(fit, fnAapc,Z0,ACS.range) {
  if(!is.null(Z0)) Z0 = t(data.frame(as.vector(Z0)));
	beta = fit$coefficient[1:fit$p, 1];
	gamma = NULL;
	if(fit$numparaCovar>0) gamma = fit$coefficient[fit$p+1:fit$numparaCovar, 1];

	deltabeta = fit$coefficient[1:fit$p, 2] * 1e-4;
	deltagamma = numeric(0);
	if(fit$numparaCovar>0) deltagamma = fit$coefficient[fit$p+1:fit$numparaCovar, 2] * 1e-4;

	est = fnAapc(fit,beta,gamma,Z0);

	jacob = matrix(NA, nrow(fit$coefficient), length(est));
	for (i in 1:length(beta)) {
		beta1 = beta;
		beta1[i] = beta[i] + deltabeta[i];
		est1 = fnAapc(fit, beta1,gamma,Z0);
		jacob[i, ] = (est1 - est) / deltabeta[i];
	}
	if(fit$numparaCovar>0)for (i in 1:length(gamma)) {
		gamma1 = gamma;
		gamma1[i] = gamma[i] + deltagamma[i];
		est1 = fnAapc(fit, beta,gamma1,Z0);
		jacob[fit$p+i, ] = (est1 - est) / deltagamma[i];
	}
	varAapc = t(jacob) %*% fit$covariance %*% jacob;
	if (length(est) == 1) seAapc = sqrt(varAapc)
	else seAapc = sqrt(diag(varAapc));
	
	if(!is.null(ACS.range)){
	  fit$apc<-fit$apc[1,]
	  fit$apc$start.year<-ACS.range[1]
	  fit$apc$end.year<-ACS.range[2]
	}
	result = fit$apc;
	result$estimate = est;
	result$std.error = seAapc;
	result$lowCI = result$estimate + (qnorm(0.025,0,1)*result$std.error);
	result$upCI = result$estimate + (qnorm(0.975,0,1)*result$std.error);
	result = data.frame(result);
	rownames(result) = NULL;
	return(result);
}

.getSegX=function(jpoints,yearvalue){
 nJP=length(jpoints);
 lineSeg=numeric(0);
	if (nJP > 0) {
		for (i in 1:nJP) {
			seg = yearvalue - jpoints[i];
			seg[seg < 0] = 0;
			lineSeg = c(lineSeg, seg);
		}
	}
	return(lineSeg);
}
.getIntervalX=function(numIntervals,intervalvalue){
 IntervalX=numeric(numIntervals);
 IntervalX[intervalvalue]=1;
	return(IntervalX);
}
.calcDerivPredInt <- function(cox_fit,jpoints,betagamma,yearvalue,intervalvalue,z0){
 #res=numeric(cox_fit$p+cox_fit$numparaCovar)
 #jpoints=cox_fit$jp

 numIntervals=cox_fit$p-length(jpoints)-1
 SegX=.getSegX(jpoints,yearvalue)

 xX = yearvalue
 IntervalX=.getIntervalX(numIntervals,intervalvalue)
 res=c(SegX,xX,IntervalX,z0)

 logsurv=-exp(t(betagamma)%*%res)
 if(length(logsurv) == 1){     ## FZ 07/25/2019
   logsurv <- rep(logsurv, length(res))
 }   
 res = res * logsurv * exp(logsurv)
 return(res)
}

.calcDerivPredCum=function(cox_fit,jpoints,betagamma,yearvalue,intervalvalue,z0){
 res=numeric(cox_fit$p+cox_fit$numparaCovar);
 #jpoints=cox_fit$jp;
 numIntervals=cox_fit$p-length(jpoints)-1;
 SegX=.getSegX(jpoints,yearvalue);
 xX = yearvalue;
 logsurv=0;
 if(intervalvalue>0) for(IntIndex in 1:intervalvalue){
   IntervalXi=.getIntervalX(numIntervals,IntIndex);
   X1Zi=c(SegX,xX,IntervalXi,z0);
   logsurvi=-exp(t(betagamma)%*%X1Zi);   
   if(length(logsurvi) == 1){     ## FZ 07/25/2019
     logsurvi <- rep(logsurvi, length(res))
   }   
   logsurv=logsurv+logsurvi;
   res = res + X1Zi* logsurvi;
 }
 res = res * exp(logsurv);
 return(res);
}
.is.substring = function(x, string,ignore.case = TRUE,...){
  x=as.character(x);
  string=as.character(string);
#  return(length(grep(x, string,ignore.case = TRUE,...)) == 1);
  return(length(grep(x, string,ignore.case=ignore.case,...)) == 1);
}
.is.substring.vector = function(x, strvector,ignore.case = TRUE,...){
  x=as.character(x);
  strvector=as.character(strvector);
  res=logical();
  if(length(strvector)>0)
  for(Index in 1:length(strvector)){
    res[Index]=.is.substring(x,strvector[Index],ignore.case=ignore.case,...)
  }
#  return(length(grep(x, string,ignore.case = TRUE,...)) == 1);
  return(res);
}
.search.substr=function(value, valuelist,ignore.case = FALSE,...){
  res=NULL;
  if(length(valuelist)>0){
	  substrtrue = .is.substring.vector(value,valuelist,ignore.case = ignore.case,...);
	  res = (1:length(valuelist))[substrtrue];
  }
  return(res);
}
.covert.numeric=function(data,cols){
  if(length(cols)>0)for(Index in 1:length(cols)){
    data[,cols[Index]]=as.numeric(data[,cols[Index]]);
  }
  return(data);
}
.getLevels.keeporder=function(x){
  x=as.vector(x);
  res=NULL;
  if(length(x)>0){
    # BW 12/03/2020 Added to remove warning
    xv    <- x
    index <- 1:length(x)
    temp=data.frame(xv=x,index=index);
    templevels=as.character(levels(factor(x)));
    Tempdata=NULL;
    for(Index in 1:length(templevels)){
      crntData=subset(temp,xv==templevels[Index]);
      minrowindex=as.character(min(as.numeric(crntData[,"index"])));
      crntrow=subset(crntData,index==minrowindex);
      Tempdata=rbind(Tempdata,crntrow);
    }
    Tempdata=Tempdata[order(Tempdata[,"index"]),,drop=FALSE];
    res=as.character(Tempdata[,"xv",drop=TRUE]);
    #print(Tempdata);
  }
  return(res);
}
.isequal.vector=function(x,y){
  res=FALSE;
  x=as.vector(x);
  y=as.vector(y);
  if(length(x) == length(y)){
    truevalues = (x == y);
    if(length(x) == length(truevalues[truevalues==TRUE])){
      res=TRUE;
    }
  }
  return(res);
}
.getLevels.keeporder.data.frame=function(x){
  x=as.data.frame(x);
  #X: column: by variables
  #   row: records, i.e. levels
  res=x[0,];
  if(nrow(x)>0 && ncol(x)){
    res=rbind(res,x[1,]);
    if(nrow(x)>1)for(RowIndex in 2:nrow(x)){
      NewItem=TRUE;
      crntrow=as.vector(t(x[RowIndex,,drop=FALSE]));
      for(Index in 1:nrow(res)){
        crntItem=as.vector(t(res[Index,,drop=FALSE]));
        if(.isequal.vector(crntrow,crntItem)){
          NewItem=FALSE;
          break;
        }
      }
      if(NewItem){
        res=rbind(res,crntrow);
      }
    }
  }
  return(res);
}
.getsubset.byvar=function(byvarnames, byvarvalues,mydata,truefalse=TRUE){
  res=data.frame();
  if(dim(mydata)[1]*dim(mydata)[2]>0){
    res = mydata;
    if(length(byvarnames)>0){
      truefalse=rep(truefalse,length(byvarnames));
      for(Index in 1:length(byvarnames)){
        rowindicators=res[,byvarnames[Index],drop=TRUE]==byvarvalues[Index];
        if(!truefalse[Index]){
          rowindicators = (!(rowindicators));
        }
        res = res[rowindicators,,drop=FALSE];
      }
    }
  }
  return(res);
}
.getRowAsVector=function(x,RowIndex){
  return(as.vector(t(x[RowIndex,,drop=FALSE])));
}
.deleteLastInterval.onecorhot = function(yearcol,intervalcol,mydata){
  res=NULL;
  if(dim(mydata)[1]*dim(mydata)[2]>0){
    YearRange=as.numeric(.getLevels.keeporder(mydata[,yearcol,drop=TRUE]));
    YearRange=YearRange[order(YearRange)];
    for(YearValue in YearRange){
      cntYearData=.getsubset.byvar(byvarnames=yearcol,byvarvalues=YearValue,mydata);
      #get interval values
      IntervalRange=as.numeric(.getLevels.keeporder(cntYearData[,intervalcol,drop=TRUE]));
      IntervalRange=IntervalRange[order(IntervalRange)];
      
      DelIntervalValue=IntervalRange[length(IntervalRange)];
      
      #delete the rows associated with the biggest interval values
      res=rbind(res,.getsubset.byvar(byvarnames=intervalcol,byvarvalues=DelIntervalValue,cntYearData,truefalse=FALSE));
    }
  }
  return(res);
}



#############################
# exported functions
#############################

aapc_OLD = function(fit, type="AbsChgSur", interval=5) {
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
			pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
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
			pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
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
			pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			
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
	    pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
	    pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
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
			pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
			pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
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
}
aapc.multiints_OLD<-function(fit, type="AbsChgSur",int.select=NULL,ACS.range=NULL,ACS.out=NULL){
  Z0<-NULL
  #calculate for each segment: AAPC measure estimate and SE
  #output columns: start.year end.year     estimate    std.error PredInterval
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
      pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);

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
      pred1 = fit$Predict(years + epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred2 = fit$Predict(years - epsilon, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
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
      pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);      
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
      pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma)
      pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma)

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
      pred1 = fit$Predict(years , fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
      pred2 = fit$Predict(years + 1, fup, beta_input=beta,Z0=Z0, gamma_input=gamma);
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
      aapc.results<-aapc.results.range
    }
  }
  return(aapc.results)
}

#delete the records of last intervals of all years:
#input:
#   byvarnames: the names of by variables. It could be a list of zero length, i.e. c()
#   yearcol,intervalcol: column names of year variable and interval variable respectively
#   mydata: a data.frame that contains data
deleteLastInterval = function(data,byvarnames,yearcol,intervalcol){
  byVarData = data[,byvarnames,drop=FALSE];
  res = NULL;
  if(dim(byVarData)[1]*dim(byVarData)[2]>0){
    byvarlevels = .getLevels.keeporder.data.frame(byVarData);
    for(Index in 1:nrow(byvarlevels)){
      crntlevels = .getRowAsVector(byvarlevels,Index);
      crntdata = .getsubset.byvar(byvarnames,crntlevels,data);
      res = rbind(res, .deleteLastInterval.onecorhot(yearcol,intervalcol,crntdata));
    }
  }else{
    res=.deleteLastInterval.onecorhot(yearcol,intervalcol,data);
  }
  return(res);
}

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
joinpoint_OLD <- function(data, subset=NULL, na.action = na.fail, year="Year", interval="Interval", 
                      number.event="Died", number.alive="Alive_at_Start",number.loss="Lost_to_Followup", 
                      expected.rate="Expected_Survival_Interval", observedrelsurv = NULL, model.form = NULL,
                      maxnum.jp = 0, proj.year.num = 5, # Edited 04/07/2020 FZ;Added 8/24/2015 DM - Number of projection years, range 0-30, default 5 set below in code
                      op=list(), delLastIntvl=FALSE){
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

#codes of this function begins --------------
#
# inner function section: define inner function here
#
	# .indv_bic() computes the BIC for the given join points.
	.indv_bic = function(jpoints,Year,nAlive, nDied, nLost, ExpSurv,X,Z,Interval) {

		lineSeg = NULL;
		nJP = length(jpoints);
		if (nJP > 0) {
			for (i in 1:nJP) {
				seg = Year - jpoints[i];
				seg[seg < 0] = 0;
				lineSeg = cbind(lineSeg, seg);
			}
			colnames(lineSeg) = paste("jp", 1:nJP, sep = "_");
		}
		
		X1 = cbind(lineSeg, X)
		cox_fit = .CoxFit(X1, nAlive, nDied, nLost, ExpSurv,Z);
		
		# run checks on proj.year.num;
		if (is.null(proj.year.num)) { 
			proj.year.num = 5;
		}else{
			if (proj.year.num<0 | proj.year.num>30) proj.year.num = 5;
		}
		
		# predict() computes the predicted cumulative survival rates.
		predict = function(years = NULL, intervals = NULL, beta_input = NULL,Z0=NULL, gamma_input = NULL) {
			if (is.null(years)) years = min(Year):(max(Year) + proj.year.num)
			if (is.null(intervals)) intervals = (1:max(Interval))
			if (is.null(Z0)) Z0 = cox_fit$Z0;
			if (is.null(beta_input)) beta_input = cox_fit$coefficients[1:cox_fit$p, 1];
			if(cox_fit$numparaCovar>0){
  			if (is.null(gamma_input)) gamma_input = cox_fit$coefficients[cox_fit$p+1:cox_fit$numparaCovar, 1];
			} else{
  			  gamma_input=NULL;
			}
			maxInterval = max(intervals);
			stopifnot(maxInterval <= max(Interval));
			
			beta2 = list();
			#cox_fit$beta = cox_fit$coefficients[, 1];
			#intercept = cox_fit$beta[nJP + 1];
			intercept = 0;
			if (nJP > 0) beta2$Seg = beta_input[1:nJP]
			beta2$Year = beta_input[nJP + 1];
			#beta2$Interval = c(0, beta_input[(nJP+2):length(beta_input)]);
			#beta2$Interval = beta2$Interval + intercept;
			beta2$Interval = beta_input[(nJP+2):cox_fit$p];
			if(cox_fit$numparaCovar>0){
                       beta2$betagamma = c(beta_input,gamma_input);
                     } else{
                       beta2$betagamma = c(beta_input);
                     }
			
			nrowcovars=1;
			if (!is.null(Z0)) {
                       nrowcovars = nrow(Z0)
                     } else {
                       Z0=0;
                     };
			result = data.frame(matrix(0, length(years) * maxInterval * nrowcovars, 4+2+length(cox_fit$covarnames)));
			names(result) = c("Year", "Interval", "pred_int", "pred_cum","pred_int_se","pred_cum_se",cox_fit$covarnames);

			for(ZIndex in 1:nrowcovars){
                       surv_int = matrix(NA, length(years), maxInterval);
  			  surv_cum = surv_int;
  			  surv_xbeta = surv_int;
			  z0=numeric(0);
  			  zgamma=0;
  		 	  if (!is.null(gamma_input) && cox_fit$numparaCovar>0) {
                         zgamma=Z0[ZIndex,,drop=FALSE] %*% gamma_input;
                         z0=Z0[ZIndex,,drop=FALSE];
                       }
        
                       for (i in 1:length(years)) {
  				for (j in 1:maxInterval) {
  					xbeta2 = beta2$Interval[j] + years[i] * beta2$Year;
  					if (nJP > 0) for (k in 1:nJP) {
  						xseg = years[i] - jpoints[k];
  						if (xseg < 0) xseg = 0
  						xbeta2 = xbeta2 + beta2$Seg[k] * xseg;
  					}
  					surv_xbeta[i, j] = xbeta2+zgamma;
  					surv_int[i, j] = exp(-exp(xbeta2+zgamma));
  					if (j == 1) surv_cum[i, j] = surv_int[i, j]
  					else surv_cum[i, j] = surv_cum[i, j-1] * surv_int[i, j]
  					idx = (ZIndex-1)*length(years)*maxInterval + (i - 1) * maxInterval + j;
  					result[idx, 1] = years[i];
  					result[idx, 2] = j;
  					result[idx, 3] = surv_int[i, j];
  					result[idx, 4] = surv_cum[i, j];

  					DerivPredInt = .calcDerivPredInt(cox_fit,jpoints,beta2$betagamma,yearvalue=years[i],intervalvalue=j,z0=z0);
  					result[idx, 5] = sqrt(t(DerivPredInt)%*% cox_fit$covariance %*%DerivPredInt);
  					DerivPredCum = .calcDerivPredCum(cox_fit,jpoints,beta2$betagamma,yearvalue=years[i],intervalvalue=j,z0=z0);
  					result[idx, 6] = sqrt(t(DerivPredCum)%*% cox_fit$covariance %*%DerivPredCum);
  					if(cox_fit$numparaCovar>0){
  					  result[idx, 7:ncol(result)] = Z0[ZIndex,,drop=TRUE]; 
  					}
  				}
  			  }
			}
			#result1 = subset(result, Interval %in% intervals);
			requested = match(result$Interval, intervals,0L);
			result1 = result[requested>0,];
			return(result1);

              } # END: predict
 		
		getApc = function() {
			result = matrix(NA, nJP+1, 3);
			result = data.frame(result);
			dimnames(result)[[2]] = c("start.year", "end.year", "estimate");
			#dimnames(result)[[1]] = paste("Seg", 1:(nJP+1), sep = " ");
			rownames(result) = NULL;
			result$estimate[1] = cox_fit$coefficients[nJP+1, 1];
			result$start.year[1] = min(Year);
			result$end.year[nJP + 1] = max(Year);
			if (nJP > 0) {
				for (i in 1:nJP) {
					result$estimate[i + 1] = result$estimate[i] + cox_fit$coefficients[i, 1];
					result$end.year[i] = jpoints[i];
					result$start.year[i + 1] = jpoints[i];
				}
			}
			return(result);
		}
	
		cox_fit$Predict = predict;
		cox_fit$apc = getApc();
		return(cox_fit);
	}

  .handle.op=function(x,default=3){
    if(is.null(x)){
      res=default;
    }else{
      res=x;
    }
    return(res);  
  }
  .get.integerlist=function(n,multiplier){
    if(n>0){
      res=(1:n)*multiplier;
    }else{
      res=integer(0);
    }
    return(res);
  }
	.getJP = function(nJP,Year,...) {
    if(nJP>0){
  #op=list(), #options 
  #             numbetwn: number of skipped obs between joinpoints exclusive (not count for the joinpoints);
  #             numfromstart: number of skipped obs from the first obs to joinpoints exclusive (not count for the joinpoint);
  #             numtoend: number of skipped obs from the first obs to joinpoints exclusive (not count for the joinpoint);
    op$numbetwn=.handle.op(op$numbetwn,2);
    op$numfromstart=.handle.op(op$numfromstart,3);
    op$numtoend=.handle.op(op$numtoend,5);

		intervalSize = op$numbetwn+1;
		#what implemented before:  so called intervalSize fixed as 3
    #   initial joinpoints: numfromstart=3, numbetwn=2, numtoend=5
		#   find next: numbetwn=2, numtoend=3
		#   exit: numtoend=5
		stopifnot(max(Year) >= min(Year) + op$numfromstart + (nJP-1) * intervalSize + op$numtoend);# here it is consistent with what's claimed for op 
		jpoints = min(Year) + op$numfromstart + c(0,.get.integerlist(nJP-1, intervalSize));# here it is consistent with what's claimed for op 
		endFlag = FALSE;
		#numtoend.temp=op$numtoend-1;#It is not consistent with what's claimed for op, in future replace op$numtoend for numtoend.temp below
		numtoend.temp=op$numtoend;# here it is consistent with what's claimed for op 
		nextjp = function(oldjp) {
			oldjp[nJP] = oldjp[nJP] + 1;
			if (nJP >= 2) for (k in nJP:2) {
				if (oldjp[k] > max(Year) - intervalSize * (nJP - k) - numtoend.temp) {
					oldjp[k - 1] = oldjp[k - 1] + 1;
					oldjp[k:nJP] = oldjp[k - 1] + intervalSize * (1:(nJP - k + 1));
				}
			}
			if (oldjp[1] > max(Year) - intervalSize * (nJP-1)-op$numtoend) endFlag <<- TRUE;# here it is consistent with what's claimed for op 
			return(oldjp);
		}
		result = list(); 
		result$jp = NA;
		result$bic = Inf;
		nIter = 0;
		jp.debug = FALSE;  # setting the debug option
		while (!endFlag){
			new_bic = .indv_bic(jpoints,Year=Year,...);
			if (new_bic$bic < result$bic) {
				result = new_bic;
				result$jp = jpoints;
			}
			jpoints = nextjp(jpoints);
			nIter = nIter + 1;
			if (jp.debug & nIter > 3) break
		}
		}else{
  		result = .indv_bic(NULL,Year=Year,...);
  		result$jp = NULL;
		}

		return(result);
	}


	.getBestJP = function(nJP,...) {
              FitList=list();
		bestFit = .getJP(0,...);
		FitList[[length(FitList)+1]]=bestFit;
		if (nJP > 0) {
			for (k in 1:nJP) {
				message(paste0("Computing estimates when number of join points is ", k))
				newFit = .getJP(k,...);
				FitList[[length(FitList)+1]]=newFit;
				if (newFit$bic < bestFit$bic) bestFit = newFit
			}
		}
		return(list(bestFit=bestFit, FitList=FitList));
	}

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

  fitres = .getBestJP(maxnum.jp,Year,nAlive, nDied, nLost, ExpSurv,X,Z,Interval);

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

  pred = cox_fit$Predict(Z0=iZ0);
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
      pred = FitList[[Index+1]]$Predict(Z0=iZ0);
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
  class(cox_fit) <- "joinpoint";
  return(cox_fit);

}



