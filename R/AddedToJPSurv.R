#########################################################################################################
#### This R file includes more functions for JPSurv tool.
#########################################################################################################
#### Functions created by Fanni Zhang.
####
#### Functions were created as requested for JPSurv web app new features (discussed in meeting 06/05/2019)
#### -download.data
#### -plot.dying.year.annotate
#### -plot.surv.year.annotate
#### -plot.surv.int.multiyears
####
#### Function input.valid was added. (meeting 03/31/2020)
#### Dec 08 2020: Fix R CMD check warnings with obsNA, trend.label in plot functions
#########################################################################################################


#########################################################################################################
# input.valid: a function that validates the selected cohort data and returns an indicator 
#              along with a warning message if the selected cohort is not valid.
# - valid=1 if the cohort selection is available; otherwise, valid=0.
# arguments:
# input - input dataset read in by function joinpoint.seerdata.
# subset - an optional string specifying a subset of observations used in the fitting process.
#########################################################################################################

input.valid<-function(input,subset=NULL){
  if(is.null(subset)){
    input.sub<-input
  }else{
    input.sub<-subset(input,eval(parse(text=subset)))
  }
  valid<-1
  if(dim(input.sub)[1]==0){
    error.str<-paste("No data available for the cohort selection '",subset,"' in the input data file.",sep="")
    warning(error.str)
    valid<-0
  }
  return(valid)
}

#########################################################################################################
# download.data: a function that returns a merged data set for download including 
# - the selected cohort input data
# - the accompanying data for plot.surv.year/plot.dying.year
# arguments:
# input - input dataset read in by function joinpoint.seerdata.
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# type - two options for this argument: "graph" for graph data and "full" for full data.
# subset - an optional string specifying a subset of observations used in the fitting process.
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# int.select - the interval values selected for the plot if the downloadtype="graph". The default is NULL.
#########################################################################################################

download.data<-function(input,fit,nJP,yearvar,downloadtype,subset=NULL,interval="Interval",int.select=NULL){
  if(is.null(subset)){
    input.sub<-input
  }else{
    input.sub<-subset(input,eval(parse(text=subset)))
  }
  if(downloadtype=="graph"){
    output.sub<-fit$FitList[[nJP+1]]$predicted
    if(!is.null(int.select)){
      ints.pred<-unique(output.sub[,interval])
      if(length(which(!int.select %in% ints.pred))>0){
        intstr<- paste(int.select,collapse=",",sep="")
        intstr.na<-paste(int.select[which(!int.select %in% ints.pred)],collapse=",",sep="")
        warn.str<-paste("int.select=c(",intstr,"). Data not available for interval(s) ",
                        intstr.na," in the selected cohort.",sep="")
        warning(warn.str)
      }
      output.sub<-output.sub[which(output.sub[,interval] %in% int.select),]
    }
  }else if(downloadtype=="full"){
    output.sub<-fit$FitList[[nJP+1]]$fullpredicted
  }else{
    stop("There is no such type of data for download. Please define the type using either 'graph' or 'full'.")
  }
  yearint.input<-paste(input.sub[,yearvar],input.sub[,interval],sep="_")
  yearint.output<-paste(output.sub[,yearvar],output.sub[,interval],sep="_")
  rows.match<-match(yearint.output,yearint.input)
  
  cols.input<-colnames(input.sub)
  if("Page_type" %in% cols.input){
    cols.input<-cols.input[-which(cols.input=="Page_type")] 
  }
  if("Expected_Survival_Interval" %in% cols.input){
    cols.input<-cols.input[-which(cols.input=="Expected_Survival_Interval")] 
  }
  cols.keep<-cols.input
  del.cols<-c("Observed_SE_Interval","Observed_SE_Cum")
  hold.cols<-c(yearvar,interval)
  if(length(which(!del.cols %in% cols.input))==0){
    cols.keep<-cols.input[-which(cols.input %in% del.cols)] 
  }
  if(length(which(!hold.cols %in% cols.keep))==0){
    cols.input<-cols.keep[-which(cols.keep %in% hold.cols)] 
  }
  merge.sub<-output.sub
  merge.sub[,cols.input]<-input.sub[rows.match,cols.input]
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_int")]<-"Predicted_Survival_Int"
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_cum")]<-"Predicted_Survival_Cum"
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_int_se")]<-"Predicted_Survival_Int_SE"
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_cum_se")]<-"Predicted_Survival_Cum_SE"
  merge.sub[,"Predicted_ProbDeath_Int"]<-1-merge.sub$Predicted_Survival_Int
  merge.sub[,"Predicted_ProbDeath_Int_SE"]<-merge.sub$Predicted_Survival_Int_SE
  cols.pred<-c("Predicted_Survival_Int","Predicted_ProbDeath_Int","Predicted_Survival_Cum",
               "Predicted_Survival_Int_SE","Predicted_Survival_Cum_SE","Predicted_ProbDeath_Int_SE")
  merge.data<-merge.sub[,c(cols.keep,cols.pred)]
  if(!is.null(subset)){
    subset.list<-strsplit(subset, " & ")
    subset.list[[1]]<-subset.list[[1]][which(grepl("==", subset.list[[1]])==TRUE)]
    cohort.list<-strsplit(subset.list[[1]], "==")
    cohort.vars<-sapply(cohort.list, "[", 1)
    cohort.vars<-gsub(" ","",cohort.vars)
    cohort.values<-sapply(cohort.list, "[", 2)
    cohort.values<-gsub(" ","",cohort.values)
    cohort.vars<-cohort.vars[which(!cohort.vars %in% c(yearvar,interval))]
    cohort.values<-cohort.values[which(!cohort.vars %in% c(yearvar,interval))]
    if(length(cohort.vars)==1){
      merge.data[,cohort.vars]<-replicate(dim(merge.data)[1],cohort.values)
    }else if(length(cohort.vars)>1){
      merge.data[,cohort.vars]<-t(replicate(dim(merge.data)[1],cohort.values))
    }
  }
  
  merge.data[,yearvar]<-as.numeric(merge.data[,yearvar])
  return(merge.data)
}

#########################################################################################################
# plot.dying.year.annotate: a function that returns a plot for Percent Change in the Annual Probability of
# Dying of Cancer by Diagnosis year using ggplot
# The annotation feature is available for nJP<=3 and the number of multiple intervals selected <=3
# arguments:
# plotdata - the graph data returned by function download.data with downloadtype="graph".
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obsintvar - the variable name for observed interval survival. The default is "Relative_Survival_Interval"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predintvar - the variable name for predicted interval survival. The default is "Predicted_ProbDeath_Int".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# annotation - the indicator for the annotation feature. The default is 0 (no annotation on the plot).Two plots
#              with and without annotation will be returned in a list when annotation=1.
# topanno - the indicator for showing the top curve annotation. The default is 1 (annotation for the top curve).
# trend - the indicator for returning the trend measure tables. The default is 0 (no trend tables returned).
#########################################################################################################

### define a plot function for Percent Change in the Annual Probability of Dying of Cancer by Diagnosis year with annotations
Plot.dying.year.annotate<-function(plotdata,fit,nJP,yearvar,obsintvar="Relative_Survival_Interval",predintvar="Predicted_ProbDeath_Int",interval="Interval",annotation=0,topanno=1,trend=0,title=NULL){
  if(is.null(title)){
    title.rch<-"Annual Probability of Dying of Cancer by Diagnosis Year"
  }else{
    title.rch<-title
  }
  interval.values<-as.numeric(unique(plotdata[,"Interval"]))
  interval.labels<-paste((interval.values-1)," to ",interval.values," years since diagnosis",sep="")
  if(interval.values[1]==1){
    interval.labels[1]<-"0 to 1 year since diagnosis"
  }
  trends.out<-list()
  rows.obs<-which(is.na(plotdata[,obsintvar]))
  plotdata[,"obsNA"]<-0
  plotdata[rows.obs,"obsNA"]<-1
  plotdata[,"obsintvar_1"]<-1-plotdata[,obsintvar]
  obsintvar_1minus<-"obsintvar_1"
  if(trend==1 | annotation==1){ 
    trends.out <- aapc.multiints(fit$FitList[[nJP+1]], type="RelChgHaz", int.select=interval.values)
    if(length(interval.values)<=3 & nJP<=3){
      ### haz results are the same for all interval values
      jp.loc<-fit$FitList[[nJP+1]]$jp
      annot.strs<-paste(sprintf("%.1f",100*trends.out[[1]]$estimate),"%",sep="")
      x.values<-list()
      plotdata[,"trend.label"]<-""
      for(i in 1:length(interval.values)){
        int.i<-interval.values[i]
        haz.apc<-trends.out[[i]]
        plotdata.i<-plotdata[which(plotdata[,interval]==int.i),]
        if(max(plotdata.i[,yearvar])==max(haz.apc$end.year) | nJP==0){
          end.year<-haz.apc$end.year
          end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
          x.values.i<-1/2*(haz.apc$start.year+end.year)
        }
        if(max(plotdata.i[,yearvar])<max(haz.apc$end.year)){
          if(nJP==1){
            if(max(plotdata.i[,yearvar])>jp.loc){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }else{
              haz.apc<-haz.apc[1:nJP,]
              x.values.i<-1/2*(haz.apc$start.year+haz.apc$end.year)
            }
          } ## nJP=1 end
          if(nJP==2){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>min(jp.loc)){
              haz.apc<-haz.apc[1:nJP,]
              x.values.i<-1/2*(haz.apc$start.year+haz.apc$end.year)
            } else{
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=2 end
          if(nJP==3){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              haz.apc<-haz.apc[1:nJP,]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else {
              haz.apc<-haz.apc[1:(nJP-2),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=3 end
        } 
        x.values[[i]]<-round(x.values.i)
        plotdata$trend.label[which((plotdata[,yearvar] %in% x.values[[i]]) & plotdata[,interval]==interval.values[[i]])]<-annot.strs[1:length(x.values[[i]])]  ### 11/25
      }## interval.values iteration end
      if(annotation==1){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        if(topanno==1){
          ymeans.byint<-aggregate(plotdata[,"Predicted_ProbDeath_Int"], list(plotdata[,interval]), mean)
          topint<-which(ymeans.byint[,2]==max(ymeans.byint[,2]))
          plotdata$trend.label[which(plotdata[,interval]!=interval.values[topint])]<-""
          plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
            geom_point(data=(subset(plotdata, plotdata[, "obsNA"]==0)), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obsintvar_1minus)), 
                                                             group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=TRUE),max(plotdata[,yearvar],na.rm=TRUE),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            #geom_text_repel(aes(label = trend.label),nudge_x=0.01,nudge_y=0.01,show.legend = FALSE)
            geom_text_repel(aes(label = plotdata[, "trend.label"]),show.legend = FALSE)
        }else{
          plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
            geom_point(data=(subset(plotdata, plotdata[, "obsNA"]==0)), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obsintvar_1minus)), 
                                                             group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=TRUE),max(plotdata[,yearvar],na.rm=TRUE),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            #geom_text_repel(aes(label = trend.label),nudge_x=0.01,nudge_y=0.01,show.legend = FALSE)
            geom_text_repel(aes(label = plotdata[, "trend.label"]),show.legend = FALSE)
        }
      }
    }
    if(length(interval.values)>3 | nJP>3){
      if(annotation==1){
        stop("This annotation feature is not available when the number of joinpoints>3 or the number of multiple intervals selected>3.")
      }
    }
  }
  plot0<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
    geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
    geom_point(data=(subset(plotdata, plotdata[, "obsNA"]==0)), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obsintvar_1minus)), 
                                                     group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
    xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
    theme_bw()+
    ggtitle(title.rch) +
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
    scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=TRUE),max(plotdata[,yearvar],na.rm=TRUE),5))+
    scale_color_hue(labels = interval.labels)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
  
  if(!annotation %in% c(0,1)){
    stop("The annotation indicator should be either 0 or 1.")
  }
  if(!topanno %in% c(0,1)){
    stop("The annotation indicator for the top curve only should be either 0 or 1.")
  }
  
  if(trend==0){
    if(annotation==0 | length(interval.values)>3 | nJP>3){
      out<-list(plot0)
      names(out)<-"plot_no_anno"
    }else{
      out<-list(plot,plot0)
      names(out)<-c("plot_anno","plot_no_anno")
    }
  }else if(trend==1){
    if(annotation==0 | length(interval.values)>3 | nJP>3){
      out<-list(trends.out,plot0)
      names(out)<-c("trends","plot_no_anno")
    }else{
      out<-list(trends.out,plot,plot0)
      names(out)<-c("trends","plot_anno","plot_no_anno")
    }
  }else{
    stop("The trend indicator should be either 0 or 1.")
  }
  return(out)
}

#########################################################################################################
# plot.surv.year.annotate: a function that returns a plot for Average Change in Cumulative Survival by diagnosis year
# using ggplot
# The annotation feature is available for nJP<=3 and the number of multiple intervals selected <=3
# arguments:
# plotdata - the graph data returned by function download.data with downloadtype="graph".
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obscumvar - the variable name for observed relative cumulative survival. The default is "Relative_Survival_Cum"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predcumvar - the variable name for predicted cumulative survival. The default is "Predicted_Survival_Cum".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# annotation - the indicator for the annotation feature. The default is 0 (no annotation on the plot).Two plots
#              with and without annotation will be returned in a list when annotation=1.
# trend - the indicator for returning the trend measure tables. The default is 0 (no trend tables returned).
#########################################################################################################

### define a plot function for Average Change in Cumulative Survival by diagnosis year with annotations
Plot.surv.year.annotate<-function(plotdata,fit,nJP,yearvar,obscumvar="Relative_Survival_Cum",predcumvar="Predicted_Survival_Cum",interval="Interval",annotation=0,trend=0,title=NULL){
  interval.values<-as.numeric(unique(plotdata[,interval]))
  trends.out<-list()
  rows.obs<-which(is.na(plotdata[,obscumvar]))
  plotdata[,"obsNA"]<-0
  plotdata[rows.obs,"obsNA"]<-1
  if(is.null(title)){
    if(grepl("Relative", obscumvar)==TRUE){
      title.acs<-"Relative Survival by Diagnosis Year"
      ylabel<-"Relative Survival (%)"
      interval.labels<-paste(interval.values,"-year Relative Survival",sep="")
    }else if(grepl("Cause", obscumvar)==TRUE){
      title.acs<-"Cause-Specific Survival by Diagnosis Year"
      ylabel<-"Cause-Specific Survival (%)"
      interval.labels<-paste(interval.values,"-year Cause-Specific Survival",sep="")
    }else{
      title.acs<-"Survival by Diagnosis Year"
      ylabel<-"Survival (%)"
      interval.labels<-paste(interval.values,"-year Survival",sep="")
    }
  }else{
    title.acs<-title
    if(grepl("Relative", obscumvar)==TRUE){
      ylabel<-"Relative Survival (%)"
      interval.labels<-paste(interval.values,"-year Relative Survival",sep="")
    }else if(grepl("Cause", obscumvar)==TRUE){
      ylabel<-"Cause-Specific Survival (%)"
      interval.labels<-paste(interval.values,"-year Cause-Specific Survival",sep="")
    }else{
      ylabel<-"Survival (%)"
      interval.labels<-paste(interval.values,"-year Survival",sep="")
    }
  }

  
  if(trend==1 | annotation==1){
    trends.out <- aapc.multiints(fit$FitList[[nJP+1]], type="AbsChgSur", int.select=interval.values)
    if(length(interval.values)<=3 & nJP<=3){
      ### haz results are the same for all interval values
      jp.loc<-fit$FitList[[nJP+1]]$jp
      x.values<-list()
      annot.strs<-list()
      plotdata[,"trend.label"]<-""
      for(i in 1:length(interval.values)){
        int.i<-interval.values[i]
        cs.aaac<-trends.out[[i]]
        annot.strs[[i]]<-paste(sprintf("%.1f",100*cs.aaac$estimate),sep="")
        plotdata.i<-plotdata[which(plotdata[,interval]==int.i),]
        if(max(plotdata.i[,yearvar])==max(cs.aaac$end.year) | nJP==0){
          end.year<-cs.aaac$end.year
          end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
          x.values.i<-1/2*(cs.aaac$start.year+end.year)
        }
        if(max(plotdata.i[,yearvar])<max(cs.aaac$end.year)){
          if(nJP==1){
            if(max(plotdata.i[,yearvar])>jp.loc){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }else{
              cs.aaac<-cs.aaac[1:nJP,]
              x.values.i<-1/2*(cs.aaac$start.year+cs.aaac$end.year)
            }
          } ## nJP=1 end
          if(nJP==2){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>min(jp.loc)){
              cs.aaac<-cs.aaac[1:nJP,]
              x.values.i<-1/2*(cs.aaac$start.year+cs.aaac$end.year)
            } else{
              cs.aaac<-cs.aaac[1:(nJP-1),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }
          } ## nJP=2 end
          if(nJP==3){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              cs.aaac<-cs.aaac[1:nJP,]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              cs.aaac<-cs.aaac[1:(nJP-1),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else {
              cs.aaac<-cs.aaac[1:(nJP-2),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }
          } ## nJP=3 end
        } 
        x.values[[i]]<-round(x.values.i)
        plotdata$trend.label[which((plotdata[,yearvar] %in% x.values[[i]]) & plotdata[,interval]==interval.values[[i]])]<-annot.strs[[i]][1:length(x.values[[i]])]
      }## interval.values iteration end
      if(annotation==1){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
          geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
          geom_point(data=(subset(plotdata, plotdata[, "obsNA"]==0)), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obscumvar)), 
                                                           group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
          xlab('Year at diagnosis') + ylab(ylabel)+
          theme_bw()+
          ggtitle(title.acs) +
          coord_cartesian(ylim=c(0,1))+
          scale_y_continuous(breaks=seq(0,1,0.1), labels=scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=TRUE),max(plotdata[,yearvar],na.rm=TRUE),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          #geom_text_repel(aes(label = trend.label),nudge_x=0.01,nudge_y=0.01,show.legend = FALSE)
          geom_text_repel(aes(label = plotdata[, "trend.label"]),show.legend = FALSE)
      }
    }
    if(length(interval.values)>3 | nJP>3){
      if(annotation==1){
        stop("This annotation feature is not available when the number of joinpoints>3 or the number of multiple intervals selected>3.")
      }
    }
  }
  
  plot0 <- ggplot(data=plotdata, 
                  aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],
                    group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
           geom_line(data=plotdata, 
                  aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), 
                      colour=as.factor(plotdata[,interval])), linetype="solid") +
           geom_point(data=(subset(plotdata, plotdata[, "obsNA"]==0)), 
                  aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obscumvar)), 
                      group=as.factor(eval(parse(text=interval))), 
                      colour=as.factor(eval(parse(text=interval))))) +
           xlab('Year at diagnosis') + ylab(ylabel) +
           theme_bw() +
           ggtitle(title.acs) +
           coord_cartesian(ylim=c(0,1)) +
           scale_y_continuous(breaks=seq(0,1,0.1), labels=scales::percent) +
           scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=TRUE),max(plotdata[,yearvar],na.rm=TRUE),5)) +
           scale_color_hue(labels = interval.labels) +
           theme(legend.position="bottom") +
           theme(legend.title=element_blank()) +
           theme(plot.title = element_text(hjust = 0.5))

  if(!annotation %in% c(0,1)){
    stop("The annotation indicator should be either 0 or 1.")
  }
  if(trend==0){
    if(annotation==0 | length(interval.values)>3 | nJP>3){
      out<-list(plot0)
      names(out)<-"plot_no_anno"
    }else{
      out<-list(plot,plot0)
      names(out)<-c("plot_anno","plot_no_anno")
    }
  }else if(trend==1){
    if(annotation==0 | length(interval.values)>3 | nJP>3){
      out<-list(trends.out,plot0)
      names(out)<-c("trends","plot_no_anno")
    }else{
      out<-list(trends.out,plot,plot0)
      names(out)<-c("trends","plot_anno","plot_no_anno")
    }
  }else{
    stop("The trend indicator should be either 0 or 1.")
  }
  return(out)
}

#########################################################################################################
# plot.surv.int.multiyears: a function that returns a plot for Cumulative Survival by Interval supporting
# multiple selected years
# arguments:
# plotdata - 1. the graph data returned by function download.data with downloadtype="graph" and all intervals needed for int.select.
#            2. the full data returned by function download.data with downloadtype="full".
#            Either of the above would work.
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obscumvar - the variable name for observed relative cumulative survival. The default is "Relative_Survival_Cum"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predcumvar - the variable name for predicted cumulative survival. The default is "Predicted_Survival_Cum".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# year.select - the year values selected for the plot. The default is NULL.
#########################################################################################################

### define a plot function for Cumulative Survival by Interval
Plot.surv.int.multiyears<-function(plotdata,fit,nJP,yearvar,obscumvar="Relative_Survival_Cum",predcumvar="Predicted_Survival_Cum",interval="Interval",year.select=NULL){
  if(min(year.select)<min(plotdata[,yearvar],na.rm=TRUE)){
    stop("The min value of selected years should be greater than the min value of years in the plotdata.")      ### 11/25
  }else if(max(year.select)>max(plotdata[,yearvar],na.rm=TRUE)){
    stop("The max value of selected years should be less than the max value of years in the plotdata.") 
  }else{
    year.select.str<- paste(year.select,collapse=", ",sep="")
    if(length(year.select)<=3){
      if(grepl("Relative", obscumvar)==TRUE){
        title.survint<-paste("Relative Survival by Interval for ",year.select.str,sep="")
        ylabel<-"Relative Survival (%)"
      }else if(grepl("Cause", obscumvar)==TRUE){
        title.survint<-paste("Cause-Specific Survival by Interval for ",year.select.str,sep="")
        ylabel<-"Cause-Specific Survival (%)"
      }else{
        title.survint<-paste("Survival by Interval for ",year.select.str,sep="")
        ylabel<-"Survival (%)"
      }
    }else{
      if(grepl("Relative", obscumvar)==TRUE){
        title.survint<-"Relative Survival by Interval"
        ylabel<-"Relative Survival (%)"
      }else if(grepl("Cause", obscumvar)==TRUE){
        title.survint<-"Cause-Specific Survival by Interval"
        ylabel<-"Cause-Specific Survival (%)"
      }else{
        title.survint<-"Survival by Interval"
        ylabel<-"Survival (%)"
      }
    }
    x.combo<-function(x){
      paste(x,collapse="_",sep="")
    }
    cohort.combo <- rep(NA,dim(plotdata)[1])
    cohort.combo <- apply(data.frame(plotdata[,1:(which(colnames(plotdata)==interval)-1)]), 1, x.combo)
    int.2d<-sprintf("%02d",plotdata[,interval])
    cohort.combo.int<-paste(cohort.combo,"_",int.2d,sep="")
    plotdata.add<-plotdata[which(plotdata[,interval]==1),]
    plotdata.add[,interval]<-0
    plotdata.add[,c(obscumvar,predcumvar)]<-1
    cols.keep<-c(colnames(plotdata)[1:which(colnames(plotdata)==interval)],obscumvar,predcumvar)
    plotdata.add[,colnames(plotdata)[which(!colnames(plotdata) %in% cols.keep)]]<-NA
    plotdata[,"cohort.combo.int"]<-cohort.combo.int
    cohort.combo.add <- apply(data.frame(plotdata.add[,1:(which(colnames(plotdata.add)==interval)-1)]), 1, x.combo)
    cohort.combo.add.int0<-paste(cohort.combo.add,"_00",sep="")
    plotdata.add[,"cohort.combo.int"]<-cohort.combo.add.int0
    plotdata.merge<-rbind(plotdata,plotdata.add)
    plotdata.sort<-plotdata.merge[order(plotdata.merge$cohort.combo.int),]
    plotdata.sub<-plotdata.sort[which((plotdata.sort[,yearvar] %in% year.select) & !is.na(plotdata.sort[,obscumvar])),]
    int.max<-max(plotdata.sub[,interval],na.rm=TRUE)
    if(int.max<=12){
      xbreaks.by<-1
    }else{
      xbreaks.by<-2
    }
    plot<-ggplot(data=plotdata.sub, aes(x=plotdata.sub[,interval], y=plotdata.sub[,predcumvar],group=as.factor(plotdata.sub[,yearvar]), colour=as.factor(plotdata.sub[,yearvar]))) + 
      geom_line(data=plotdata.sub, aes(x=plotdata.sub[,interval], y=plotdata.sub[,predcumvar],group=as.factor(plotdata.sub[,yearvar]), colour=as.factor(plotdata.sub[,yearvar])),linetype="solid")+
      geom_point(data=plotdata.sub, aes(x=plotdata.sub[,interval], y=plotdata.sub[,obscumvar], group=as.factor(plotdata.sub[,yearvar]), colour=as.factor(plotdata.sub[,yearvar])))+
      xlab('Interval') + ylab(ylabel)+
      theme_bw()+
      ggtitle(title.survint) +
      coord_cartesian(ylim=c(0,1))+
      scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
      scale_x_continuous(breaks=seq(min(plotdata.sub[,interval],na.rm=TRUE),max(plotdata.sub[,interval],na.rm=TRUE),xbreaks.by))+
      scale_color_hue(labels = year.select)+
      theme(legend.position="bottom")+
      theme(legend.title=element_blank())+
      theme(plot.title = element_text(hjust = 0.5))
  }
}

#########################################################################################################
#### Functions Edited by Daniel Miller
# Combination of data input and model execution.  
joinpoint.seerdata = function(seerfilename, newvarnames,UseVarLabelsInData=FALSE) {
  
  #read data from seer files
  seerdata = .read.SeerStat(seerfilename,UseVarLabelsInData=UseVarLabelsInData);
  
  SessionAttr = attr(seerdata,"DICInfo");
  SessionType = SessionAttr$SessionOptionInfo[3,2]; #Either "Relative Survival" or "Cause-Specific Survival"
  if(SessionType == "Cause-Specific Survival"){
    Expected_Survival_Interval = rep(1,length(seerdata$CauseSpecific_Survival_Cum)); 
    seerdata = cbind(seerdata,Expected_Survival_Interval);
  }
  
  varnames=c(newvarnames);
  varnames=factor(varnames);
  varnames=levels(varnames);
  #tempdata = attr(seerdata,"assignColNames")(seerdata,varnames);
  
  #colnames(seerdata) = colnames(tempdata);
  
  return(seerdata);
  
}

# Function to return the information from the .dic file for parsing purposes.
# Added 12/12/2014 - DM
dictionary.overview = function(DICfilename) {
  
  DICdata = .read.SeerStatDIC(DICfilename);
  
  DICoverview = list(
    SystemInfo = DICdata$SystemInfo,
    SessionOptionInfo = DICdata$SessionOptionInfo,
    ExportOptionInfo = DICdata$ExportOptionInfo,
    VarAllInfo = DICdata$VarAllInfo,
    VarFormatSecList = DICdata$VarFormatSecList,
    VarLabelInfo=DICdata$VarLabelInfo,
    VarWithoutFormatItem = DICdata$VarWithoutFormatItem
  );
  
  return(DICoverview);
  
}

# Function to return the value of the Life Page variable in the loaded dictionary file.
# Added 1/23/2015 - DM
find.life.page.val = function(pagetypeDF) {
  
  row_num <- nrow(pagetypeDF);
  for (i in 1:row_num ) {
    if(pagetypeDF[i,2] == "Life Page"){
      lpvalueSTR <- pagetypeDF[i,1];
    };
  };
  
  lpvalue <- as.numeric(lpvalueSTR);
  
  return(lpvalue);
  
}
