
library(Bolstad2); library(ROCR); library(survival); library(timeROC); library(devtools); library(dismo);

### Scoring function as provided by the Prostate Cancer DREAM Challenge Administrators (March 16, 2015, https://www.synapse.org/#!Synapse:syn2813558/wiki/70844)
# question 1A
# riskScoreGlobal is the global risk score
# riskScore12, riskScore18, riskScore24 are the risk scores at 12, 18 and 24 months
# time is called time_event in the CoreTable meaning the last known follow up time in days
# death is last known follow up status (F=survival, T=death)
# all input parameters are vectors
# returned value is a vector containing:
#  * concordance index
#  * AUC of ROC at 12, 18, and 24 months
#  * integrated AUC (integrated over all time points)
## from Justin: https://mail.google.com/mail/u/1/#inbox/14c8f82a26d23b80
score_q1a<-function(time, death, riskScoreGlobal, riskScore12, riskScore18, riskScore24)
{

  if (missing(riskScore12)) {
    riskScore12 <- riskScoreGlobal
  }
  if (missing(riskScore18)) {
    riskScore18 <- riskScoreGlobal
  }
  if (missing(riskScore24)) {
    riskScore24 <- riskScoreGlobal
  }

  auc12 <- timeROC(T=time,
                  delta=death,
                  marker=riskScore12,
                  cause=1,
                  weighting="marginal",
                  times=12 * 30.5,
                  iid=FALSE)$AUC[2]

  auc18 <- timeROC(T=time,
                   delta=death,
                   marker=riskScore18,
                   cause=1,
                   weighting="marginal",
                   times=18 * 30.5,
                   iid=FALSE)$AUC[2]

  auc24 <- timeROC(T=time,
                   delta=death,
                   marker=riskScore24,
                   cause=1,
                   weighting="marginal",
                   times=24 * 30.5,
                   iid=FALSE)$AUC[2]

  # compute global concordance index
  surv <- Surv(time,death)
  cIndex <- survConcordance(surv ~ riskScoreGlobal)$concordance

  # compute iAUC from 0 to 30 months
  times <- seq(0,30,by=1) * 30.5
  aucs <- timeROC(T=time,
                  delta=death,
                  marker=riskScoreGlobal,
                  cause=1,
                  weighting="marginal",
                  times=times,
                  iid=FALSE)$AUC

  # Simpsons rules for integrating under curve
  iAUC <- sintegral(times, aucs)$int / (max(times) - min(times))

  return (list(cIndex=cIndex, auc12=auc12, auc18=auc18, auc24=auc24, iAUC=iAUC))
}

crossValidate = function(model, train_data, vFold=1000, c_range=534){
	aucs1=c()
	aucs12=c()
	aucs18=c()
	aucs24=c()
	cIdxs=c()
	n=length(train_data[,1])
	for( ii in 1:vFold){
		index <- sample(1:n)[1:c_range]
		sub_train <- train_data[-index,]
		sub_test <- train_data[index,]
		cox_sub = coxph(
			eval(parse(text=model)),
			data = sub_train, na.action = na.omit
		)
		risk_scores <- predict( cox_sub, type = "risk", newdata = sub_test )
		
		out_s = score_q1a( time = sub_test$time_event, death = sub_test$status, riskScoreGlobal = risk_scores)
		aucs1=c(aucs1, out_s$iAUC)
		aucs12=c(aucs12, out_s$auc12)
		aucs18=c(aucs18, out_s$auc18)
		aucs24=c(aucs24, out_s$auc24)
		cIdxs=c(cIdxs, out_s$cIndex)
	}
	cat(
		model, 
		"\tAUC12=", mean(na.omit(aucs12)), 
		"  AUC18=", mean(na.omit(aucs18)), 
		"  AUC24=", mean(na.omit(aucs24)), 
		"  cIndex=", mean(na.omit(cIdxs)), 
		"  iAUC=", mean(na.omit(aucs1)), 
		"\n",sep=""
	)
	return(data.frame(aucs1,aucs12,aucs18,aucs24,cIdxs))
}

determineFeatureToRemove <- function(train_data, initial_model, vFold=10000,c_range=534){
	exclusion=c()
	models=c()
	aucs1=c()
	aucs12=c()
	aucs18=c()
	aucs24=c()
	cIdxs=c()
	n=length(train_data[,1])
	for(removedFeature in initial_model){
		modelString="0"
		for(i in 1:length(initial_model)){
			if(initial_model[i] != removedFeature){
				modelString=paste(modelString," + ",initial_model[i],sep="")
			}
		}
		modelString=paste("Surv(time = time_event, event = status) ~ ", modelString, sep="")
		aucs1=c()
		aucs2=c()
		for( ii in 1:vFold ){
			index <- sample(1:n)[1:c_range]
			sub_train <- train_data[-index,]
			sub_test <- train_data[index,]
			cox_sub = coxph(
				eval(parse(text=modelString)),
				data = sub_train, na.action = na.omit
			)
			risk_scores <- predict( cox_sub, type = "risk", newdata = sub_test )
			
			out_s = score_q1a( time = sub_test$time_event, death = sub_test$status, riskScoreGlobal = risk_scores)
			aucs1=c(aucs1, out_s$iAUC)
			aucs12=c(aucs12, out_s$auc12)
			aucs18=c(aucs18, out_s$auc18)
			aucs24=c(aucs24, out_s$auc24)
			cIdxs=c(cIdxs, out_s$cIndex)
		}
		cat(
			modelString, 
			"\tAUC12=", mean(na.omit(aucs12)), 
			"  AUC18=", mean(na.omit(aucs18)), 
			"  AUC24=", mean(na.omit(aucs24)), 
			"  cIndex=", mean(na.omit(cIdxs)), 
			"  iAUC=", mean(na.omit(aucs1)), 
			"\n",sep=""
		)
	}
}

determineFeatureToAdd <- function(train_data, initial_model, initial_exclusion, vFold=500, c_range=534){
	exclude=c(initial_exclusion)
	feature_set=setdiff(names(train_data),exclude)
	
	maxModel=""
	maxCI=0
	
	model_str=c()
	mean_iaucs=c( )
	mean_aucs12=c( )
	mean_aucs18=c( )
	mean_aucs24=c( )
	mean_cIdxs=c( )
	
	n=length(train_data[,1])
	
	for(addedFeature in feature_set){
		cat("PROCESSING FEATURE:",addedFeature,"... ")
		modelString=paste(initial_model[1],sep="")
		for(i in 2:length(initial_model)){
			modelString=paste(modelString," + ",initial_model[i],sep="")
		}
		modelString=paste(modelString," + ",addedFeature,sep="")
		
		subM0=paste("Surv(time = time_event, event = status) ~ ",modelString,sep="")
		Submodels=c(subM0)
		frequencyTable = sort(table(train_data[addedFeature]))
		abundanceRate = frequencyTable[length(frequencyTable)] / sum(frequencyTable)
		cat("abundance rate =",abundanceRate,"\n")
		if(abundanceRate<0.77){
			for(i in 1:length(initial_model)){
				submodel=paste(modelString," + ",initial_model[i],"*",addedFeature,sep="")
				submodel=paste("Surv(time = time_event, event = status) ~ ",submodel,sep="")
				Submodels=c(Submodels, submodel)
			}
		}
		else{
			cat("--> HIGH ABUNDANCE: SKIPPING INTERACTIONS!\n",sep="")
		}
		for(submodel in Submodels){
			aucs1=c()
			aucs12=c()
			aucs18=c()
			aucs24=c()
			cIdxs=c()
			for( ii in 1:vFold){
				index <- sample(1:n)[1:c_range]
				sub_train <- train_data[-index,]
				sub_test <- train_data[index,]
				cox_sub = coxph(
					eval(parse(text=submodel)),
					data = sub_train, na.action = na.omit
				)
				risk_scores <- predict( cox_sub, type = "risk", newdata = sub_test )
				out_s = score_q1a( time = sub_test$time_event, death = sub_test$status, riskScoreGlobal = risk_scores)
				aucs1=c(aucs1, out_s$iAUC)
				aucs12=c(aucs12, out_s$auc12)
				aucs18=c(aucs18, out_s$auc18)
				aucs24=c(aucs24, out_s$auc24)
				cIdxs=c(cIdxs, out_s$cIndex)
			}
			
			model_str=c(model_str, submodel)
			
			miauc = mean(na.omit(aucs1))
			m12auc = mean(na.omit(aucs12))
			m18auc = mean(na.omit(aucs18))
			m24auc = mean(na.omit(aucs24))
			mcidxs = mean(na.omit(cIdxs))
			
			mean_iaucs=c(mean_iaucs, miauc)
			mean_aucs12=c(mean_aucs12, m12auc)
			mean_aucs18=c(mean_aucs18, m18auc)
			mean_aucs24=c(mean_aucs24, m24auc)
			mean_cIdxs=c(mean_cIdxs, mcidxs)
			
			if(maxCI<mcidxs){
				maxCI=mcidxs
				maxModel=submodel
				cat("New max model!!!\n\t",maxModel,"\n\tcIdx=",maxCI,"  iAUC=",miauc,"\n",sep="")
			}
		}
	}
	return( cbind(model_str, mean_iaucs, mean_aucs12, mean_aucs18, mean_aucs24, mean_cIdxs) )
}
