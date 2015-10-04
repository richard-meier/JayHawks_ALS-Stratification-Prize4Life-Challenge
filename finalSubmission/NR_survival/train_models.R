library(survival)
library(pec)

load(".../slope_features_default_reg.rdata") 
outFile = ".../clusterModels_SNR.RData"

# hardcoded classification

classifyPatient = function(priorSlope, onset_delta){
	label = -1;
	if(priorSlope >= -0.001987586){
		if(onset_delta < 0.132676146){
			label=1
		} else{
			label=2
		}
	} else{
		label=3
	}
}

getGroupLabels = function(table){
	labels = c()
	for(i in 1:length(table[,1])){
		row = table[i,]
		ps = row$priorSlope[1]
		od = row$onset_delta[1]
		label = classifyPatient(priorSlope=ps, onset_delta=od)
		labels[i] = label
	}
	return(labels)
}

doesPatientInFeatureTabHaveSurvival = !is.na(feature_tab$time_event)
feature_tab$CART_label = getGroupLabels(feature_tab)
train_tab = feature_tab[doesPatientInFeatureTabHaveSurvival,]

timePoints = c(12, 18, 24) * 30.5

grp_tab1 = train_tab[train_tab$CART_label==1,]
modelForCluster_1 = coxph(Surv(time = time_event, event = status) ~ priorSlope + onsetAge + ALS_Staging_Total + priorSlope*ALS_Staging_Total + Age + respiratory, grp_tab1)
probs = predictSurvProb(modelForCluster_1, newdata=grp_tab1, times=timePoints)
avgProbs1 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

grp_tab2 = train_tab[train_tab$CART_label==2,]
modelForCluster_2 = coxph(Surv(time = time_event, event = status) ~ priorSlope + onsetAge + Q1_Speech + Gender, grp_tab2)
probs = predictSurvProb(modelForCluster_2, newdata=grp_tab2, times=timePoints)
avgProbs2 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

grp_tab3 = train_tab[train_tab$CART_label==3,]
modelForCluster_3 = coxph(Surv(time = time_event, event = status) ~ priorSlope + onsetAge + hands + MEDHx_Thyroid, grp_tab3)
probs = predictSurvProb(modelForCluster_3, newdata=grp_tab3, times=timePoints)
avgProbs3 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

save(
	modelForCluster_1, modelForCluster_2, modelForCluster_3, 
	grp_tab1, grp_tab2, grp_tab3,
	avgProbs1, avgProbs2, avgProbs3,
	file=outFile
)
