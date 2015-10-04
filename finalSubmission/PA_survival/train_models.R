library(survival)
library(pec)

load(".../slope_features_default_V05.rdata") 
outFile = ".../clusterModels_SPA.RData"

# hardcoded classification

classifyPatient = function(fvc_percent, onsetAge, ALSFRS_Total){
	label = -1;
	if(fvc_percent >= -0.464656830){
		if(onsetAge < -0.007705498){
			label=1
		} else{
			if(ALSFRS_Total >= 0.077406712){
				label=2
			} else{
				label=3
			}
		}
	} else{
		label=4
	}
}

getGroupLabels = function(table){
	labels = c()
	for(i in 1:length(table[,1])){
		row = table[i,]
		fp = row$fvc_percent[1]
		oa = row$onsetAge[1]
		at = row$ALSFRS_Total[1]
		label = classifyPatient(fvc_percent=fp, onsetAge=oa, ALSFRS_Total=at)
		labels[i] = label
	}
	return(labels)
}

doesPatientInFeatureTabHaveSurvival = !is.na(feature_tab$time_event)
feature_tab$CART_label = getGroupLabels(feature_tab)
train_tab = feature_tab[doesPatientInFeatureTabHaveSurvival,]

timePoints = c(12, 18, 24) * 30.5

grp_tab1 = train_tab[train_tab$CART_label==1,]
modelForCluster_1 = coxph(Surv(time = time_event, event = status) ~ fvc + priorSlope + onsetAge + fvc_normal + Chloride, data=grp_tab1)
probs = predictSurvProb(modelForCluster_1, newdata=grp_tab1, times=timePoints)
avgProbs1 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

grp_tab2 = train_tab[train_tab$CART_label==2,]
modelForCluster_2 = coxph(Surv(time = time_event, event = status) ~ trunk + priorSlope + onsetAge + fvc_percent, data=grp_tab2)
probs = predictSurvProb(modelForCluster_2, newdata=grp_tab2, times=timePoints)
avgProbs2 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

grp_tab3 = train_tab[train_tab$CART_label==3,]
modelForCluster_3 = coxph(Surv(time = time_event, event = status) ~ fvc + priorSlope + onsetAge + fvc_percent + Gender, data=grp_tab3)
probs = predictSurvProb(modelForCluster_3, newdata=grp_tab3, times=timePoints)
avgProbs3 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

grp_tab4 = train_tab[train_tab$CART_label==4,]
modelForCluster_4 = coxph(Surv(time = time_event, event = status) ~ fvc + priorSlope + onsetAge + Hemoglobin + Chloride, data=grp_tab4)
probs = predictSurvProb(modelForCluster_4, newdata=grp_tab4, times=timePoints)
avgProbs4 = c(mean(probs[,1]),mean(probs[,2]),mean(probs[,3]))

save(
	modelForCluster_1, modelForCluster_2, modelForCluster_3, modelForCluster_4, 
	grp_tab1, grp_tab2, grp_tab3, grp_tab4,
	avgProbs1, avgProbs2, avgProbs3, avgProbs4,
	file=outFile
)
