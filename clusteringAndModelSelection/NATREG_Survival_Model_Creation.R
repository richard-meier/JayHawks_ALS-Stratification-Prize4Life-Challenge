source(file=".../SurvivalCrossValidationUtilities.R")
load(".../slope_features_default_reg.rdata") 

### hardcoded classification
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

doesPatientInFeatureTabHaveSurvivalData = !is.na(feature_tab$time_event)
feature_tab$CART_label = getGroupLabels(feature_tab)
train_tab = feature_tab[doesPatientInFeatureTabHaveSurvivalData,]



### perform stepwise forward / backward selection
currentOut <- list()
for(grp in c(1,2,3)){
	print(grp); cat(" ##########################################\n")
	grp_tab = train_tab[train_tab$CART_label==grp,]
	init = NULL
	
	pool=NULL
	p_rg=0
	
	if(grp == 1){
		init = c("priorSlope","onsetAge", "ALS_Staging_Total", "priorSlope*ALS_Staging_Total")
	}
	if(grp == 2){
		init = c("priorSlope","onsetAge", "Q1_Speech")
	}
	if(grp == 3){
		init = c("priorSlope","onsetAge", "hands")
	}
	
	N0 = length(grp_tab[,1])
	cat("current group size = ",N0,"\n")
	crg = round(0.333334 * N0)
	cat("test set size = ",crg,"\n")
	
	## EITHER INCLUSION
	exclude = c("SubjectID","ALSFRS_slope","time_event","status","CART_label")
	currentAdd = determineFeatureToAdd(train_data=grp_tab, initial_model=init, initial_exclusion=exclude, vFold=5000, c_range=crg)
	currentOut[[grp]] = currentAdd
	
	## OR EXCLUSION
	#determineFeatureToRemove(train_data=grp_tab, initial_model=init, vFold=5000, c_range=crg)
}

par(mfrow=c(2,2))
for(grp in c(1,2,3)){
	print(grp); cat(" ##########################################\n")
	test=data.frame(currentOut[[grp]],stringsAsFactors=FALSE)
	test$mean_cIdxs=as.numeric(test$mean_cIdxs)
	test=test[with(test,order(mean_cIdxs, decreasing = TRUE)),]
	print(head(test,n=20))
}
