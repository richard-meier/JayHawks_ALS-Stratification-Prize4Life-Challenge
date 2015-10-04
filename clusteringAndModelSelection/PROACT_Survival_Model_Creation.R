source(file=".../SurvivalCrossValidationUtilities.R")
load(".../slope_features_default_V05.rdata") 

### hardcoded classification
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

doesPatientInFeatureTabHaveSurvivalData = !is.na(feature_tab$time_event)
feature_tab$CART_label = getGroupLabels(feature_tab)
train_tab = feature_tab[doesPatientInFeatureTabHaveSurvivalData,]



### perform stepwise forward / backward selection
currentOut <- list()
for(grp in c(1,2,3,4)){
	print(grp); cat(" ##########################################\n")
	grp_tab = train_tab[train_tab$CART_label==grp,]
	init = NULL
	
	pool=NULL
	p_rg=0
	
	if(grp == 1){
		init = c("fvc","priorSlope","onsetAge")
	}
	if(grp == 2){
		init = c("trunk","priorSlope","onsetAge")
	}
	if(grp == 3){
		init = c("fvc","priorSlope","onsetAge")
	}
	if(grp == 4){
		init = c("fvc","priorSlope","onsetAge")
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

for(grp in c(1,2,3,4)){
	print(grp); cat(" ##########################################\n")
	test=data.frame(currentOut[[grp]],stringsAsFactors=FALSE)
	test$mean_cIdxs=as.numeric(test$mean_cIdxs)
	test=test[with(test,order(mean_cIdxs, decreasing = TRUE)),]
	print(head(test,n=20))
}
