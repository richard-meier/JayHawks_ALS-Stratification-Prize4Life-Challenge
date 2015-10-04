source(file=".../AlsSlopeCrossValidationUtilities.R")
load(".../slope_features_default_reg.rdata")

### hardcoded classification
classifyPatient = function(onset_delta, priorSlope){
	label = -1;
	if(onset_delta >= 0.1464){
		if(priorSlope < 0.4279){
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
		od = row$onset_delta[1]
		ps = row$priorSlope[1]
		label = classifyPatient(onset_delta=od, priorSlope=ps)
		labels[i] = label
	}
	return(labels)
}

doesPatientInFeatureTabHaveSlope = !is.na(feature_tab$ALSFRS_slope)
feature_tab$CART_label = getGroupLabels(feature_tab)
train_tab = feature_tab[doesPatientInFeatureTabHaveSlope,]



### perform stepwise forward / backward selection
currentOut <- list()
for(grp in c(1,2,3)){
	cat(grp, " ###########################################################\n")
	grp_tab = train_tab[train_tab$CART_label==grp,]
	grp_tab = grp_tab[!is.na(grp_tab$ALSFRS_slope),]
	init = NULL
	
	pool=NULL
	p_rg=0
	
	if(grp == 1){
		
	  init = c("ALS_Staging_Com", "ALSFRS_Total", "Q1_Speech", "Q9_Climbing_Stairs", "hands")
	  }
	if(grp == 2){
		init = c("Q1_Speech", "Gender", "Q4_Handwriting", "MEDHx_Diabetes", "ALS_Staging_Total", "ALS_Staging_Move") 
	}
	if(grp == 3){
		init = c("onset_site", "ALS_Staging_Eat", "priorSlope", "Q6_Dressing_and_Hygiene", "trunk", "Gender")
	}
	
	N0 = length(grp_tab[,1])
	cat("current group size = ",N0,"\n")
	crg = round(0.333334 * N0)
	cat("test set size = ",crg,"\n")
	
	## EITHER INCLUSION
	exclude = c("SubjectID","ALSFRS_slope","time_event","status","CART_label")
	currentAdd = determineFeatureToAdd(train_data=grp_tab, initial_model=init, initial_exclusion=exclude, vFold=5000, c_range=crg, addPool=pool, p_range=p_rg)
	currentOut[[grp]] = currentAdd
	
	## OR EXCLUSION
	#determineFeatureToRemove(train_data=grp_tab, initial_model=init, vFold=5000, c_range=crg)
}

par(mfrow=c(2,2))
for(grp in c(1,2,3)){
	print(grp); cat(" ##########################################\n")
	test=data.frame(currentOut[[grp]],stringsAsFactors=FALSE)
	test$model_RMSE_mean=as.numeric(test$model_RMSE_mean)
	test=test[with(test,order(model_RMSE_mean, decreasing = FALSE)),]
	print(head(test,n=20))
}
