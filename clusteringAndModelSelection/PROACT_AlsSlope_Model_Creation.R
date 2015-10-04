source(file=".../AlsSlopeCrossValidationUtilities.R")
load(".../slope_features_default_V05.rdata") 

### hardcoded classification
classifyPatient = function(onset_delta, Q1_Speech, trunk, hands){
	label = -1;
	if(onset_delta >= 0.04964584){
		if(Q1_Speech < 0.66757882){
			if(onset_delta>=0.053412932){
				if(hands < 2.01044607){
					label = 111
				} else{
					label = 112
				}
			} else{
				label = 12
			}
		} else{
			if(trunk < 0.58163565){
				label = 21
			} else{
				label = 22
			}
		}
	} else{
		if(onset_delta >= -0.03113981){
			label = 3
		} else{
			label = 4
		}
	}
}

getGroupLabels = function(table){
	labels = c()
	for(i in 1:length(table[,1])){
		row = table[i,]
		od = row$onset_delta[1]
		Q1 = row$Q1_Speech[1]
		tr = row$trunk[1]
		ha = row$hands[1]
		label = classifyPatient(onset_delta=od, Q1_Speech=Q1, trunk=tr, hands=ha)
		labels[i] = label
	}
	return(labels)
}

doesPatientInFeatureTabHaveSlope = !is.na(feature_tab$ALSFRS_slope)
feature_tab$CART_label = getGroupLabels(feature_tab)
train_tab = feature_tab[doesPatientInFeatureTabHaveSlope,]



### perform stepwise forward / backward selection
currentOut <- list()
for(grp in c(111,112,12,21,22,3,4)){
	cat(grp, " ###########################################################\n")
	grp_tab = train_tab[train_tab$CART_label==grp,]
	grp_tab = grp_tab[!is.na(grp_tab$ALSFRS_slope),]
	init = NULL
	
	pool=NULL
	p_rg=0
	
	if(grp == 111){
		init = c("Q9_Climbing_Stairs","treatment_group","Phosphorus")
	}
	if(grp == 112){
		init = c("Gender","Creatinine","Alkaline.Phosphatase","White.Blood.Cell..WBC.","Creatinine:White.Blood.Cell..WBC.")
	}
	if(grp == 12){
		init = c("ALSFRS_Total","CK","ALSFRS_Total*CK", "fvc_normal")
	}
	if(grp == 21){
		init = c("onset_delta","bp_systolic","Q10_Respiratory","Alkaline.Phosphatase")
	}
	if(grp == 22){
		init = c("onset_delta","Q5_Cutting","fvc_percent")
	}
	if(grp == 3){
		init = c("mouth","trunk","fvc","mouth:trunk")
	}
	if(grp == 4){
		init = c("trunk","Race","bp_diastolic","if_use_Riluzole")
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
for(grp in c(111,112,12,21,22,3,4)){
	cat(grp, " ###########################################################\n")
	test=data.frame(currentOut[[grp]],stringsAsFactors=FALSE)
	test$model_RMSE_mean=as.numeric(test$model_RMSE_mean)
	test=test[with(test,order(model_RMSE_mean, decreasing = FALSE)),]
	print(head(test,n=20))
}
