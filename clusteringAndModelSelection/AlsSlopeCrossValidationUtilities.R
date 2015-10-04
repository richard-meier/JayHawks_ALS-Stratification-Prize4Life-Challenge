
rmse <- function(error) {
    sqrt(mean(error^2))
}

determineFeatureToRemove <- function(train_data, initial_model, vFold=10000, c_range=534, addPool=NULL, p_range=0){
	n = length(train_data[,1])
	exclusion=c()
	models=c()
	RMSE_means=c()
	for(removedFeature in initial_model){
		modelString=""
		for(i in 1:length(initial_model)){
			if(initial_model[i] != removedFeature){
				modelString=paste(modelString,initial_model[i]," + ",sep="")
			}
		}
		modelString = substr(modelString, start=1, stop=nchar(modelString)-3)
		modelString=paste("ALSFRS_slope ~ ", modelString, sep="")
		aucs=c()
		for( ii in 1:vFold ){
			index <- sample(1:n)[1:c_range]
			sub_train <- train_data[-index,]
			categoricalVoid_train = drawVoidSet(sub_train, c("Race", "onset_site", "Gender", "treatment_group", "if_use_Riluzole"))
			sub_train <- rbind(sub_train, categoricalVoid_train)
			sub_test <- train_data[index,]
			if(!is.null(addPool)){
				np = length(addPool[,1])
				sub_p = addPool[ sample(1:np)[1:p_range], ]
				sub_test = rbind(sub_test, sub_p)
			}
			fit = lm(
				eval(parse(text=modelString)),
				data = sub_train, na.action = na.omit
			)
			pred <- predict( fit, newdata = sub_test, interval="predict")
			obs <- sub_test$ALSFRS_slope
			error = obs - pred[,1]
			aucs=c(aucs, rmse(error))
		}
		RMSE = mean(aucs)
		cat(modelString,"\tRMSE_mean=",RMSE,"\n",sep="")
		RMSE_means=c(RMSE_means,RMSE)
		models=c(models,modelString)
		exclusion=c(exclusion,removedFeature)
	}
	return(cbind(exclusion,models,RMSE_means))
}

drawVoidSet <- function(original_data, features){
	out = original_data[original_data$Gender=="NonExistentDummy",] # draw empty table
	n = length(original_data[,1])
	uID = -100
	for(feat in features){
		featVect = unlist(original_data[feat])
		present = names(table(featVect)[table(featVect)>0])
		full = levels(featVect)
		missing = setdiff(full, present)
		for( lvl in missing ){
			sub_rows = original_data[sample(1:n)[1],]
			uID = uID -1
			sub_rows$SubjectID[1]=uID
			sub_rows[1,feat] = lvl
			sub_rows$ALSFRS_slope[1] = sub_rows$ALSFRS_slope[1]
			out=rbind(out,sub_rows)
		}
	}
	return(out)
}

determineFeatureToAdd <- function(train_data, initial_model, initial_exclusion, vFold=500, c_range=534, addPool=NULL, p_range=0){
	exclude=c(initial_exclusion)
	feature_set=setdiff(names(train_data),exclude)
	
	maxModel=""
	minRMSE=9999999999
	
	model_str=c()
	model_RMSE_mean=c()
	
	n = length(train_data[,1])
	
	for(addedFeature in feature_set){
		cat("PROCESSING FEATURE:",addedFeature,"... ")
		modelString=paste(initial_model[1],sep="")
		for(i in 2:length(initial_model)){
			modelString=paste(modelString," + ",initial_model[i],sep="")
		}
		modelString=paste(modelString," + ",addedFeature,sep="")
		
		subM0=paste("ALSFRS_slope ~ ",modelString,sep="")
		Submodels=c(subM0)
		frequencyTable = sort(table(train_data[addedFeature]))
		abundanceRate = frequencyTable[length(frequencyTable)] / sum(frequencyTable)
		cat("abundance rate =",abundanceRate,"\n")
		if(abundanceRate<0.77){
			for(i in 1:length(initial_model)){
				submodel=paste(modelString," + ",initial_model[i],"*",addedFeature,sep="")
				submodel=paste("ALSFRS_slope ~ ",submodel,sep="")
				Submodels=c(Submodels, submodel)
			}
		}
		else{
			cat("--> HIGH ABUNDANCE: SKIPPING INTERACTIONS!\n",sep="")
		}
		for(submodel in Submodels){
			aucs=c()
			for( ii in 1:vFold){
				index <- sample(1:n)[1:c_range]
				sub_train <- train_data[-index,]
				categoricalVoid_train = drawVoidSet(sub_train, c("Race", "onset_site", "Gender", "treatment_group", "if_use_Riluzole"))
				sub_train <- rbind(sub_train, categoricalVoid_train)
				sub_test <- train_data[index,]
				if(!is.null(addPool)){
					np = length(addPool[,1])
					sub_p = addPool[ sample(1:np)[1:p_range], ]
					sub_test = rbind(sub_test, sub_p)
				}
				fit = lm(
					eval(parse(text=submodel)),
					data = sub_train, na.action = na.omit
				)
				pred <- predict( fit, newdata = sub_test, interval="predict" )
				obs <- sub_test$ALSFRS_slope
				error = obs - pred[,1]
				aucs=c(aucs, rmse(error))
			}
			RMSE = mean(aucs)
			model_str=c(model_str, submodel)
			model_RMSE_mean=c(model_RMSE_mean, RMSE)
			if(minRMSE>RMSE){
				minRMSE=RMSE
				maxModel=submodel
				cat("New best model!!!\n\t",maxModel,"\n\tmeanRMSE=",minRMSE,"\n",sep="")
			}
		}
	}
	return( cbind(model_str, model_RMSE_mean) )
}
