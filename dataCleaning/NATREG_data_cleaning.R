
if(!require(plyr)){
install.packages("plyr")
}
library(plyr)
#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute)

### define output paths
pathToSingularDataTable = ".../singularDataTable_reg.csv"
pathToSurvivalDataTable = ".../surv_response_registry_training.txt"
pathToDummyCodedFeatureTable = ".../slope_features_dummy_coded_reg.rdata"
pathToNormallyCodedFeatureTable = ".../slope_features_default_reg.rdata"

### load data
data.allforms<-read.delim(".../train/all_forms_registry_training.txt",sep="|", header=F)
data.ALSslope<-read.delim(".../train/ALSFRS_slope_registry_train.txt",sep="|", header=T)
colnames(data.allforms) = c("SubjectID","form_name","feature_name","feature_value","feature_unit","feature_delta")
colnames(data.ALSslope) = c("patients", "ALSFRS_slope")
slope_patients = as.character(unlist(data.ALSslope["patients"]))

data = data.allforms[data.allforms$form_name != "Adverse Event" & data.allforms$form_name != "Concomitant Medication",]
intervals = data.allforms[data.allforms$form_name == "Adverse Event" | data.allforms$form_name == "Concomitant Medication",]

feature_name = levels(factor(data$feature_name))
dataALS = vector("list")
for(i in feature_name){
  print(i)
  temp = data[data$feature_name == i,]
  temp = temp[ , -which(names(temp) %in% c("form_name","feature_name"))]
  for(j in 2:length(temp)){
    temp[,j] = factor(temp[,j])
  }
  dataALS[[i]] = temp
}

features_with_less_than_25perc_missing = c(
	"Age","ALS_Staging_Breath","ALS_Staging_Com","ALS_Staging_Eat","ALS_Staging_Move","ALS_Staging_Total","ALSFRS_R_Total",
	"ALSFRS_Total","family_ALS_hist","Gender","hands","leg","MEDHx_Diabetes","MEDHx_Hypertens","MEDHx_Thyroid","mouth",
	"Mutation","onset_delta","onset_site","Q1_Speech","Q2_Salivation","Q3_Swallowing","Q4_Handwriting","Q5_Cutting",
	"Q6_Dressing_and_Hygiene","Q7_Turning_in_Bed","Q8_Walking","Q9_Climbing_Stairs","R1_Dyspnea","R2_Orthopnea",
	"R3_Respiratory_Insufficiency","respiratory","respiratory_R","trunk"
)
dataALS_full=dataALS



### obtain all relevant patients and store them as a set
patientSet = c()
for(i in features_with_less_than_25perc_missing){
	featTable = dataALS_full[[i]]
	patient_sub = names(table((featTable$SubjectID)))
	patientSet = union(patientSet, patient_sub)
}



### delete entries consisting of a single dash character '-'
for(i in features_with_less_than_25perc_missing){
	featTable = dataALS_full[[i]]
	selection = featTable$feature_value=="-"
	selection[is.na(selection)]=FALSE
	selection = !selection
	featTable = featTable[selection,]
	#cat(i, length(featTable[,1]),"\n")
	dataALS_full[[i]] = featTable
}



### create early ALSFRS_slope

aTotal = dataALS_full[["ALSFRS_Total"]]
oDelta = dataALS_full[["onset_delta"]]

progressCnt = 0
earlySlope=c()
for(pat in patientSet){
	if( (progressCnt %% 500 ) == 0 ){
		cat("Processed: ",progressCnt, "\n")
	}
	progressCnt = progressCnt + 1;
	slope = NA
	if(!is.null(oDelta) && !is.null(aTotal)){
		sod = oDelta[oDelta$SubjectID==pat,]
		onset_delta = unlist(as.numeric(as.matrix(sod$feature_value)))[1]
		atot = aTotal[aTotal$SubjectID==pat,]
		atot$feature_delta = unlist(as.numeric(as.matrix(atot$feature_delta)))
		atot = atot[atot$feature_delta<92,]
		atot = atot[atot$feature_delta==max(atot$feature_delta),]
		most_recent_score = unlist(as.numeric(as.matrix(atot$feature_value)))[1]
		mrs_delta = unlist(as.numeric(as.matrix(atot$feature_delta)))[1]
		slope = (most_recent_score - 39) / (mrs_delta - onset_delta)
		slope = slope * 365/12
	}
	earlySlope = c(earlySlope, slope)
}


### create onset age
pAge = dataALS_full[["Age"]]
oDelta = dataALS_full[["onset_delta"]]

median_onset_delta = median( unlist(as.numeric(as.matrix( oDelta$feature_value ))) )
testAge = c()
onsetAge = c()
for(pat in patientSet){
	oAge = NA
	age = NA
	if(!is.null(oDelta) && !is.null(pAge)){
		sod = oDelta[oDelta$SubjectID==pat,]
		pag = pAge[pAge$SubjectID==pat,]
		onset_delta = unlist(as.numeric(as.matrix(sod$feature_value)))[1]
		age = unlist(as.numeric(as.matrix(pag$feature_value)))[1]
		if(!is.na(onset_delta) && !is.na(age)){
			oAge = age + onset_delta/365
		} else if(is.na(onset_delta) && !is.na(age)){
			oAge = age + median_onset_delta/365
		}
	}
	testAge=c(testAge,age)
	onsetAge = c(onsetAge, oAge)
}



categorical_pcures = c("MEDHx_Thyroid","Gender","MEDHx_Diabetes","MEDHx_Hypertens","family_ALS_hist","onset_site","Mutation")
numerical_pcures = setdiff(features_with_less_than_25perc_missing, categorical_pcures)
numerical_pcures = c(numerical_pcures, "priorSlope", "onsetAge")

sample_mode <- function(x){
	xv = as.vector(unlist(x))
	freqs = table(xv)
	smode = names(freqs)[freqs == max(freqs)]
	return(smode[1])
}

transformLongitudinalToSingular <- function(dALS, featName, patients, type, method=mean){
	featTable = dALS[[featName]]
	feature_singular = c()
	for(patient in patients){
		sub_tab = featTable[featTable$SubjectID==patient,]
		sub_tab = sub_tab[c("SubjectID","feature_value","feature_delta")]
		sub_tab = na.omit(sub_tab)
		
		patientIsMissing = FALSE
		if(is.null(sub_tab)){
			patientIsMissing = TRUE
		} else if(length(sub_tab[,1])<1) {
			patientIsMissing = TRUE
		}

		if(patientIsMissing){
			feature_singular=c(feature_singular, NA)
		} else {
			
			xy = data.frame(sub_tab$feature_delta, sub_tab$feature_value)
			colnames(xy)=c("delta","value")
			xy["delta"] = as.numeric(as.matrix(xy$delta))
			
			if(type=="numeric"){
				xy["value"] = as.numeric(as.matrix(xy$value))
				xy = na.omit(xy)
			} else {
				xy["value"] = as.character(as.matrix(xy$value))
				xy = na.omit(xy)
			}

			xy = xy[xy$delta<92,]
			if(length(xy[,1])<1){
				feature_singular=c(feature_singular, NA)
			} else{
				pfeat_sing = method(xy$value)
				feature_singular=c(feature_singular, pfeat_sing)
			}
		}
	}
	out = data.frame(patients, feature_singular, stringsAsFactors=FALSE)
	colnames(out)=c("SubjectID", featName)
	return(out)
}



### transform all relevant longitudinals into singular variables
out=data.frame(patientSet, earlySlope, onsetAge)
colnames(out) = c("SubjectID","priorSlope", "onsetAge")
for(i in features_with_less_than_25perc_missing){
	cat(i,"... ")
	if(i %in% categorical_pcures){
		tmp = transformLongitudinalToSingular(dALS=dataALS_full, featName=i, type="categorical", patients=patientSet, method=sample_mode)
		out = merge(out, tmp, by="SubjectID", all = TRUE)
	} else{
		tmp = NULL
		if(i == "ALSFRS_Total"){
			print("using min collapse!\n")
			tmp = transformLongitudinalToSingular(dALS=dataALS_full, featName=i, type="numeric", patients=patientSet, method=min)
		} else{
			tmp = transformLongitudinalToSingular(dALS=dataALS_full, featName=i, type="numeric", patients=patientSet, method=mean)
		}
		out = merge(out, tmp, by="SubjectID", all = TRUE)
	}
	cat("DONE!\n")
}
cat("\n")
write.csv(out, file=pathToSingularDataTable, row.names=FALSE)

out=read.csv(pathToSingularDataTable)



### impute singular values: KNN based

v_num_names = make.names(numerical_pcures)
v_cat_names = make.names(categorical_pcures)

slopeTab_tmp = data.frame(data.ALSslope$patients,data.ALSslope$ALSFRS_slope)
colnames(slopeTab_tmp) = c("SubjectID","ALSFRS_slope")
slopeTab_tmp = merge(slopeTab_tmp, out, by="SubjectID", all = TRUE) # THIS GOT CHANGED !!!

numTab = slopeTab_tmp[v_num_names]
numMat = as.matrix(numTab)
ki = impute.knn(data=numMat, k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)

kidat = ki$data
rightSkewed = c("Gamma.glutamyltransferase", "CK", "Red.Blood.Cells..RBC.", "Urine.Ph", "AST.SGOT.", "ALT.SGPT.", "Bilirubin..Total.")
leftSkewed = c("onset_delta", "hands", "leg", "mouth", "respiratory", "trunk", "ALSFRS_Total")
for(i in 1:length(kidat[1,])){
	currentFeature = kidat[,i]
	f = colnames(kidat)[i]
	if(f %in% rightSkewed){
		print("R")
		l_offs = abs(min(currentFeature)) + 1
		currentFeature = scale(log(currentFeature+l_offs))
	} else if(f %in% leftSkewed){
		print("L")
		currentFeature = scale(currentFeature^3)
	} else{
		currentFeature = scale(currentFeature)
	}
	kidat[,i] = currentFeature
}

catTab = slopeTab_tmp[v_cat_names]
for(i in colnames(catTab)){
	cur_col = unlist(catTab[i])
	cur_mode = sample_mode(na.omit(cur_col))
	cur_col[is.na(cur_col)] = cur_mode
	catTab[i] = as.factor(cur_col)
}

slopeTab_alt = data.frame(slopeTab_tmp$SubjectID, slopeTab_tmp$ALSFRS_slope, kidat, catTab)
colnames(slopeTab_alt)[1:2] = c("SubjectID","ALSFRS_slope")




### dummy coding

appendPrefix=function(i_names, prefix){
  if(length(i_names)==1){
	return( paste(prefix, i_names, sep="") )
  }
  newNames=c()
  for(i in 1:length(i_names)){
    featureName=paste(prefix, i_names[i], sep="")
    newNames=c(newNames, featureName)
  }
  return(newNames)
}

transformed_tab = slopeTab_alt[c("SubjectID", "ALSFRS_slope", v_num_names)]
v_dmy_names = c()
for(i in v_cat_names){
	col_vals = unlist(slopeTab_alt[i])
	modelString = paste("~", i, sep=" ")
	t_mat = model.matrix( eval(parse(text=modelString)), data = slopeTab_alt )
	t_mat = t_mat[,2:length(t_mat[1,])]
	lvls = levels(col_vals)
	lvls = lvls[2:length(lvls)]
	t_dat = data.frame(t_mat)
	colnames(t_dat) = appendPrefix(lvls, prefix=paste(i,".",sep=""))
	v_dmy_names = c(v_dmy_names, colnames(t_dat))
	transformed_tab = data.frame(transformed_tab, t_dat)
}
numerical_feature_names = v_num_names
categorical_feature_contrast_dummies = v_dmy_names


### prepare cleaned data files
data.surv<-read.delim(pathToSurvivalDataTable, sep="|", header=T)

transformed_tab = merge(transformed_tab, data.surv, by="SubjectID", all = TRUE)
save(transformed_tab, file=pathToDummyCodedFeatureTable)

feature_tab = slopeTab_alt[c("SubjectID", "ALSFRS_slope", v_num_names, v_cat_names)]
feature_tab = merge(feature_tab, data.surv, by="SubjectID", all = TRUE)
save(feature_tab, file=pathToNormallyCodedFeatureTable)

