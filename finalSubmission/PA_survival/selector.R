library(survival)

args<-commandArgs(trailingOnly = TRUE)      #to pass the selector.sh arguments to R
input_file_path<-args[1]
output_file_path<-args[2]

train_not_imputed = read.csv("/home/rama.raghavan/als_submission/ProActSurvival/singularDataTable_V05.csv")
#train_not_imputed = read.csv("R:/BIO RSCH/genomics/bfridley/ALS_Challenge/rpgm/richard/final_submission/singularDataTable_V05.csv")

# Read the input file
test_patient <-read.delim(input_file_path,sep="|",header=TRUE)
write.table(test_patient,file=paste("/home/rama.raghavan/als_submission/ProActSurvival/",basename(input_file_path),sep=""),sep="|")

test = test_patient[test_patient$form_name != "Adverse Event" & test_patient$form_name != "Concomitant Medication",]

feature_name = levels(factor(test_patient$feature_name))
dataALS_full = vector("list")
for(i in feature_name){
  temp = test[test$feature_name == i,]
  temp = temp[ , -which(names(temp) %in% c("form_name","feature_name"))]
  for(j in 2:length(temp)){
    temp[,j] = factor(temp[,j])
  }
  dataALS_full[[i]] = temp
}

patientSet = c(test_patient$SubjectID[1])

features_with_less_than_25perc_missing = c("Lymphocytes","Basophils","Monocytes","Total Cholesterol","Gamma-glutamyltransferase","CK","height","Red Blood Cells (RBC)","White Blood Cell (WBC)","Urine Ph","Bicarbonate","if_use_Riluzole","Q10_Respiratory","respiratory_rate","Calcium","Phosphorus","Platelets","Alkaline Phosphatase","bp_diastolic","bp_systolic","pulse","treatment_group","Hematocrit","Hemoglobin","Chloride","fvc","fvc_normal","fvc_percent","Potassium","Q5a_Cutting_without_Gastrostomy","Sodium","AST(SGOT)","Blood Urea Nitrogen (BUN)","Creatinine","ALT(SGPT)","Bilirubin (Total)","onset_delta","Age","ALSFRS_Total","Gender","hands","leg","mouth","onset_site","Q1_Speech","Q2_Salivation","Q3_Swallowing","Q4_Handwriting","Q5_Cutting","Q6_Dressing_and_Hygiene","Q7_Turning_in_Bed","Q8_Walking","Q9_Climbing_Stairs","Race","respiratory","trunk","weight")

### delete entries consisting of a single dash character '-'
for(i in features_with_less_than_25perc_missing){
	featTable = dataALS_full[[i]]
	selection = featTable$feature_value=="-"
	selection[is.na(selection)]=FALSE
	selection = !selection
	featTable = featTable[selection,]
	dataALS_full[[i]] = featTable
}


### create early ALSFRS_slope

aTotal = dataALS_full[["ALSFRS_Total"]]
oDelta = dataALS_full[["onset_delta"]]

earlySlope=c()
for(pat in patientSet){
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
oDelta_full = na.omit(train_not_imputed$onset_delta)

median_onset_delta = median( oDelta_full )
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



########################################################################################################

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
####################################################################################################################



### transform all relevant longitudinals into singular variables
categorical_pcures = c("Race", "onset_site", "Gender", "treatment_group", "if_use_Riluzole")
numerical_pcures = setdiff(features_with_less_than_25perc_missing, categorical_pcures)
numerical_pcures = c(numerical_pcures, "priorSlope", "onsetAge")

collapsed=data.frame(patientSet, earlySlope, onsetAge)

colnames(collapsed) = c("SubjectID", "priorSlope", "onsetAge")
for(i in features_with_less_than_25perc_missing){
	cat(i,"... ")
    tmp = NULL
	if(i %in% categorical_pcures){
		tmp = transformLongitudinalToSingular(dALS=dataALS_full, featName=i, type="categorical", patients=patientSet, method=sample_mode)
		collapsed = merge(collapsed, tmp, by="SubjectID", all = TRUE)
	} else{
		if(i == "ALSFRS_Total"){
			cat("using min collapse! ")
			tmp = transformLongitudinalToSingular(dALS=dataALS_full, featName=i, type="numeric", patients=patientSet, method=min)
		} else{
			tmp = transformLongitudinalToSingular(dALS=dataALS_full, featName=i, type="numeric", patients=patientSet, method=mean)
		}
		collapsed = merge(collapsed, tmp, by="SubjectID", all = TRUE)
	}
	cat("DONE!\n")
}



### impute singular values: KNN based

v_num_names = make.names(numerical_pcures)
v_cat_names = make.names(categorical_pcures)

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute)

colnames(collapsed) = make.names(colnames(collapsed))

collapsed$SubjectID = collapsed$SubjectID * (-1)

dat_merge = merge(train_not_imputed, collapsed, all = TRUE)
merged_numMat = as.matrix(dat_merge[v_num_names])

ki = impute.knn(data=merged_numMat, k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)
imputed_numerical_mat = ki$data

rightSkewed = c("Gamma.glutamyltransferase", "CK", "Red.Blood.Cells..RBC.", "Urine.Ph", "AST.SGOT.", "ALT.SGPT.", "Bilirubin..Total.")
leftSkewed = c("onset_delta", "hands", "leg", "mouth", "respiratory", "trunk", "ALSFRS_Total")
for(i in 1:length(imputed_numerical_mat[1,])){
	currentFeature = imputed_numerical_mat[,i]
	f = colnames(imputed_numerical_mat)[i]
	if(f %in% rightSkewed){
		l_offs = abs(min(currentFeature)) + 1
		currentFeature = scale(log(currentFeature+l_offs))
	} else if(f %in% leftSkewed){
		currentFeature = scale(currentFeature^3)
	} else{
		currentFeature = scale(currentFeature)
	}
	imputed_numerical_mat[,i] = currentFeature
}

row.names(imputed_numerical_mat) = dat_merge$SubjectID

catTab = dat_merge[v_cat_names]
for(i in colnames(catTab)){
	cur_col = unlist(catTab[i])
	cur_mode = sample_mode(na.omit(cur_col))
	cur_col[is.na(cur_col)] = cur_mode
	catTab[i] = as.factor(cur_col)
}
row.names(catTab) = dat_merge$SubjectID

patient_imputed = data.frame(row.names(imputed_numerical_mat), imputed_numerical_mat, catTab)
colnames(patient_imputed)[1] = c("SubjectID")

rownames(patient_imputed) = NULL
patient_imputed[nrow(patient_imputed)+1,] = patient_imputed[1,]
patient_imputed = patient_imputed[-1,]

patient_imputed = patient_imputed[nrow(patient_imputed),]


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

predictedGroupLabel = getGroupLabels(patient_imputed)

dat <-read.delim(input_file_path,sep="|",header=TRUE)


#(1) Surv(time = time_event, event = status) ~ fvc + priorSlope + onsetAge + fvc_normal + Chloride
#(2) Surv(time = time_event, event = status) ~ trunk + priorSlope + onsetAge + fvc_percent + Race + trunk*Race
#(3) Surv(time = time_event, event = status) ~ fvc + priorSlope + onsetAge + fvc_percent + Gender
#(4) Surv(time = time_event, event = status) ~ fvc + priorSlope + onsetAge + Hemoglobin + Chloride

# note that any model that contains priorSlope will pull onset_delta and ALSFRS_Total, since those two are used to calculate priorSlope
# note that any model that contains onsetAge will pull onset_delta and Age, since those two are used to calculate onsetAge
if(predictedGroupLabel == 1){
    feature_name = c("fvc","fvc_normal","Chloride","ALSFRS_Total","Age","onset_delta")
} else if(predictedGroupLabel == 2){
    feature_name = c("trunk","fvc_percent","Race","ALSFRS_Total","Age","onset_delta")
} else if(predictedGroupLabel == 3){
    feature_name = c("fvc","fvc_percent","Gender","ALSFRS_Total","Age","onset_delta")
} else if(predictedGroupLabel == 4){
    feature_name = c("fvc","Hemoglobin","Chloride","ALSFRS_Total","Age","onset_delta")
}



ind_cl<-predictedGroupLabel    

outdat = dat[dat$feature_name %in% feature_name,]

#write the output file
outdat<-rbind(c(paste("cluster: ",ind_cl,sep=""),rep("",dim(outdat)[2]-1)),as.matrix(outdat))
write.table(outdat,file=paste("/home/rama.raghavan/als_submission/ProActProgression/",basename(input_file_path),sep=""),sep="|",row.names=F,col.names=F,quote=F)
write.table(outdat,output_file_path,sep="|",row.names=F,col.names=F,quote=F)


