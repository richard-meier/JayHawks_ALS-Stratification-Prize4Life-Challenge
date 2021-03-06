args<-commandArgs(trailingOnly = TRUE)      #to pass the selector.sh arguments to R
input_file_path<-args[1]
output_file_path<-args[2]

train_not_imputed = read.csv("/home/rama.raghavan/als_submission/NatRegProgression/singularDataTable_reg.csv")
load("/home/rama.raghavan/als_submission/NatRegProgression/clusterModels_ANR.RData")

# Read the input file
test_patient <-read.delim(input_file_path,sep="|",header=TRUE)

write.table(test_patient,file=paste("/home/rama.raghavan/als_submission/NatRegProgression/z_",basename(input_file_path),sep=""),sep="|",row.names=TRUE)

cat("Cluster value =================>",names(test_patient)[1],"\n")
cl_name_header = as.character(names(test_patient)[1])
cluster = as.numeric(substr(cl_name_header, 10, nchar(cl_name_header)))
cat(cluster,"\n")

colnames(test_patient) <- c("SubjectID","form_name","feature_name","feature_value","feature_unit","feature_delta")

test_patient$SubjectID = unlist(as.character(as.matrix(test_patient$SubjectID)))

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

features_with_less_than_25perc_missing = c("Age","ALS_Staging_Breath","ALS_Staging_Com","ALS_Staging_Eat","ALS_Staging_Move","ALS_Staging_Total","ALSFRS_R_Total",
 "ALSFRS_Total","family_ALS_hist","Gender","hands","leg","MEDHx_Diabetes","MEDHx_Hypertens","MEDHx_Thyroid","mouth",
 "Mutation","onset_delta","onset_site","Q1_Speech","Q2_Salivation","Q3_Swallowing","Q4_Handwriting","Q5_Cutting",
 "Q6_Dressing_and_Hygiene","Q7_Turning_in_Bed","Q8_Walking","Q9_Climbing_Stairs","R1_Dyspnea","R2_Orthopnea",
 "R3_Respiratory_Insufficiency","respiratory","respiratory_R","trunk")
 
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
categorical_pcures = c("MEDHx_Thyroid","Gender","MEDHx_Diabetes","MEDHx_Hypertens","family_ALS_hist","onset_site","Mutation")

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
collapsed$SubjectID = paste("!!!",collapsed$SubjectID,sep="")
train_not_imputed$SubjectID = unlist(as.character(as.matrix(train_not_imputed$SubjectID)))

dat_merge = merge(train_not_imputed, collapsed, all = TRUE, sort=FALSE)
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
patient_imputed = patient_imputed[nrow(patient_imputed),]


timePoints = c(12, 18, 24) * 30.5
if(cluster == 1){
	result = predict(modelForCluster_1, newdata = patient_imputed, se.fit=TRUE)
	avgSlp = avgSlp1
} else if(cluster == 2){
    result = predict(modelForCluster_2, newdata = patient_imputed, se.fit=TRUE)
	avgSlp = avgSlp2
} else if(cluster == 3){
    result = predict(modelForCluster_3, newdata = patient_imputed, se.fit=TRUE)
	avgSlp = avgSlp3
}

res_out=c(result$fit[1], result$se.fit[1])
if(is.na(res_out[1])){
	res_out[1] = avgSlp
	res_out[2] = 0
	cat("WARNING! Slope is NA and is replaced!!!")
}

outdat<-c(paste(result$fit,result$se.fit,sep="|"))

write.table(outdat,output_file_path,sep="|",row.names=F,col.names=F,quote=F)
