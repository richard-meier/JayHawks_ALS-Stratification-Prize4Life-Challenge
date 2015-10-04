#library(survival)
#library(pec)

load(".../slope_features_default_V05.rdata") 
outFile = ".../clusterModels_APA.RData"

# hardcoded classification

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

timePoints = c(12, 18, 24) * 30.5

#(111)	ALSFRS_slope ~ Q9_Climbing_Stairs + treatment_group + Phosphorus
#(112)	ALSFRS_slope ~ Gender + Creatinine + Alkaline.Phosphatase + White.Blood.Cell..WBC. + Creatinine*White.Blood.Cell..WBC. + Potassium + White.Blood.Cell..WBC.*Potassium
#(12)	ALSFRS_slope ~ ALSFRS_Total + CK + ALSFRS_Total*CK + fvc_normal
#(21)	ALSFRS_slope ~ onset_delta + bp_systolic + Q10_Respiratory + Alkaline.Phosphatase
#(22)	ALSFRS_slope ~ onset_delta + Q5_Cutting + fvc_percent
#(3)	ALSFRS_slope ~ mouth + trunk + fvc + mouth*trunk
#(4)	ALSFRS_slope ~ trunk + Race + bp_diastolic + if_use_Riluzole

grp_tab111 = train_tab[train_tab$CART_label==111,]
modelForCluster_111 = lm(ALSFRS_slope ~ Q9_Climbing_Stairs + treatment_group + Phosphorus, data=grp_tab111)
slp = predict(modelForCluster_111, newdata = grp_tab111, interval="predict")
avgSlp111 = mean(slp[,1])

grp_tab112 = train_tab[train_tab$CART_label==112,]
modelForCluster_112 = lm(ALSFRS_slope ~ Gender + Creatinine + Alkaline.Phosphatase + White.Blood.Cell..WBC. + Creatinine*White.Blood.Cell..WBC. + Potassium + White.Blood.Cell..WBC.*Potassium, data=grp_tab112)
slp = predict(modelForCluster_112, newdata = grp_tab112, interval="predict")
avgSlp112 = mean(slp[,1])

grp_tab12 = train_tab[train_tab$CART_label==12,]
modelForCluster_12 = lm(ALSFRS_slope ~ ALSFRS_Total + CK + ALSFRS_Total*CK + fvc_normal, data=grp_tab12)
slp = predict(modelForCluster_12, newdata = grp_tab12, interval="predict")
avgSlp12 = mean(slp[,1])

grp_tab21 = train_tab[train_tab$CART_label==21,]
modelForCluster_21 = lm(ALSFRS_slope ~ onset_delta + bp_systolic + Q10_Respiratory + Alkaline.Phosphatase, data=grp_tab21)
slp = predict(modelForCluster_21, newdata = grp_tab21, interval="predict")
avgSlp21 = mean(slp[,1])

grp_tab22 = train_tab[train_tab$CART_label==22,]
modelForCluster_22 = lm(ALSFRS_slope ~ onset_delta + Q5_Cutting + fvc_percent, data=grp_tab22)
slp = predict(modelForCluster_22, newdata = grp_tab22, interval="predict")
avgSlp22 = mean(slp[,1])

grp_tab3 = train_tab[train_tab$CART_label==3,]
modelForCluster_3 = lm(ALSFRS_slope ~ mouth + trunk + fvc + mouth*trunk, data=grp_tab3)
slp = predict(modelForCluster_3, newdata = grp_tab3, interval="predict")
avgSlp3 = mean(slp[,1])

grp_tab4 = train_tab[train_tab$CART_label==4,]
white_sub = grp_tab4[grp_tab4$Race == "White",]
dummyPatient = white_sub[1,]
dummyPatient$Race = "American Indian"
grp_tab4 = rbind(grp_tab4, dummyPatient)
modelForCluster_4 = lm(ALSFRS_slope ~ trunk + Race + bp_diastolic + if_use_Riluzole, data=grp_tab4)
slp = predict(modelForCluster_4, newdata = grp_tab4, interval="predict")
avgSlp4 = mean(slp[,1])

save(
	modelForCluster_111, modelForCluster_112, modelForCluster_12, modelForCluster_21, modelForCluster_22, modelForCluster_3, modelForCluster_4,
	grp_tab111, grp_tab112, grp_tab12, grp_tab21, grp_tab22, grp_tab3, grp_tab4,
	avgSlp111, avgSlp112, avgSlp12, avgSlp21, avgSlp22, avgSlp3, avgSlp4,
	file=outFile
)
