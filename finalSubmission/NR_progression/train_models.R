
load(".../slope_features_default_reg.rdata") 
outFile = ".../clusterModels_ANR.RData"

# hardcode classification

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

timePoints = c(12, 18, 24) * 30.5

#(1) ALSFRS_slope ~ ALS_Staging_Com + ALSFRS_Total + Q1_Speech + Q9_Climbing_Stairs + hands
#(2) ALSFRS_slope ~ Q1_Speech + Gender + Q4_Handwriting + MEDHx_Diabetes + ALS_Staging_Total + ALS_Staging_Move
#(3) ALSFRS_slope ~ onset_site+ALS_Staging_Eat+priorSlope+Q6_Dressing_and_Hygiene+trunk + Gender
 
grp_tab1 = train_tab[train_tab$CART_label==1,]
modelForCluster_1 = lm(ALSFRS_slope ~ ALS_Staging_Com + ALSFRS_Total + Q1_Speech + Q9_Climbing_Stairs + hands, data=grp_tab1)
slp = predict(modelForCluster_1, newdata = grp_tab1, interval="predict")
avgSlp1 = mean(slp[,1])

grp_tab2 = train_tab[train_tab$CART_label==2,]
modelForCluster_2 = lm(ALSFRS_slope ~ Q1_Speech + Gender + Q4_Handwriting + MEDHx_Diabetes + ALS_Staging_Total + ALS_Staging_Move, data=grp_tab2)
slp = predict(modelForCluster_2, newdata = grp_tab2, interval="predict")
avgSlp2 = mean(slp[,1])

grp_tab3 = train_tab[train_tab$CART_label==3,]
modelForCluster_3 = lm(ALSFRS_slope ~ onset_site + ALS_Staging_Eat + priorSlope + Q6_Dressing_and_Hygiene + trunk + Gender, data=grp_tab3)
slp = predict(modelForCluster_3, newdata = grp_tab3, interval="predict")
avgSlp3 = mean(slp[,1])

save(
	modelForCluster_1, modelForCluster_2, modelForCluster_3,
	grp_tab1, grp_tab2, grp_tab3,
	avgSlp1, avgSlp2, avgSlp3,
	file=outFile
)
