library(rpart)
load(".../slope_features_default_reg.rdata")
names(feature_tab)

fit1<-rpart(
	ALSFRS_slope ~ ALS_Staging_Breath + ALS_Staging_Com + ALS_Staging_Eat + ALS_Staging_Move + ALS_Staging_Total + ALSFRS_R_Total + 
	ALSFRS_Total + family_ALS_hist + Gender + hands + leg + MEDHx_Diabetes + MEDHx_Hypertens + MEDHx_Thyroid + mouth + Mutation + 
	onset_delta + onset_site + Q1_Speech + Q2_Salivation + Q3_Swallowing + Q7_Turning_in_Bed + Q8_Walking + Q9_Climbing_Stairs +           
    Q4_Handwriting + Q5_Cutting + Q6_Dressing_and_Hygiene + R1_Dyspnea + R2_Orthopnea + R3_Respiratory_Insufficiency + respiratory + 
	respiratory_R + trunk + priorSlope + onsetAge, 
	method="anova", control=rpart.control(minsplit=25, cp=0.005),data=regslope2
)
plot(fit1, uniform=TRUE, main="Classification Tree for ALSFRS_slope")
text(fit1, use.n=TRUE, all=TRUE, cex=.8)
