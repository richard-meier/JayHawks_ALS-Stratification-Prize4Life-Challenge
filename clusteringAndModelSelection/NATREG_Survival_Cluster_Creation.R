library(rpart)
library(OIsurv)
library(partykit)
load(".../slope_features_default_reg.rdata") 
names(feature_tab)

par(mfrow=c(2,1))
tfit<-rpart(
	Surv(time_event, status) ~ Age + ALS_Staging_Breath + ALS_Staging_Com + ALS_Staging_Eat + ALS_Staging_Move + ALS_Staging_Total + 
	ALSFRS_R_Total + ALSFRS_Total + hands + leg + mouth + onset_delta + Q1_Speech + Q2_Salivation + Q3_Swallowing + Q4_Handwriting + 
	Q5_Cutting + Q6_Dressing_and_Hygiene + Q7_Turning_in_Bed + Q8_Walking + Q9_Climbing_Stairs + R1_Dyspnea + R2_Orthopnea + 
	R3_Respiratory_Insufficiency + respiratory + respiratory_R + trunk + priorSlope + onsetAge + MEDHx_Thyroid + Gender + 
	MEDHx_Diabetes + MEDHx_Hypertens + family_ALS_hist + onset_site + Mutation, 
	control=rpart.control(minsplit=400, cp=0.005), data=feature_tab
)
plot(tfit)
text(tfit)

tfit2<-as.party(tfit)
plot(tfit2)
