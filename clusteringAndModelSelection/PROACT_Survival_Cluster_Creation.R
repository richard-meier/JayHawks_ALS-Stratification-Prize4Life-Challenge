library(rpart)
library(OIsurv)
library(partykit)
load(".../slope_features_default_V05.rdata") 
names(feature_tab)

par(mfrow=c(2,1))
tfit<-rpart(
	Surv(time_event, status) ~ Lymphocytes + Basophils + Monocytes + Total.Cholesterol + Gamma.glutamyltransferase + CK +
	height + Red.Blood.Cells..RBC. + White.Blood.Cell..WBC. + Urine.Ph + Bicarbonate + Q10_Respiratory + respiratory_rate + 
	Calcium + Phosphorus + Platelets + Alkaline.Phosphatase + bp_diastolic + bp_systolic + pulse + Hematocrit + Hemoglobin + 
	Chloride + fvc + fvc_normal + fvc_percent + Potassium + Q5a_Cutting_without_Gastrostomy + Sodium + AST.SGOT. + 
	Blood.Urea.Nitrogen..BUN. + Creatinine + ALT.SGPT. + Bilirubin..Total. + onset_delta + onsetAge + Age + ALSFRS_Total + 
	hands + leg + mouth + Q1_Speech + Q2_Salivation + Q3_Swallowing + Q4_Handwriting + Q5_Cutting + Q6_Dressing_and_Hygiene + 
	Q7_Turning_in_Bed + Q8_Walking + Q9_Climbing_Stairs + respiratory + trunk + weight + Race + onset_site + Gender + 
	treatment_group + if_use_Riluzole, 
	control=rpart.control(minsplit=3000, cp=0.005), data=feature_tab
)
plot(tfit)
text(tfit)

tfit2<-as.party(tfit)
plot(tfit2)
