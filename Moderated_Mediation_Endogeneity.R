
#################################################
# IMPORT REQUIRED PACKAGES
library(psych)
library(GPArotation)
library(lavaan)
library(Hmisc)
library(polycor)
library(QuantPsyc)
library(car)
library(lmSupport)
library(pwr)
library(AER)
library(Hmisc)

#################################################
# IMPORT DATASET
## first the full data set
ITSMdatafull <- read.csv("C:/Users/jwulf/OneDrive/Research/IWI Forschung/JMIS/Userinputs_v13_n256.csv", na.strings=c("-77"," "))

##prepare individual instruments
ITSMperformance <- ITSMdatafull[,c('Perform1', 'Perform2', 'Perform3', 	'Perform4', 'Perform5')]
##REplace NAs with means
for(i in 1:ncol(ITSMperformance)){
  ITSMperformance[is.na(ITSMperformance[,i]), i] <- mean(ITSMperformance[,i], na.rm = TRUE)
}

ITSMalign <- ITSMdatafull[,c('Align1',	'Align2',	'Align3','Align4' ,	'Align5' ,	'Align6')]
for(i in 1:ncol(ITSMalign)){
  ITSMalign[is.na(ITSMalign[,i]), i] <- mean(ITSMalign[,i], na.rm = TRUE)
}


ITSMdata1 <- ITSMdatafull[,c('Ind1_ProdVsService',	'Ind2_PhysVsInform')]
for(i in 1:ncol(ITSMdata1)){
  ITSMdata1[is.na(ITSMdata1[,i]), i] <- mean(ITSMdata1[,i], na.rm = TRUE)
}

ITSMdata2 <-  ITSMdatafull[,c('IndReg1', 'IndReg2')]
for(i in 1:ncol(ITSMdata2)){
  ITSMdata2[is.na(ITSMdata2[,i]), i] <- mean(ITSMdata2[,i], na.rm = TRUE)
}

ITSMdata7 <- ITSMdatafull[,c('ITStrat_Inno1','ITStrat_Inno2','ITStrat_Inno3')]
for(i in 1:ncol(ITSMdata7)){
  ITSMdata7[is.na(ITSMdata7[,i]), i] <- mean(ITSMdata7[,i], na.rm = TRUE)
}

ITSMdata8 <- ITSMdatafull[,c('ITStrat_Cons2','ITStrat_Cons3')]
for(i in 1:ncol(ITSMdata8)){
  ITSMdata8[is.na(ITSMdata8[,i]), i] <- mean(ITSMdata8[,i], na.rm = TRUE)
}

# ITSMdata9 <- ITSMdatafull[,c('AttitudITIL1', 'AttitudITIL2', 'AttitudITIL3')]
# for(i in 1:ncol(ITSMdata9)){
#   ITSMdata9[is.na(ITSMdata9[,i]), i] <- mean(ITSMdata9[,i], na.rm = TRUE)
# }
##prepare EFA

ITSMdataEFA <- ITSMdatafull[,c('ITIL_SS3_FinMgt', 'ITIL_SS4_DmdMgt', 'ITIL_SS5_BRelMgt', 'ITIL_SD3_AvMgt', 'ITIL_SD4_CapMgt', 'ITIL_SD5_ContMgt', 'ITIL_ST2_ChangeM', 'ITIL_ST3_ConfMgt', 'ITIL_ST4_RelMgt', 'ITIL_ST5_SValTestM', 'ITIL_ST6_ChngEval', 'ITIL_SO2_IncMgt', 'ITIL_SO3_ReqFul', 'ITIL_SO4_ProbMgt')]

###append Align
ITSMdataEFA <- cbind(ITSMdataEFA, ITSMdatafull[,c('Align1',	'Align2',	'Align3','Align4' ,	'Align5' ,	'Align6')])

###append Performance
ITSMdataEFA <- cbind(ITSMdataEFA, ITSMdatafull[,c('Perform1', 'Perform2', 'Perform3', 	'Perform4', 'Perform5')])

###append IndustryType
ITSMdataEFA <- cbind(ITSMdataEFA, ITSMdatafull[,c('Ind1_ProdVsService',	'Ind2_PhysVsInform')])

###append Regulation
ITSMdataEFA <- cbind(ITSMdataEFA, ITSMdatafull[,c('IndReg1', 'IndReg2')])

###append Inno
ITSMdataEFA <- cbind(ITSMdataEFA, ITSMdatafull[,c('ITStrat_Inno1','ITStrat_Inno2','ITStrat_Inno3')])

###append Cons
ITSMdataEFA <- cbind(ITSMdataEFA, ITSMdatafull[,c('ITStrat_Cons2','ITStrat_Cons3')])

##REplace NAs with means
for(i in 1:ncol(ITSMdataEFA)){
  ITSMdataEFA[is.na(ITSMdataEFA[,i]), i] <- mean(ITSMdataEFA[,i], na.rm = TRUE)
}

##Now the EFA factor analysis
corMatEFA <- cor(ITSMdataEFA)
EFAfull_results <- fa(r = corMatEFA, nfactors = 9, rotate = "bentlerQ", fm = "minres", scores = "regression")
View(unclass(EFAfull_results$loadings))
View(unclass(EFAfull_results$Structure))
View(t(EFAfull_results$values))
EFAfull_results_tests <- VSS(corMatEFA, n = 9, rotate = "bentlerQ", fm = "ml", n.obs=256)
nfactors(corMatEFA, n = 9, rotate = "bentlerQ", fm = "ml", n.obs=256)



###################################################
# CALCULATE INTERNAL CONSISTENCIES
###################################################


EFAfull_results_scoring <- factor.scores(ITSMdataEFA,EFAfull_results)
EFAfull_results_scores <- EFAfull_results_scoring$scores
alphaITSM <- psych::alpha(EFAfull_results_scores[,c('MR2','MR8','MR9')])
alphaITSM$total$std.alpha

alphaPlan <- psych::alpha(ITSMdataEFA[,c('ITIL_SS3_FinMgt', 'ITIL_SS4_DmdMgt', 'ITIL_SS5_BRelMgt', 'ITIL_SD3_AvMgt', 'ITIL_SD4_CapMgt', 'ITIL_SD5_ContMgt')])
alphaPlan$total$std.alpha

alphaTransit <- psych::alpha(ITSMdataEFA[,c('ITIL_ST4_RelMgt', 'ITIL_ST5_SValTestM', 'ITIL_ST6_ChngEval')])
alphaTransit$total$std.alpha

alphaOperate <- psych::alpha(ITSMdataEFA[,c('ITIL_ST2_ChangeM', 'ITIL_ST3_ConfMgt','ITIL_SO2_IncMgt', 'ITIL_SO3_ReqFul', 'ITIL_SO4_ProbMgt')])
alphaOperate$total$std.alpha

alphaAlign <- psych::alpha(ITSMdataEFA[,c('Align1','Align2','Align3','Align4' ,'Align5' ,'Align6')])
alphaAlign$total$std.alpha

alphaPerform <- psych::alpha(ITSMdataEFA[,c('Perform1', 'Perform2', 'Perform3', 	'Perform4', 'Perform5')])
alphaPerform$total$std.alpha

alphaInd <- psych::alpha(ITSMdataEFA[,c('Ind1_ProdVsService',	'Ind2_PhysVsInform')])
alphaInd$total$std.alpha

alphaReg <- psych::alpha(ITSMdataEFA[,c('IndReg1','IndReg2')])
alphaReg$total$std.alpha

alphaStratInno <- psych::alpha(ITSMdataEFA[,c('ITStrat_Inno1','ITStrat_Inno2','ITStrat_Inno3')])
alphaStratInno$total$std.alpha

alphaStratCons <- psych::alpha(ITSMdataEFA[,c('ITStrat_Cons2','ITStrat_Cons3')])
alphaStratCons$total$std.alpha

###################################################
# Assign Scores and calculate descriptive statistics
###################################################

# mean and sd over all included ITSM pratices
mean(as.matrix(ITSMdatafull[,c('ITIL_SS3_FinMgt', 'ITIL_SS4_DmdMgt', 'ITIL_SS5_BRelMgt', 'ITIL_SD3_AvMgt', 'ITIL_SD4_CapMgt', 'ITIL_SD5_ContMgt', 'ITIL_ST2_ChangeM', 'ITIL_ST3_ConfMgt', 'ITIL_ST4_RelMgt', 'ITIL_ST5_SValTestM', 'ITIL_ST6_ChngEval', 'ITIL_SO2_IncMgt', 'ITIL_SO3_ReqFul', 'ITIL_SO4_ProbMgt')]))
sd(as.matrix(ITSMdatafull[,c('ITIL_SS3_FinMgt', 'ITIL_SS4_DmdMgt', 'ITIL_SS5_BRelMgt', 'ITIL_SD3_AvMgt', 'ITIL_SD4_CapMgt', 'ITIL_SD5_ContMgt', 'ITIL_ST2_ChangeM', 'ITIL_ST3_ConfMgt', 'ITIL_ST4_RelMgt', 'ITIL_ST5_SValTestM', 'ITIL_ST6_ChngEval', 'ITIL_SO2_IncMgt', 'ITIL_SO3_ReqFul', 'ITIL_SO4_ProbMgt')]))

#SP
mean(as.matrix(ITSMdatafull[,c('ITIL_SS3_FinMgt', 'ITIL_SS4_DmdMgt', 'ITIL_SS5_BRelMgt', 'ITIL_SD3_AvMgt', 'ITIL_SD4_CapMgt', 'ITIL_SD5_ContMgt')]))
sd(as.matrix(ITSMdatafull[,c('ITIL_SS3_FinMgt', 'ITIL_SS4_DmdMgt', 'ITIL_SS5_BRelMgt', 'ITIL_SD3_AvMgt', 'ITIL_SD4_CapMgt', 'ITIL_SD5_ContMgt')]))

#ST
mean(as.matrix(ITSMdatafull[,c('ITIL_ST4_RelMgt', 'ITIL_ST5_SValTestM', 'ITIL_ST6_ChngEval')]))
sd(as.matrix(ITSMdatafull[,c('ITIL_ST4_RelMgt', 'ITIL_ST5_SValTestM', 'ITIL_ST6_ChngEval')]))

#SO
mean(as.matrix(ITSMdatafull[,c('ITIL_ST2_ChangeM', 'ITIL_ST3_ConfMgt', 'ITIL_SO2_IncMgt', 'ITIL_SO3_ReqFul', 'ITIL_SO4_ProbMgt')]))
sd(as.matrix(ITSMdatafull[,c('ITIL_ST2_ChangeM', 'ITIL_ST3_ConfMgt', 'ITIL_SO2_IncMgt', 'ITIL_SO3_ReqFul', 'ITIL_SO4_ProbMgt')]))

mean(as.matrix(ITSMdataEFA[,c('Align1','Align2','Align3','Align4' ,'Align5' ,'Align6')]))
sd(as.matrix(ITSMdataEFA[,c('Align1','Align2','Align3','Align4' ,'Align5' ,'Align6')]))

mean(as.matrix(ITSMdataEFA[,c('Perform1', 'Perform2', 'Perform3', 	'Perform4', 'Perform5')]))
sd(as.matrix(ITSMdataEFA[,c('Perform1', 'Perform2', 'Perform3', 	'Perform4', 'Perform5')]))

mean(as.matrix(ITSMdataEFA[,c('Ind1_ProdVsService',	'Ind2_PhysVsInform')]))
sd(as.matrix(ITSMdataEFA[,c('Ind1_ProdVsService',	'Ind2_PhysVsInform')]))

mean(as.matrix(ITSMdataEFA[,c('IndReg1','IndReg2')]))
sd(as.matrix(ITSMdataEFA[,c('IndReg1','IndReg2')]))

mean(as.matrix(ITSMdataEFA[,c('ITStrat_Inno1','ITStrat_Inno2','ITStrat_Inno3')]))
sd(as.matrix(ITSMdataEFA[,c('ITStrat_Inno1','ITStrat_Inno2','ITStrat_Inno3')]))

mean(as.matrix(ITSMdataEFA[,c('ITStrat_Cons2','ITStrat_Cons3')]))
sd(as.matrix(ITSMdataEFA[,c('ITStrat_Cons2','ITStrat_Cons3')]))


###################################################
# CALCULATE FACTOR SCORES FOR ITSMcap
## First the factor scores for all instruments

##Second order ITSM construct
####ADAPT M'S!!
corMatITSM <- cor(EFAfull_results_scores[,c('MR2','MR8','MR9')])
efa_results_3itsm <- fa(r = corMatITSM, nfactors = 1, rotate = "bentlerQ", fm = "minres")
View(unclass(efa_results_3itsm$loadings))
EFAfull_results_scoring_itsm <- factor.scores(EFAfull_results_scores[,c('MR2','MR8','MR9')],efa_results_3itsm)
ITSMdatafull$ITSMCapScore <- EFAfull_results_scoring_itsm$scores

###Scores for performance
corMatp <- cor(ITSMperformance)
efa_resultsp <- fa(r = corMatp, nfactors = 1, rotate = "BentlerQ", fm = "minres")
scoresp <- factor.scores(ITSMperformance, efa_resultsp)
ITSMdatafull$PerformanceScore <- scoresp$scores[,1]

###Scores for alignment
corMata <- cor(ITSMalign)
efa_resultsa <- fa(r = corMata, nfactors = 1, rotate = "BentlerQ", fm = "minres")
scoresa <- factor.scores(ITSMalign, efa_resultsa)
ITSMdatafull$AlignScore <- scoresa$scores[,1]

###Scores for the others
corMat1 <- cor(ITSMdata1)
efa_results1 <- fa(r = corMat1, nfactors = 1, rotate = "BentlerQ", fm = "minres")
scores1 <- factor.scores(ITSMdata1, efa_results1)
ITSMdatafull$IndScore <- scores1$scores[,1]

corMat2 <- cor(ITSMdata2)
efa_results2 <- fa(r = corMat2, nfactors = 1, rotate = "BentlerQ", fm = "minres")
scores2 <- factor.scores(ITSMdata2, efa_results2)
ITSMdatafull$RegScore <- scores2$scores[,1]

corMat7 <- cor(ITSMdata7)
efa_results7 <- fa(r = corMat7, nfactors = 1, rotate = "BentlerQ", fm = "minres")
scores7 <- factor.scores(ITSMdata7, efa_results7)
ITSMdatafull$ITStrat_Inno <- scores7$scores[,1]

corMat8 <- cor(ITSMdata8)
efa_results8 <- fa(r = corMat8, nfactors = 1, rotate = "BentlerQ", fm = "minres")
scores8 <- factor.scores(ITSMdata8, efa_results8)
ITSMdatafull$ITStrat_Cons <- scores8$scores[,1]

##Build StratOrient
ITSMdatafull$ITStratOrient <- ITSMdatafull$ITStrat_Cons - ITSMdatafull$ITStrat_Inno


###and here the correlation coefficients
corresults <- cor(ITSMdatafull[,c('ITSMCapScore', 'PerformanceScore', 'AlignScore','ITStratOrient', 'Log_Empl_Client', 'IndScore', 'RegScore', 'PositionBusSP01', 'Position_Vertical', 'PositionTime')], use="pairwise.complete.obs")
View(corresults)

corresults2 <- rcorr(as.matrix(ITSMdatafull[,c('ITSMCapScore', 'PerformanceScore', 'AlignScore','ITStratOrient', 'Log_Empl_Client', 'IndScore', 'RegScore', 'PositionBusSP01', 'Position_Vertical', 'PositionTime')]))
View(corresults2$P)

###build the application indicators: xxx0 is existence, xxx1 is ITIL coverage
ITSMdatafull$SPf1_ProjMgmt0 <- ifelse(ITSMdatafull$SPf1_ProjMgmt==1, 0, 1)
ITSMdatafull$SPf1_ProjMgmt1 <- ifelse((ITSMdatafull$SPf1_ProjMgmt==1) | (ITSMdatafull$SPf1_ProjMgmt==2), 0, 1)
ITSMdatafull$SPf1_ProjMgmt0 <- as.factor(ITSMdatafull$SPf1_ProjMgmt0)
ITSMdatafull$SPf1_ProjMgmt1 <- as.factor(ITSMdatafull$SPf1_ProjMgmt1)
ITSMdatafull$SPf2_SoftwDev0 <- ifelse(ITSMdatafull$SPf2_SoftwDev==1, 0, 1)
ITSMdatafull$SPf2_SoftwDev1 <- ifelse((ITSMdatafull$SPf2_SoftwDev==1) | (ITSMdatafull$SPf2_SoftwDev==2), 0, 1)
ITSMdatafull$SPf2_SoftwDev0 <- as.factor(ITSMdatafull$SPf2_SoftwDev0)
ITSMdatafull$SPf2_SoftwDev1 <- as.factor(ITSMdatafull$SPf2_SoftwDev1)
ITSMdatafull$SPf3_ApplCust0 <- ifelse(ITSMdatafull$SPf3_ApplCust==1, 0, 1)
ITSMdatafull$SPf3_ApplCust1 <- ifelse((ITSMdatafull$SPf3_ApplCust==1) | (ITSMdatafull$SPf3_ApplCust==2), 0, 1)
ITSMdatafull$SPf3_ApplCust0 <- as.factor(ITSMdatafull$SPf3_ApplCust0)
ITSMdatafull$SPf3_ApplCust1 <- as.factor(ITSMdatafull$SPf3_ApplCust1)
ITSMdatafull$SPf4_ApplMgmt0 <- ifelse(ITSMdatafull$SPf4_ApplMgmt==1, 0, 1)
ITSMdatafull$SPf4_ApplMgmt1 <- ifelse((ITSMdatafull$SPf4_ApplMgmt==1) | (ITSMdatafull$SPf4_ApplMgmt==2), 0, 1)
ITSMdatafull$SPf4_ApplMgmt0 <- as.factor(ITSMdatafull$SPf4_ApplMgmt0)
ITSMdatafull$SPf4_ApplMgmt1 <- as.factor(ITSMdatafull$SPf4_ApplMgmt1)
ITSMdatafull$SPf5_DataCentr0 <- ifelse(ITSMdatafull$SPf5_DataCentr==1, 0, 1)
ITSMdatafull$SPf5_DataCentr1 <- ifelse((ITSMdatafull$SPf5_DataCentr==1) | (ITSMdatafull$SPf5_DataCentr==2), 0, 1)
ITSMdatafull$SPf5_DataCentr0 <- as.factor(ITSMdatafull$SPf5_DataCentr0)
ITSMdatafull$SPf5_DataCentr1 <- as.factor(ITSMdatafull$SPf5_DataCentr1)
ITSMdatafull$SPf6_NetwCom0 <- ifelse(ITSMdatafull$SPf6_NetwCom==1, 0, 1)
ITSMdatafull$SPf6_NetwCom1 <- ifelse((ITSMdatafull$SPf6_NetwCom==1) | (ITSMdatafull$SPf6_NetwCom==2), 0, 1)
ITSMdatafull$SPf6_NetwCom0 <- as.factor(ITSMdatafull$SPf6_NetwCom0)
ITSMdatafull$SPf6_NetwCom1 <- as.factor(ITSMdatafull$SPf6_NetwCom1)
ITSMdatafull$SPf7_DesktPrint0 <- ifelse(ITSMdatafull$SPf7_DesktPrint==1, 0, 1)
ITSMdatafull$SPf7_DesktPrint1 <- ifelse((ITSMdatafull$SPf7_DesktPrint==1) | (ITSMdatafull$SPf7_DesktPrint==2), 0, 1)
ITSMdatafull$SPf7_DesktPrint0 <- as.factor(ITSMdatafull$SPf7_DesktPrint0)
ITSMdatafull$SPf7_DesktPrint1 <- as.factor(ITSMdatafull$SPf7_DesktPrint1)
ITSMdatafull$SPf8_SupHlpDsk0 <- ifelse(ITSMdatafull$SPf8_SupHlpDsk==1, 0, 1)
ITSMdatafull$SPf8_SupHlpDsk1 <- ifelse((ITSMdatafull$SPf8_SupHlpDsk==1) | (ITSMdatafull$SPf8_SupHlpDsk==2), 0, 1)
ITSMdatafull$SPf8_SupHlpDsk0 <- as.factor(ITSMdatafull$SPf8_SupHlpDsk0)
ITSMdatafull$SPf8_SupHlpDsk1 <- as.factor(ITSMdatafull$SPf8_SupHlpDsk1)
ITSMdatafull$SPf9_Training0 <- ifelse(ITSMdatafull$SPf9_Training==1, 0, 1)
ITSMdatafull$SPf9_Training1 <- ifelse((ITSMdatafull$SPf9_Training==1) | (ITSMdatafull$SPf9_Training==2), 0, 1)
ITSMdatafull$SPf9_Training0 <- as.factor(ITSMdatafull$SPf9_Training0)
ITSMdatafull$SPf9_Training1 <- as.factor(ITSMdatafull$SPf9_Training1)
ITSMdatafull$SPf10_Other0 <- ifelse(ITSMdatafull$SPf10_Other==1, 0, 1)
ITSMdatafull$SPf10_Other1 <- ifelse((ITSMdatafull$SPf10_Other==1) | (ITSMdatafull$SPf10_Other==2), 0, 1)
ITSMdatafull$SPf10_Other0 <- as.factor(ITSMdatafull$SPf10_Other0)
ITSMdatafull$SPf10_Other1 <- as.factor(ITSMdatafull$SPf10_Other1)




###fill missing data
for(i in 1:ncol(ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')])){
  #convert to number
  ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][,i] <- as.numeric(as.character(ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][,i]))
  #substitute with median
  ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][is.na(ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][,i]), i] <- median(ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][,i], na.rm = TRUE)
  #convert to factor
  ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][,i] <- as.factor(ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1')][,i])
  
}

###check variance
describe(ITSMdatafull[,c('SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1')])

for(i in 1:ncol(ITSMdatafull[,c('NoOfServ', 'Log_Empl_Client', 'Log_Empl_SP')])){
  ITSMdatafull[,c('NoOfServ', 'Log_Empl_Client', 'Log_Empl_SP')][is.na(ITSMdatafull[,c('NoOfServ', 'Log_Empl_Client', 'Log_Empl_SP')][,i]), i] <- mean(ITSMdatafull[,c('NoOfServ', 'Log_Empl_Client', 'Log_Empl_SP')][,i], na.rm = TRUE)
}

ITSMdatafull$IndustryFactor <- ITSMdatafull$Industry
###set empty to select
ITSMdatafull$IndustryFactor[is.na(ITSMdatafull$IndustryFactor)] <- 1
###set other to select
ITSMdatafull$IndustryFactor[ITSMdatafull$IndustryFactor==19] <- 1

#create new Industry ISIC Factor
ITSMdatafull$IndustryISIC <- as.character(ITSMdatafull$IndustryFactor)
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="1"] <- "Other"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="2"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="3"] <- "Financial and Insurance Activities"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="4"] <- "Utilities, Construction, and Trade"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="5"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="6"] <- "Utilities, Construction, and Trade"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="7"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="8"] <- "Utilities, Construction, and Trade"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="9"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="10"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="11"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="12"] <- "Professional, Adminstrative, and Public Services"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="13"] <- "Professional, Adminstrative, and Public Services"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="14"] <- "Information and Communication"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="15"] <- "Manufacturing"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="16"] <- "Transportation, Accomodation and Food Services"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="17"] <- "Transportation, Accomodation and Food Services"
ITSMdatafull$IndustryISIC[ITSMdatafull$IndustryISIC=="18"] <- "Financial and Insurance Activities"
ITSMdatafull$IndustryISIC <- as.factor(ITSMdatafull$IndustryISIC)

###we want "Other" not to show
ITSMdatafull <- within(ITSMdatafull, IndustryISIC <- relevel(IndustryISIC, ref = "Other"))


###check correlation
ITSMdatafull$ITSMCapScoreN <- as.numeric(ITSMdatafull$ITSMCapScore)
corApps <- hetcor(ITSMdatafull[,c('AlignScore','ITSMCapScoreN','PerformanceScore', 'SPf1_ProjMgmt1', 'SPf2_SoftwDev1', 'SPf3_ApplCust1', 'SPf4_ApplMgmt1', 'SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','SPf10_Other1', 'NoOfServ' , 'Log_Empl_Client' , 'RegScore')])
View(corApps$correlations)
View(corApps$tests)

aov1 = aov(AlignScore ~ IndustryISIC, data=ITSMdatafull)
summary(aov1)
aov2 = aov(ITSMCapScoreN ~ IndustryISIC, data=ITSMdatafull)
summary(aov2)
aov3 = aov(PerformanceScore ~ IndustryISIC, data=ITSMdatafull)
summary(aov3)

###more means and sds
mean(ITSMdatafull$Position_Vertical, na.rm=TRUE)
sd(ITSMdatafull$Position_Vertical, na.rm=TRUE)

mean(ITSMdatafull$PositionTime, na.rm=TRUE)
sd(ITSMdatafull$PositionTime, na.rm=TRUE)

mean(ITSMdatafull$Log_Empl_Client, na.rm=TRUE)
sd(ITSMdatafull$Log_Empl_Client, na.rm=TRUE)


#####include country information before regression
UserInputs_countryselected <- read.csv("C:/Users/jwulf/OneDrive/Research/IWI Forschung/JMIS/UserInputs_continent_time_countryselected.csv", sep=";")

##merge with the other data
ITSMdatafull_merged <- merge(ITSMdatafull,UserInputs_countryselected, by.x = "ID", by.y = "ID", all=FALSE, all.x=TRUE)
ITSMdatafull <- ITSMdatafull_merged

#we need to join Europe AND Africa
ITSMdatafull$Continent <- as.character(ITSMdatafull$Continent)
ITSMdatafull$Continent <- ifelse((ITSMdatafull$Continent == "Europe") | (ITSMdatafull$Continent == "Africa"), "EMEA",ITSMdatafull$Continent)
ITSMdatafull$Continent <- as.factor(ITSMdatafull$Continent)
ITSMdatafull$Continent <- relevel(ITSMdatafull$Continent, ref = "EMEA")



########Stepwise models ####################################################

######################DV: Alignment
###A1 and A2 model for alignment
align1 <- lm(AlignScore ~ ITSMCapScore, data = ITSMdatafull)
summary(align1) 
lm.beta(align1)
# 
align2 <- lm(AlignScore ~ ITSMCapScore*ITStratOrient  , data = ITSMdatafull)
summary(align2)
lm.beta(align2)
vif(align2)
modelCompare(align1, align2)

###stepwise models for performance
#only significant controls with size
perform1 <- lm(PerformanceScore ~  IndustryISIC+ Log_Empl_Client + RegScore + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)

summary(perform1)
lm.beta(perform1)

#only with significant controls with size
perform2 <- lm(PerformanceScore ~  ITSMCapScore+ Log_Empl_Client + IndustryISIC + RegScore + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)

summary(perform2)
lm.beta(perform2)
modelCompare(perform1, perform2)
# 

##E3a third predictor and mediator
#only significant controls with size
perform3a <- lm(PerformanceScore ~  ITSMCapScore + AlignScore+ Log_Empl_Client + IndustryISIC + RegScore + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)

summary(perform3a)
lm.beta(perform3a)
vif(perform3a)
modelCompare(perform2, perform3a)

perform4 <- lm(PerformanceScore ~  ITSMCapScore*ITStratOrient + AlignScore+ Log_Empl_Client  + IndustryISIC + RegScore + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)

summary(perform4)
lm.beta(perform4)
vif(perform4)
modelCompare(perform3a, perform4)



##A5. full

#only significant controls with size
perform5 <- lm(PerformanceScore ~  ITSMCapScore*ITStratOrient + AlignScore*ITStratOrient+ Log_Empl_Client + IndustryISIC + RegScore + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)

summary(perform5)
lm.beta(perform5)
vif(perform5)
modelCompare(perform4, perform5)


##########################################################
#######the same models with full controls as robustness check

#E1 first only controls
perform1f <- lm(PerformanceScore ~  IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
summary(perform1f)
lm.beta(perform1f)

#E2 second controls with predictor (cap)
perform2f <- lm(PerformanceScore ~  ITSMCapScore + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
summary(perform2f)
lm.beta(perform2f)
modelCompare(perform1f, perform2f)
# 

##E3a third predictor and mediator
perform3af <- lm(PerformanceScore ~  ITSMCapScore + AlignScore + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
summary(perform3af)
lm.beta(perform3af)
vif(perform3af)
modelCompare(perform2f, perform3af)

##A4. direct moderation model

perform4f <- lm(PerformanceScore ~  ITSMCapScore*ITStratOrient + AlignScore + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
summary(perform4f)
lm.beta(perform4f)
vif(perform4f)
modelCompare(perform3af, perform4f)

##A5. full
perform5f <- lm(PerformanceScore ~  ITSMCapScore*ITStratOrient + AlignScore*ITStratOrient + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes+ Continent, data = ITSMdatafull)
summary(perform5f)
lm.beta(perform5f)
vif(perform5f)
modelCompare(perform4f, perform5f)



###linear regression for first stage
firstStage <- lm(ITSMCapScore ~ SPf1_ProjMgmt1+SPf2_SoftwDev1+SPf3_ApplCust1+SPf4_ApplMgmt1+SPf5_DataCentr1+SPf6_NetwCom1+SPf7_DesktPrint1+SPf8_SupHlpDsk1+SPf9_Training1+Log_Empl_SP, data = ITSMdatafull)

summary(firstStage)
ITSMdatafull$firstResiduals <- residuals(firstStage)


#####Endogeneity for Alignment

secondStage <- lm(AlignScore ~ ITSMCapScore + firstResiduals, data = ITSMdatafull, weights = (1/abs(firstResiduals)^2))
summary(secondStage)
ncvTest(secondStage)

##check endogeneity instrument diagnostics
##interpretation guidelines here: https://stats.stackexchange.com/questions/134789/interpretation-of-ivreg-diagnostics-in-r
alignivreg <- ivreg(AlignScore ~ ITSMCapScore | SPf1_ProjMgmt1+SPf2_SoftwDev1+SPf3_ApplCust1+SPf4_ApplMgmt1+SPf5_DataCentr1+SPf6_NetwCom1+SPf7_DesktPrint1+SPf8_SupHlpDsk1+SPf9_Training1+Log_Empl_SP, data = ITSMdatafull)
summary(alignivreg, diagnostics= TRUE)

#####Endogeneity for Performance
###
secondStage <- lm(PerformanceScore ~ ITSMCapScore + firstResiduals, data = ITSMdatafull, weights = (1/abs(firstResiduals)))
summary(secondStage)
ncvTest(secondStage)


###now with controls
secondStagePerform2 <- lm(PerformanceScore ~  ITSMCapScore + firstResiduals + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime, data = ITSMdatafull, weights = (1/abs(firstResiduals)))
summary(secondStagePerform2)
ncvTest(secondStagePerform2)

secondStagePerform3 <- lm(PerformanceScore ~ ITSMCapScore*ITStratOrient + firstResiduals, data = ITSMdatafull)
summary(secondStagePerform3)
ncvTest(secondStagePerform3)

###check endogeneity statistics
performivreg <- ivreg(PerformanceScore ~ ITSMCapScore*ITStratOrient | SPf1_ProjMgmt1+SPf2_SoftwDev1+SPf3_ApplCust1+SPf4_ApplMgmt1+SPf5_DataCentr1+SPf6_NetwCom1+SPf7_DesktPrint1+SPf8_SupHlpDsk1+SPf9_Training1+Log_Empl_SP, data = ITSMdatafull)
summary(performivreg, diagnostics= TRUE)

##test ivmodel
library(ivmodel)

#Y: outcome
Y<-ITSMdatafull$PerformanceScore

#D: endogenous variable (JUST ONE!!)
D<-ITSMdatafull$ITSMCapScore

#Z: data frame instruments
Zname<-c('SPf1_ProjMgmt1','SPf2_SoftwDev1','SPf3_ApplCust1','SPf4_ApplMgmt1','SPf5_DataCentr1','SPf6_NetwCom1','SPf7_DesktPrint1','SPf8_SupHlpDsk1','SPf9_Training1','Log_Empl_SP')
Z<-ITSMdatafull[,Zname]

#X: exogeneous controls (include the other predictors here...)
ITSMdatafull$ITSMCapScoreXITStratOrient<-ITSMdatafull$ITSMCapScore*ITSMdatafull$ITStratOrient
X<-ITSMdatafull[,c('ITSMCapScoreXITStratOrient','ITStratOrient')]

ivmodel <- ivmodel(Y=Y, D=D, Z=Z, X=X)
AR.test(ivmodel, beta0=0, alpha=0.05)



################################TESTS for moderated mediation#################

###Control Prediction: Client Size, INdustry, Position HOr, Position Vertical, Position Tenure, Marker
for(i in 1:ncol(ITSMdatafull[,c('Position_Vertical', 'PositionBusSP01', 'PositionTime')])){
  ITSMdatafull[,c('Position_Vertical', 'PositionBusSP01', 'PositionTime')][is.na(ITSMdatafull[,c('Position_Vertical', 'PositionBusSP01', 'PositionTime')][,i]), i] <- median(ITSMdatafull[,c('Position_Vertical', 'PositionBusSP01', 'PositionTime')][,i], na.rm = TRUE)
}

###################direct effect moderation model
###the bootstrap function for testing direct effect moderation model with mean centered IT strat or
bootfunction = function (dataset, random) {
  d = dataset[random, ] ###randomize by row
  d$ITStratOrientp1 = I(d$ITStratOrient-mean(d$ITStratOrient))+sd(d$ITStratOrient)
  d$ITStratOrientm1 = I(d$ITStratOrient-mean(d$ITStratOrient))-sd(d$ITStratOrient)
  apath = lm(AlignScore ~ I(ITSMCapScore-mean(ITSMCapScore)) , data = d)
  bpathp1 =  lm(PerformanceScore ~  I(ITSMCapScore-mean(ITSMCapScore))*ITStratOrientp1 + AlignScore + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime, data = d)
  bpathm1 =  lm(PerformanceScore ~  I(ITSMCapScore-mean(ITSMCapScore))*ITStratOrientm1 + AlignScore + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime, data = d)
  
  indirectp1 = apath$coefficients[2]*bpathp1$coefficients[4]
  indirectm1 = apath$coefficients[2]*bpathm1$coefficients[4]
  indirectdiff = indirectp1-indirectm1
  
  directp1=bpathp1$coefficients[2]+bpathp1$coefficients[17]
  directm1=bpathm1$coefficients[2]+bpathm1$coefficients[17]
  directdiff = directp1-directm1
  
  totalp1=directp1+indirectp1
  totalm1=directm1+indirectm1
  totaldiff = totalp1-totalm1
  
  first=apath$coefficients[2]
  secondp1=bpathp1$coefficients[4]
  secondm1=bpathm1$coefficients[4]
  
  return(c(indirectp1,indirectm1,indirectdiff,directp1,directm1,directdiff,totalp1,totalm1,totaldiff, first, secondp1, secondm1))
}

###do the bootstrapping 
bootresults = boot(data = ITSMdatafull, statistic = bootfunction, R = 1000)
boot.ci(bootresults, conf=.99, type='bca', index=1)#indirectp1
boot.ci(bootresults, conf=.99, type='bca', index=2)#indirectm1
boot.ci(bootresults, conf=.99, type='bca', index=3)#indirectdiff
boot.ci(bootresults, conf=.99, type='bca', index=4)#directp1
boot.ci(bootresults, conf=.99, type='bca', index=5)#directm1
boot.ci(bootresults, conf=.99, type='bca', index=6)#directdiff
boot.ci(bootresults, conf=.99, type='bca', index=7)#totalp1
boot.ci(bootresults, conf=.99, type='bca', index=8)#totalm1
boot.ci(bootresults, conf=.99, type='bca', index=9)#totaldiff
boot.ci(bootresults, conf=.99, type='bca', index=10)#first
boot.ci(bootresults, conf=.99, type='bca', index=11)#secondp1
boot.ci(bootresults, conf=.99, type='bca', index=12)#secondm1

boot.ci(bootresults, conf=.95, type='bca', index=1)#indirectp1
boot.ci(bootresults, conf=.95, type='bca', index=2)#indirectm1
boot.ci(bootresults, conf=.95, type='bca', index=3)#indirectdiff
boot.ci(bootresults, conf=.95, type='bca', index=4)#directp1
boot.ci(bootresults, conf=.95, type='bca', index=5)#directm1
boot.ci(bootresults, conf=.95, type='bca', index=6)#directdiff
boot.ci(bootresults, conf=.95, type='bca', index=7)#totalp1
boot.ci(bootresults, conf=.95, type='bca', index=8)#totalm1
boot.ci(bootresults, conf=.95, type='bca', index=9)#totaldiff
boot.ci(bootresults, conf=.95, type='bca', index=10)#first
boot.ci(bootresults, conf=.95, type='bca', index=11)#secondp1
boot.ci(bootresults, conf=.95, type='bca', index=12)#secondm1

#####bootstrapping for the full effects model
modtestfull = function (dataset, random) {
  d = dataset[random, ] ###randomize by row
  apath = lm(AlignScore ~ ITSMCapScore*ITStrat_Cons + ITSMCapScore*ITStrat_Inno, data = d)
  return(coef(apath))
}








####################now we bootstrap model E4
perform8a <- lm(PerformanceScore ~  ITSMCapScore*ITStratOrient + AlignScore + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
summary(perform8a)


####in preparation for boostrapping, we want mean centered ITSMCapScore, ITStratOrient und AlignScore
ITSMdatafull$ITSMCapScoreMC <- ITSMdatafull$ITSMCapScore-mean(ITSMdatafull$ITSMCapScore)
ITSMdatafull$AlignScoreMC <- ITSMdatafull$AlignScore-mean(ITSMdatafull$AlignScore)
ITSMdatafull$ITStratOrientMC <- ITSMdatafull$ITStratOrient-mean(ITSMdatafull$ITStratOrient)


###the bootstrap function for testing mediation and direct moderated path (see Edwards and Lambert 2007)
ITStratOrientp1 = mean(ITSMdatafull$ITStratOrientMC)+sd(ITSMdatafull$ITStratOrientMC)
ITStratOrientm1 = mean(ITSMdatafull$ITStratOrientMC)-sd(ITSMdatafull$ITStratOrientMC)
boottotaleffect = function (dataset, random) {
  d = dataset[random, ] ###randomize by row
  apath = lm(AlignScoreMC ~ ITSMCapScoreMC*ITStratOrientMC , data = d)
  #bpath = lm(PerformanceScore ~  ITSMCapScoreMC*ITStratOrientMC + AlignScoreMC+ Log_Empl_Client + IndustryISIC + RegScore + PositionBusSP01 + PositionTime_6classes  + Continent, data = d)
  bpath = lm(PerformanceScore ~  ITSMCapScoreMC*ITStratOrientMC + AlignScoreMC + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = d)
  
  #first
  firstp1=coef(apath)[2]+coef(apath)[4]*ITStratOrientp1
  firstm1=coef(apath)[2]+coef(apath)[4]*ITStratOrientm1
  firstdiff=firstp1-firstm1
  
  #second
  secondp1=coef(bpath)[4]
  secondm1=coef(bpath)[4]
  seconddiff=secondp1-secondm1
  
  #direct:
  directp1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientp1
  directm1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientm1
  directdiff = directp1 - directm1
  
  #indirect:
  indirectp1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientp1)*(coef(bpath)[4])
  indirectm1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientm1)*(coef(bpath)[4])
  indirectdiff = indirectp1 - indirectm1  
  
  #total
  totalp1=directp1+indirectp1
  totalm1=directm1+indirectm1
  totaldiff = totalp1 - totalm1
  
  return(c(firstp1, firstm1, firstdiff, secondp1, secondm1, seconddiff, directp1, directm1, directdiff, indirectp1, indirectm1, indirectdiff, totalp1, totalm1, totaldiff))
}


###bootstrap until solution significant
notstop=TRUE
i<-0
while(notstop){
  i<-i+1
  bootresults5 = boot(data = ITSMdatafull, statistic = boottotaleffect, R = 1000)
  direc.ci <- boot.ci(bootresults5, conf=.95, type='norm', index=9)#directdiff
  notstop = direc.ci$normal[2]<0
  print(i)
  flush.console()
}

#.95
boot.ci(bootresults5, conf=.95, type='norm', index=1)#firstp1
boot.ci(bootresults5, conf=.95, type='norm', index=2)#firstm1
boot.ci(bootresults5, conf=.95, type='norm', index=3)#firstdiff
boot.ci(bootresults5, conf=.95, type='norm', index=4)#secondp1
boot.ci(bootresults5, conf=.95, type='norm', index=5)#secondm1
boot.ci(bootresults5, conf=.95, type='norm', index=6)#seconddiff
boot.ci(bootresults5, conf=.95, type='norm', index=7)#directp1
boot.ci(bootresults5, conf=.95, type='norm', index=8)#directm1
boot.ci(bootresults5, conf=.95, type='norm', index=9)#directdiff
boot.ci(bootresults5, conf=.95, type='norm', index=10)#indirectp1
boot.ci(bootresults5, conf=.95, type='norm', index=11)#indirectm1
boot.ci(bootresults5, conf=.95, type='norm', index=12)#indirectdiff
boot.ci(bootresults5, conf=.95, type='norm', index=13)#totalp1
boot.ci(bootresults5, conf=.95, type='norm', index=14)#totalm1
boot.ci(bootresults5, conf=.95, type='norm', index=15)#totaldiff

#.99
boot.ci(bootresults5, conf=.99, type='norm', index=1)#firstp1
boot.ci(bootresults5, conf=.99, type='norm', index=2)#firstm1
boot.ci(bootresults5, conf=.99, type='norm', index=3)#firstdiff
boot.ci(bootresults5, conf=.99, type='norm', index=4)#secondp1
boot.ci(bootresults5, conf=.99, type='norm', index=5)#secondm1
boot.ci(bootresults5, conf=.99, type='norm', index=6)#seconddiff
boot.ci(bootresults5, conf=.99, type='norm', index=7)#directp1
boot.ci(bootresults5, conf=.99, type='norm', index=8)#directm1
boot.ci(bootresults5, conf=.99, type='norm', index=9)#directdiff
boot.ci(bootresults5, conf=.99, type='norm', index=10)#indirectp1
boot.ci(bootresults5, conf=.99, type='norm', index=11)#indirectm1
boot.ci(bootresults5, conf=.99, type='norm', index=12)#indirectdiff
boot.ci(bootresults5, conf=.99, type='norm', index=13)#totalp1
boot.ci(bootresults5, conf=.99, type='norm', index=14)#totalm1
boot.ci(bootresults5, conf=.99, type='norm', index=15)#totaldiff


#now lets calculate the simple effects
apath = lm(AlignScoreMC ~ ITSMCapScoreMC*ITStratOrientMC , data = ITSMdatafull)
apath
#bpath = lm(PerformanceScore ~  ITSMCapScoreMC*ITStratOrientMC + AlignScoreMC+ Log_Empl_Client + IndustryISIC +  RegScore + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
bpath = lm(PerformanceScore ~  ITSMCapScoreMC*ITStratOrientMC + AlignScoreMC + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)

bpath
#first
firstp1=coef(apath)[2]+coef(apath)[4]*ITStratOrientp1
firstp1
firstm1=coef(apath)[2]+coef(apath)[4]*ITStratOrientm1
firstm1
firstdiff=firstp1-firstm1
firstdiff

#second
secondp1=coef(bpath)[4]
secondp1
secondm1=coef(bpath)[4]
secondm1
seconddiff=secondp1-secondm1
seconddiff

#direct:
directp1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientp1
directp1
directm1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientm1
directm1
directdiff = directp1 - directm1
directdiff

#indirect:
indirectp1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientp1)*(coef(bpath)[4])
indirectp1
indirectm1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientm1)*(coef(bpath)[4])
indirectm1
indirectdiff = indirectp1 - indirectm1  
indirectdiff

#total
totalp1=directp1+indirectp1
totalp1
totalm1=directm1+indirectm1
totalm1
totaldiff = totalp1 - totalm1
totaldiff 


####this is bootstrapping the full model

###################full moderated mediation test

perform9a <- lm(PerformanceScore ~  ITSMCapScore*ITStratOrient + AlignScore*ITStratOrient + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
summary(perform9a)



####in preparation for boostrapping, we want mean centered ITSMCapScore, ITStratOrient und AlignScore
ITSMdatafull$ITSMCapScoreMC <- ITSMdatafull$ITSMCapScore-mean(ITSMdatafull$ITSMCapScore)
ITSMdatafull$AlignScoreMC <- ITSMdatafull$AlignScore-mean(ITSMdatafull$AlignScore)
ITSMdatafull$ITStratOrientMC <- ITSMdatafull$ITStratOrient-mean(ITSMdatafull$ITStratOrient)


###the bootstrap function for testing mediation and direct moderated path (see Edwards and Lambert 2007)
ITStratOrientp1 = mean(ITSMdatafull$ITStratOrientMC)+sd(ITSMdatafull$ITStratOrientMC)
ITStratOrientm1 = mean(ITSMdatafull$ITStratOrientMC)-sd(ITSMdatafull$ITStratOrientMC)
boottotaleffect = function (dataset, random) {
  d = dataset[random, ] ###randomize by row
  apath = lm(AlignScoreMC ~ ITSMCapScoreMC*ITStratOrientMC , data = d)
  bpath = lm(PerformanceScore ~  ITSMCapScoreMC*ITStratOrientMC + AlignScoreMC*ITStratOrientMC + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = d)

  #first
  firstp1=coef(apath)[2]+coef(apath)[4]*ITStratOrientp1
  firstm1=coef(apath)[2]+coef(apath)[4]*ITStratOrientm1
  firstdiff=firstp1-firstm1
  
  #second
  secondp1=coef(bpath)[4]+coef(bpath)[20]*ITStratOrientp1
  secondm1=coef(bpath)[4]+coef(bpath)[20]*ITStratOrientm1
  seconddiff=secondp1-secondm1
  
  #direct:
  directp1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientp1
  directm1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientm1
  directdiff = directp1 - directm1
  
  #indirect:
  indirectp1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientp1)*(coef(bpath)[4]+coef(bpath)[20]*ITStratOrientp1)
  indirectm1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientm1)*(coef(bpath)[4]+coef(bpath)[20]*ITStratOrientm1)
  indirectdiff = indirectp1 - indirectm1  
  
  #total
  totalp1=directp1+indirectp1
  totalm1=directm1+indirectm1
  totaldiff = totalp1 - totalm1
  
  return(c(firstp1, firstm1, firstdiff, secondp1, secondm1, seconddiff, directp1, directm1, directdiff, indirectp1, indirectm1, indirectdiff, totalp1, totalm1, totaldiff))
}


###bootstrap until solution significant
notstop=TRUE
i<-0
while(notstop){
  i<-i+1
  bootresults5 = boot(data = ITSMdatafull, statistic = boottotaleffect, R = 1000)
  direc.ci <- boot.ci(bootresults5, conf=.9, type='norm', index=9)#directdiff
  notstop = direc.ci$normal[2]<0
  
  #direc.ci <- boot.ci(bootresults5, conf=.9, type='norm', index=9)#directdiff
  #direc.ci2 <- boot.ci(bootresults5, conf=.9, type='norm', index=3)
  #stopping = direc.ci$normal[2]>0 & direc.ci2$normal[2]<0
  #notstop = !stopping
  print(i)
  flush.console()
}

#.9
boot.ci(bootresults5, conf=.9, type='norm', index=1)#firstp1
boot.ci(bootresults5, conf=.9, type='norm', index=2)#firstm1
boot.ci(bootresults5, conf=.9, type='norm', index=3)#firstdiff
boot.ci(bootresults5, conf=.9, type='norm', index=4)#secondp1
boot.ci(bootresults5, conf=.9, type='norm', index=5)#secondm1
boot.ci(bootresults5, conf=.9, type='norm', index=6)#seconddiff
boot.ci(bootresults5, conf=.9, type='norm', index=7)#directp1
boot.ci(bootresults5, conf=.9, type='norm', index=8)#directm1
boot.ci(bootresults5, conf=.9, type='norm', index=9)#directdiff
boot.ci(bootresults5, conf=.9, type='norm', index=10)#indirectp1
boot.ci(bootresults5, conf=.9, type='norm', index=11)#indirectm1
boot.ci(bootresults5, conf=.9, type='norm', index=12)#indirectdiff
boot.ci(bootresults5, conf=.9, type='norm', index=13)#totalp1
boot.ci(bootresults5, conf=.9, type='norm', index=14)#totalm1
boot.ci(bootresults5, conf=.9, type='norm', index=15)#totaldiff

#.9
boot.ci(bootresults5, conf=.925, type='norm', index=1)#firstp1
boot.ci(bootresults5, conf=.925, type='norm', index=2)#firstm1
boot.ci(bootresults5, conf=.925, type='norm', index=3)#firstdiff
boot.ci(bootresults5, conf=.925, type='norm', index=4)#secondp1
boot.ci(bootresults5, conf=.925, type='norm', index=5)#secondm1
boot.ci(bootresults5, conf=.925, type='norm', index=6)#seconddiff
boot.ci(bootresults5, conf=.925, type='norm', index=7)#directp1
boot.ci(bootresults5, conf=.925, type='norm', index=8)#directm1
boot.ci(bootresults5, conf=.925, type='norm', index=9)#directdiff
boot.ci(bootresults5, conf=.925, type='norm', index=10)#indirectp1
boot.ci(bootresults5, conf=.925, type='norm', index=11)#indirectm1
boot.ci(bootresults5, conf=.925, type='norm', index=12)#indirectdiff
boot.ci(bootresults5, conf=.925, type='norm', index=13)#totalp1
boot.ci(bootresults5, conf=.925, type='norm', index=14)#totalm1
boot.ci(bootresults5, conf=.925, type='norm', index=15)#totaldiff

#.95
boot.ci(bootresults5, conf=.95, type='norm', index=1)#firstp1
boot.ci(bootresults5, conf=.95, type='norm', index=2)#firstm1
boot.ci(bootresults5, conf=.95, type='norm', index=3)#firstdiff
boot.ci(bootresults5, conf=.95, type='norm', index=4)#secondp1
boot.ci(bootresults5, conf=.95, type='norm', index=5)#secondm1
boot.ci(bootresults5, conf=.95, type='norm', index=6)#seconddiff
boot.ci(bootresults5, conf=.95, type='norm', index=7)#directp1
boot.ci(bootresults5, conf=.95, type='norm', index=8)#directm1
boot.ci(bootresults5, conf=.95, type='norm', index=9)#directdiff
boot.ci(bootresults5, conf=.95, type='norm', index=10)#indirectp1
boot.ci(bootresults5, conf=.95, type='norm', index=11)#indirectm1
boot.ci(bootresults5, conf=.95, type='norm', index=12)#indirectdiff
boot.ci(bootresults5, conf=.95, type='norm', index=13)#totalp1
boot.ci(bootresults5, conf=.95, type='norm', index=14)#totalm1
boot.ci(bootresults5, conf=.95, type='norm', index=15)#totaldiff

#.99
boot.ci(bootresults5, conf=.99, type='norm', index=1)#firstp1
boot.ci(bootresults5, conf=.99, type='norm', index=2)#firstm1
boot.ci(bootresults5, conf=.99, type='norm', index=3)#firstdiff
boot.ci(bootresults5, conf=.99, type='norm', index=4)#secondp1
boot.ci(bootresults5, conf=.99, type='norm', index=5)#secondm1
boot.ci(bootresults5, conf=.99, type='norm', index=6)#seconddiff
boot.ci(bootresults5, conf=.99, type='norm', index=7)#directp1
boot.ci(bootresults5, conf=.99, type='norm', index=8)#directm1
boot.ci(bootresults5, conf=.99, type='norm', index=9)#directdiff
boot.ci(bootresults5, conf=.99, type='norm', index=10)#indirectp1
boot.ci(bootresults5, conf=.99, type='norm', index=11)#indirectm1
boot.ci(bootresults5, conf=.99, type='norm', index=12)#indirectdiff
boot.ci(bootresults5, conf=.99, type='norm', index=13)#totalp1
boot.ci(bootresults5, conf=.99, type='norm', index=14)#totalm1
boot.ci(bootresults5, conf=.99, type='norm', index=15)#totaldiff


#now lets calculate the simple effects
apath = lm(AlignScoreMC ~ ITSMCapScoreMC*ITStratOrientMC , data = ITSMdatafull)
apath
bpath = lm(PerformanceScore ~  ITSMCapScoreMC*ITStratOrientMC + AlignScoreMC*ITStratOrientMC + IndustryISIC + Log_Empl_Client + IndScore + RegScore + Position_Vertical + PositionBusSP01 + PositionTime_6classes + Continent, data = ITSMdatafull)
bpath
#first
firstp1=coef(apath)[2]+coef(apath)[4]*ITStratOrientp1
firstp1
firstm1=coef(apath)[2]+coef(apath)[4]*ITStratOrientm1
firstm1
firstdiff=firstp1-firstm1
firstdiff

#second
secondp1=coef(bpath)[4]+coef(bpath)[20]*ITStratOrientp1
secondp1
secondm1=coef(bpath)[4]+coef(bpath)[20]*ITStratOrientm1
secondm1
seconddiff=secondp1-secondm1
seconddiff

#direct:
directp1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientp1
directp1
directm1=coef(bpath)[2]+coef(bpath)[19]*ITStratOrientm1
directm1
directdiff = directp1 - directm1
directdiff

#indirect:
indirectp1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientp1)*(coef(bpath)[4]+coef(bpath)[20]*ITStratOrientp1)
indirectp1
indirectm1=(coef(apath)[2]+coef(apath)[4]*ITStratOrientm1)*(coef(bpath)[4]+coef(bpath)[20]*ITStratOrientm1)
indirectm1
indirectdiff = indirectp1 - indirectm1  
indirectdiff

#total
totalp1=directp1+indirectp1
totalp1
totalm1=directm1+indirectm1
totalm1
totaldiff = totalp1 - totalm1
totaldiff 



