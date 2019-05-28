
#calculate the loadings and report on model quality
#for two models: all data (with itil_mean) and C12 data (with ITILcap)

library(lavaan)
library(semTools)
library(Hmisc)
library(QuantPsyc)

#lets get the data first

##Import Data
ExpertData <- read.csv("C:/Users/jwulf/OneDrive/Research/IWI Forschung/GFF/ITSM Data/20161004_Experts/ExpertData.csv", sep=";", na.strings="")

##IMport extensions
ExpertDataExtensions <- read.csv("C:/Users/jwulf/OneDrive/Research/IWI Forschung/GFF/ITSM Data/20161004_Experts/UserInputs (24)_selected_n541_selected_instruments.csv", sep=";", na.strings="")


#add the data we need
ExpertData$Class<-as.factor(ifelse(ExpertData$P1>.85,1,ifelse(ExpertData$P2>.85,2,0)))


##add Class Dummies
ExpertData$C1<-ifelse(ExpertData$P1>.85,1,0)
ExpertData$C2<-ifelse(ExpertData$P2>.85,1,0)
ExpertData$C0<-ifelse(ExpertData$MAXPI<=.85,1,0)

##separate Frames for classes
ExpertDataC1<-ExpertData[ExpertData$C1==1,]
ExpertDataC2<-ExpertData[ExpertData$C2==1,]

#get the csvs
class2PersLoc<- read.csv("C:/Users/jwulf/OneDrive/Research/IWI Forschung/GFF/ITSM Data/20161004_Experts/Final Analyzes/PersLocAnalyses/2classes_PersLoc.csv", sep=";", na.strings="")
class2C2PersLoc<- read.csv("C:/Users/jwulf/OneDrive/Research/IWI Forschung/GFF/ITSM Data/20161004_Experts/Final Analyzes/PersLocAnalyses/2classes2_PersLoc.csv", sep=";", na.strings="")

#merge
C1merge <- merge(x=ExpertDataC1,y=class2PersLoc,by.x="ID",by.y="c1it")
C2merge <- merge(x=ExpertDataC2,y=class2C2PersLoc,by.x="ID",by.y="c2it")

#rename class2PersLoc <= class2C2PersLoc
names(C1merge)[names(C1merge)=="c1perspar"] <- "perspar"
names(C2merge)[names(C2merge)=="c2perspar"] <- "perspar"


#create standardized person parameters
C1merge$persparSTA <- scale(C1merge$perspar)
C2merge$persparSTA <- scale(C2merge$perspar)

#bind c1merge and c2merge into a single dataframe where
ExpertsC12merge <- rbind(C1merge,C2merge)

#median impute
ExpertsC12merge$EmpC[is.na(ExpertsC12merge$EmpC)] <- median(ExpertsC12merge$EmpC, na.rm=TRUE)
ExpertsC12merge$EmplSP[is.na(ExpertsC12merge$EmplSP)] <- median(ExpertsC12merge$EmplSP, na.rm=TRUE)

#now merge with extracolumns
Experts_new <- merge(ExpertsC12merge, ExpertDataExtensions, by.x = "ID", by.y = "ID")
ExpertsC12merge <- Experts_new 


#model for efa, regression and chow
itsmmod <- 'Performance =~ Perform1+Perform2+Perform3+Perform4+Perform5
FirmSize =~ EmpC+ EmplSP
Regulation =~ Reg1 + Reg2
IndServ =~ Ind1+ Ind2'

#fit the models
cfa_fit<- cfa(itsmmod, data = ExpertsC12merge, estimator = "DWLS")

#give out reliability indicators, korrelationcoefficients and p values
reliability(cfa_fit)

#Build factor scores and correlation coeffs
faccoef<-predict(cfa_fit)


indcorr <- rcorr(faccoef)
indcorr$r
indcorr$P

#now htmt
htmt( ExpertsC12merge[,c('Perform1','Perform2','Perform3','Perform4','Perform5','EmpC','EmplSP','Reg1','Reg2','Ind1','Ind2')], itsmmod)

#exploratory factor analysis 
#this is for Harman's test 
faALL1 <- factanal(na.omit(ExpertsC12merge[,c('Perform1', 'Perform2', 'Perform3', 'Perform4', 'Perform5', 'EmpC', 'EmplSP', 'Ind1', 'Ind2', 'Reg1', 'Reg2')]), factors=1, rotation="varimax", scores="regression")
print(faALL1$loadings, cutoff = 0)

#this is the EFA
faALL <- factanal(na.omit(ExpertsC12merge[,c('Perform1', 'Perform2', 'Perform3', 'Perform4', 'Perform5', 'EmpC', 'EmplSP', 'Ind1', 'Ind2', 'Reg1', 'Reg2')]), factors=4, rotation="varimax", scores="regression")
print(faALL$loadings, cutoff = 0)


###Prepare additional variables
#now bind the cortable to expertsc12merge and split
ExpertsC12scores<-cbind(ExpertsC12merge,faccoef)

#create class2 dummy
ExpertsC12scores$class2d<-as.factor(ifelse(ExpertsC12scores$Class=="1",1,0))

#industry uncertainty
ExpertsC12scores$IndUnc <- ifelse(ExpertsC12scores$Industry==1, 22.58333333,
                                  ifelse(ExpertsC12scores$Industry==2, 17.8,
                                         ifelse(ExpertsC12scores$Industry==3, 17.1,
                                                ifelse(ExpertsC12scores$Industry==4, 22.58333333,
                                                       ifelse(ExpertsC12scores$Industry==5, 19.5,
                                                              ifelse(ExpertsC12scores$Industry==6, 20.9,
                                                                     ifelse(ExpertsC12scores$Industry==7, 17.8,
                                                                            ifelse(ExpertsC12scores$Industry==8, 30.6,
                                                                                   ifelse(ExpertsC12scores$Industry==9, 17.8,
                                                                                          ifelse(ExpertsC12scores$Industry==10, 19.5,
                                                                                                 ifelse(ExpertsC12scores$Industry==11, 29.6,
                                                                                                        ifelse(ExpertsC12scores$Industry==12, 22.58333333,
                                                                                                               ifelse(ExpertsC12scores$Industry==13, 22.58333333,
                                                                                                                      ifelse(ExpertsC12scores$Industry==14, 22.58333333,
                                                                                                                             ifelse(ExpertsC12scores$Industry==15, 29.6,
                                                                                                                                    ifelse(ExpertsC12scores$Industry==16, 22.58333333,
                                                                                                                                           ifelse(ExpertsC12scores$Industry==17, 30.6,
                                                                                                                                                  ifelse(ExpertsC12scores$Industry==18, 17.1, 22.58333333))))))))))))))))))


ExpertsC12scores$Industry<-as.factor(ExpertsC12scores$Industry)


#external vs internal
ExpertsC12scores$External <- ifelse(ExpertsC12scores$Level=="5Extern",1,0)
ExpertsC12scores$External <- as.factor(ExpertsC12scores$External)

#managementlevel as scales
ExpertsC12scores$MgmtLevel <- as.integer(ExpertsC12scores$Level)
ExpertsC12scores$MgmtLevel <- ifelse(ExpertsC12scores$MgmtLevel>4,NA,ExpertsC12scores$MgmtLevel)
ExpertsC12scores$MgmtLevel[is.na(ExpertsC12scores$MgmtLevel)] <- median(ExpertsC12scores$MgmtLevel, na.rm=TRUE)


ExpertsC12scores$Pos <- as.character(ExpertsC12scores$Pos)
ExpertsC12scores$Pos[ExpertsC12scores$Pos=="NULL"] <- "ITSP"
ExpertsC12scores$Pos <- as.factor(ExpertsC12scores$Pos)

ExpertsC12scores$SPType <- as.character(ExpertsC12scores$SPType)
ExpertsC12scores$SPType[ExpertsC12scores$SPType=="NULL"] <- "1"
ExpertsC12scores$SPType <- as.numeric(ExpertsC12scores$SPType)
ExpertsC12scores$SPType <- as.factor(ExpertsC12scores$SPType)

##set the industry
ExpertsC12scores$IndustryFactor <- ExpertsC12scores$Industry
###set empty to select
ExpertsC12scores$IndustryFactor[is.na(ExpertsC12scores$IndustryFactor)] <- 1
###set other to select
ExpertsC12scores$IndustryFactor[ExpertsC12scores$IndustryFactor==19] <- 1

#create new Industry ISIC Factor
ExpertsC12scores$IndustryISIC <- as.character(ExpertsC12scores$IndustryFactor)
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="1"] <- "Other"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="2"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="3"] <- "Financial and Insurance Activities"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="4"] <- "Utilities, Construction, and Trade"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="5"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="6"] <- "Utilities, Construction, and Trade"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="7"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="8"] <- "Utilities, Construction, and Trade"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="9"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="10"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="11"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="12"] <- "Professional, Adminstrative, and Public Services"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="13"] <- "Professional, Adminstrative, and Public Services"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="14"] <- "Information and Communication"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="15"] <- "Manufacturing"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="16"] <- "Transportation, Accomodation and Food Services"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="17"] <- "Transportation, Accomodation and Food Services"
ExpertsC12scores$IndustryISIC[ExpertsC12scores$IndustryISIC=="18"] <- "Financial and Insurance Activities"
ExpertsC12scores$IndustryISIC <- as.factor(ExpertsC12scores$IndustryISIC)
###we want "Other" not to show
ExpertsC12scores <- within(ExpertsC12scores, IndustryISIC <- relevel(IndustryISIC, ref = "Other"))



##correlations
ExpertsC12scores$Pos <- as.integer(ExpertsC12scores$Pos)
selma <- ExpertsC12scores[,c('Performance','persparSTA','FirmSize','Regulation','MgmtLevel','SPType','IndServ','Pos','IndUnc')]
rcorr(as.matrix(selma))
ExpertsC12scores$Pos <- as.factor(ExpertsC12scores$Pos)


#next split into both classes
ExpertsC1scores<-ExpertsC12scores[ExpertsC12scores$Class==1,]
ExpertsC2scores<-ExpertsC12scores[ExpertsC12scores$Class==2,]

ExpertsC2scores$Pos <- as.integer(ExpertsC2scores$Pos)
selma <- ExpertsC2scores[,c('Performance','persparSTA','FirmSize','Regulation','MgmtLevel','SPType','IndServ','Pos','IndUnc')]
rcorr(as.matrix(selma))
ExpertsC2scores$Pos <- as.factor(ExpertsC2scores$Pos)

ExpertsC1scores$Pos <- as.integer(ExpertsC1scores$Pos)
selma <- ExpertsC1scores[,c('Performance','persparSTA','FirmSize','Regulation','MgmtLevel','SPType','IndServ','Pos','IndUnc')]
rcorr(as.matrix(selma))
ExpertsC1scores$Pos <- as.factor(ExpertsC1scores$Pos)

#check if we can through at $SPType
C1SPType<-lm(ExpertsC1scores$Performance~ExpertsC1scores$SPType)
summary(C1SPType)






#this is with all controls:
C12lm6<-lm(ExpertsC12scores$Performance~ExpertsC12scores$persparSTA+ExpertsC12scores$FirmSize+ExpertsC12scores$Regulation+ExpertsC12scores$MgmtLevel+ExpertsC12scores$Pos+ExpertsC12scores$SPType+ExpertsC12scores$IndServ+ExpertsC12scores$IndUnc)
summary(C12lm6)
lm.beta(C12lm6)

##best version for SPType categorical
#this works with estimator = "DWLS"
C2lmmod6<-lm(ExpertsC2scores$Performance~ExpertsC2scores$persparSTA+ExpertsC2scores$FirmSize+ExpertsC2scores$SPType+ExpertsC2scores$IndUnc+ExpertsC2scores$IndServ+ExpertsC2scores$MgmtLevel)
summary(C2lmmod6)
lm.beta(C2lmmod6)

C1lmmod6<-lm(ExpertsC1scores$Performance~ExpertsC1scores$persparSTA+ExpertsC1scores$FirmSize+ExpertsC1scores$SPType+ExpertsC1scores$IndUnc+ExpertsC1scores$IndServ+ExpertsC1scores$MgmtLevel)
summary(C1lmmod6)
lm.beta(C1lmmod6)



###these are the full models

C1lmmod6<-lm(ExpertsC1scores$Performance~ExpertsC1scores$persparSTA+ExpertsC1scores$FirmSize+ExpertsC1scores$Regulation+ExpertsC1scores$MgmtLevel+ExpertsC1scores$Pos+ExpertsC1scores$SPType+ExpertsC1scores$IndServ+ExpertsC1scores$IndUnc)
summary(C1lmmod6)
lm.beta(C1lmmod6)

C2lmmod6<-lm(ExpertsC2scores$Performance~ExpertsC2scores$persparSTA+ExpertsC2scores$FirmSize+ExpertsC2scores$Regulation+ExpertsC2scores$MgmtLevel+ExpertsC2scores$Pos+ExpertsC2scores$SPType+ExpertsC2scores$IndServ+ExpertsC2scores$IndUnc)
summary(C2lmmod6)
lm.beta(C2lmmod6)


##############################################################################################################
####NOW Prepare the chow test as described in kennedy a guide to econometrics 2008 6E, pp. 238-239
#we assume only the intercept and the performance influence to change 
#first create the interaction term
ExpertsC12scores$class2dint <- as.integer(ExpertsC12scores$class2d)-1 #factor to dummy
ExpertsC12scores$persparSTAxclass2dint <- ExpertsC12scores$persparSTA*ExpertsC12scores$class2dint

#first unrestricted (without different intercepts)
C12lm6unrestricted<-lm(ExpertsC12scores$Performance~ExpertsC12scores$persparSTA+ExpertsC12scores$persparSTAxclass2dint+ExpertsC12scores$FirmSize+ExpertsC12scores$Regulation+ExpertsC12scores$MgmtLevel+ExpertsC12scores$Pos+ExpertsC12scores$SPType+ExpertsC12scores$IndServ+ExpertsC12scores$IndUnc)

summary(C12lm6unrestricted)
C12unrestrictedSSE <- sum(C12lm6unrestricted$residuals**2)

#now restricted (we assume only the classdummy to be changing)
C12lm6restricted<-lm(ExpertsC12scores$Performance~ExpertsC12scores$persparSTA+ExpertsC12scores$FirmSize+ExpertsC12scores$Regulation+ExpertsC12scores$MgmtLevel+ExpertsC12scores$Pos+ExpertsC12scores$SPType+ExpertsC12scores$IndServ+ExpertsC12scores$IndUnc)

summary(C12lm6restricted)
C12restrictedSSE <- sum(C12lm6restricted$residuals**2)

###And the F statistic
dof <- nrow(ExpertsC12scores)-10 #degree of freedom denominator (dof for numterator is fixed to 2)
T <- ((C12restrictedSSE-C12unrestrictedSSE)/2)/(C12unrestrictedSSE/dof)
dof
T
