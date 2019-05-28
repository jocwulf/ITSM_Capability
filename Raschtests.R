#import libraries
library(eRm)

#get the csvs
data_exp_c1 <- read.csv(file="C:/Users/jwulf/Documents/IWI Forschung/GFF/ITSM Data/20161004_Experts/Experts_classSeparation_class1v2.csv", header=FALSE, sep=";")
data_exp_c1_ITIL <- data_exp_c1[,c(71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,91,92,94,95,96)]

data_exp_c2 <- read.csv(file="C:/Users/jwulf/Documents/IWI Forschung/GFF/ITSM Data/20161004_Experts/Experts_classSeparationClass2.csv", header=FALSE, sep=";")
data_exp_c2_ITIL <- data_exp_c2[,c(71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,91,92,94,95,96)]

#calculate the pcms
res_exp_c1 <- PCM(data_exp_c1_ITIL)
res_exp_c2 <- PCM(data_exp_c2_ITIL)

#distinguish industry 1
splitvecind_c1_ind1 <- ifelse(data_exp_c1[,27]>median(data_exp_c1[,27]),1,0)
LR_res_exp_c1_ind1 <- LRtest(res_exp_c1, splitcr=splitvecind_c1_ind1)

splitvecind_c2_ind1 <- ifelse(data_exp_c2[,27]>median(data_exp_c2[,27]),1,0)
LR_res_exp_c2_ind1 <- LRtest(res_exp_c2, splitcr=splitvecind_c2_ind1)

#distinguish industry 2
splitvecind_c1_ind2 <- ifelse(data_exp_c1[,28]>median(data_exp_c1[,28]),1,0)
LR_res_exp_c1_ind2 <- LRtest(res_exp_c1, splitcr=splitvecind_c1_ind2)

splitvecind_c2_ind2 <- ifelse(data_exp_c2[,28]>median(data_exp_c2[,28]),1,0)
LR_res_exp_c2_ind2 <- LRtest(res_exp_c2, splitcr=splitvecind_c2_ind2)

#and our countries DE=1
splitvecind_c1_country <- ifelse(data_exp_c1[,114]=="DE",1,0)
LR_res_exp_c1_country <- LRtest(res_exp_c1, splitcr=splitvecind_c1_country)

splitvecind_c2_country <- ifelse(data_exp_c2[,114]=="DE",1,0)
LR_res_exp_c2_country <- LRtest(res_exp_c2, splitcr=splitvecind_c2_country)

#affiliation (business vs IT)
splitvecind_c1_aff <- ifelse(data_exp_c1[,110]=="ITSP",1,0)
LR_res_exp_c1_aff <- LRtest(res_exp_c1, splitcr=splitvecind_c1_aff)

splitvecind_c2_aff <- ifelse(data_exp_c2[,110]=="ITSP",1,0)
LR_res_exp_c2_aff <- LRtest(res_exp_c2, splitcr=splitvecind_c2_aff)


#ITSP Type
splitvec_c1_type <- ifelse(data_exp_c1$V11==1,1,0)
LR_res_exp_c1_type <- LRtest(res_exp_c1, splitcr=splitvec_c1_type)

splitvec_c2_type <- ifelse(data_exp_c2$V11==1,1,0)
LR_res_exp_c2_type <- LRtest(res_exp_c2, splitcr=splitvec_c2_type)


