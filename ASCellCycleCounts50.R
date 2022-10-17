library(dplyr)
library(tidyr)
install.packages("phateR")
library(phateR)
BiocManager::install("fCI")
library(fCI)
IncLevelvar = 0.05
FDR.variable = 0.05
counts.variable = 100




f <- function(x, ...) {
  abs(max(x) -min(x))
}

# Cellcycle_RI_JCEC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/RI.MATS.JCEC.txt",header = T)
# Cellcycle_RI_JCEC <- as.data.frame(Cellcycle_RI_JCEC)
# Cellcycle_RI_JCEC$Key <- paste(Cellcycle_RI_JCEC$ID,Cellcycle_RI_JCEC$GeneID,sep = "_")
# rownames(Cellcycle_RI_JCEC) <- Cellcycle_RI_JCEC$Key
# AS = as.data.frame(rep("RI",dim(Cellcycle_RI_JCEC)[1]))
# names(AS)[names(AS) == 'rep("RI", dim(Cellcycle_RI_JCEC)[1])'] <- "Event"
# Cellcycle_RI_JCEC = cbind(Cellcycle_RI_JCEC,AS)
# 
# 
# Cellcycle_RI_JCEC <- separate(Cellcycle_RI_JCEC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_RI_JCEC <- separate(Cellcycle_RI_JCEC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_RI_JCEC <- separate(Cellcycle_RI_JCEC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_RI_JCEC <- separate(Cellcycle_RI_JCEC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_RI_JCEC[11:40] <- as.numeric(unlist(Cellcycle_RI_JCEC[11:40]))
# 
# Cellcycle_RI_JCEC <- separate(Cellcycle_RI_JCEC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_RI_JCEC <- separate(Cellcycle_RI_JCEC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_RI_JCEC[45:58] <- as.numeric(unlist(Cellcycle_RI_JCEC[45:58]))
# 
# Cellcycle_RI_JCEC = unite(Cellcycle_RI_JCEC,"Position",c(4:12),sep = ":",remove = T)
# 
# 
# Cellcycle_RI_JCEC50 = subset(Cellcycle_RI_JCEC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
#                                IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
#                                IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)
# 
# 
# Cellcycle_RI_JCEC50$PSIDiff <- (apply(Cellcycle_RI_JCEC50[37:50], MARGIN =  1, f,))
# 
# 
# 
# 
# # Cellcycle_RIsigs_JCEC = subset(Cellcycle_RI_JCEC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
# #                                  IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
# 


Cellcycle_RI_JC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/RI.MATS.JC.txt",header = T)
Cellcycle_RI_JC <- as.data.frame(Cellcycle_RI_JC)
Cellcycle_RI_JC$Key <- paste(Cellcycle_RI_JC$ID,Cellcycle_RI_JC$GeneID,sep = "_")
rownames(Cellcycle_RI_JC) <- Cellcycle_RI_JC$Key
AS = as.data.frame(rep("RI",dim(Cellcycle_RI_JC)[1]))
names(AS)[names(AS) == 'rep("RI", dim(Cellcycle_RI_JC)[1])'] <- "Event"
Cellcycle_RI_JC = cbind(Cellcycle_RI_JC,AS)

Cellcycle_RI_JC <- separate(Cellcycle_RI_JC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_RI_JC <- separate(Cellcycle_RI_JC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_RI_JC <- separate(Cellcycle_RI_JC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_RI_JC <- separate(Cellcycle_RI_JC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_RI_JC[11:40] <- as.numeric(unlist(Cellcycle_RI_JC[11:40]))

Cellcycle_RI_JC <- separate(Cellcycle_RI_JC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
Cellcycle_RI_JC <- separate(Cellcycle_RI_JC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
Cellcycle_RI_JC[45:58] <- as.numeric(unlist(Cellcycle_RI_JC[45:58]))

Cellcycle_RI_JC = unite(Cellcycle_RI_JC,"Position",c(4:11),sep = ":",remove = T)

Cellcycle_RI_JC50 = subset(Cellcycle_RI_JC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
                               IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
                               IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)


Cellcycle_RI_JC50$PSIDiff <- (apply(Cellcycle_RI_JC50[38:51], MARGIN =  1, f,))


# Cellcycle_RIsigs_JC = subset(Cellcycle_RI_JC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
#                                  IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cellcycle_SE_JCEC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/SE.MATS.JCEC.txt",header = T)
# Cellcycle_SE_JCEC <- as.data.frame(Cellcycle_SE_JCEC)
# Cellcycle_SE_JCEC$Key <- paste(Cellcycle_SE_JCEC$ID,Cellcycle_SE_JCEC$GeneID,sep = "_")
# rownames(Cellcycle_SE_JCEC) <- Cellcycle_SE_JCEC$Key
# AS = as.data.frame(rep("SE",dim(Cellcycle_SE_JCEC)[1]))
# names(AS)[names(AS) == 'rep("SE", dim(Cellcycle_SE_JCEC)[1])'] <- "Event"
# Cellcycle_SE_JCEC = cbind(Cellcycle_SE_JCEC,AS)
# 
# Cellcycle_SE_JCEC <- separate(Cellcycle_SE_JCEC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_SE_JCEC <- separate(Cellcycle_SE_JCEC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_SE_JCEC <- separate(Cellcycle_SE_JCEC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_SE_JCEC <- separate(Cellcycle_SE_JCEC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_SE_JCEC[11:40] <- as.numeric(unlist(Cellcycle_SE_JCEC[11:40]))
# 
# Cellcycle_SE_JCEC <- separate(Cellcycle_SE_JCEC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_SE_JCEC <- separate(Cellcycle_SE_JCEC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_SE_JCEC[45:58] <- as.numeric(unlist(Cellcycle_SE_JCEC[45:58]))
# 
# Cellcycle_SE_JCEC = unite(Cellcycle_SE_JCEC,"Position",c(4:12),sep = ":",remove = T)
# 
# Cellcycle_SE_JCEC50 = subset(Cellcycle_SE_JCEC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
#                                IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
#                                IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)
# 
# 
# Cellcycle_SE_JCEC50$PSIDiff <- (apply(Cellcycle_SE_JCEC50[37:50], MARGIN =  1, f,))
# 
# 
# 
# # Cellcycle_SEsigs_JCEC = subset(Cellcycle_SE_JCEC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
# #                                  IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
# # 


Cellcycle_SE_JC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/SE.MATS.JC.txt",header = T)
Cellcycle_SE_JC <- as.data.frame(Cellcycle_SE_JC)
Cellcycle_SE_JC$Key <- paste(Cellcycle_SE_JC$ID,Cellcycle_SE_JC$GeneID,sep = "_")
rownames(Cellcycle_SE_JC) <- Cellcycle_SE_JC$Key
AS = as.data.frame(rep("SE",dim(Cellcycle_SE_JC)[1]))
names(AS)[names(AS) == 'rep("SE", dim(Cellcycle_SE_JC)[1])'] <- "Event"
Cellcycle_SE_JC = cbind(Cellcycle_SE_JC,AS)

Cellcycle_SE_JC <- separate(Cellcycle_SE_JC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_SE_JC <- separate(Cellcycle_SE_JC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_SE_JC <- separate(Cellcycle_SE_JC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_SE_JC <- separate(Cellcycle_SE_JC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_SE_JC[11:40] <- as.numeric(unlist(Cellcycle_SE_JC[11:40]))

Cellcycle_SE_JC <- separate(Cellcycle_SE_JC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
Cellcycle_SE_JC <- separate(Cellcycle_SE_JC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
Cellcycle_SE_JC[45:58] <- as.numeric(unlist(Cellcycle_SE_JC[45:58]))

Cellcycle_SE_JC = unite(Cellcycle_SE_JC,"Position",c(4:11),sep = ":",remove = T)

Cellcycle_SE_JC50 = subset(Cellcycle_SE_JC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
                             IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
                             IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)


Cellcycle_SE_JC50$PSIDiff <- (apply(Cellcycle_SE_JC50[38:51], MARGIN =  1, f,))




# Cellcycle_SEsigs_JC = subset(Cellcycle_SE_JC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
#                                IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cellcycle_A3SS_JCEC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/A3SS.MATS.JCEC.txt",header = T)
# Cellcycle_A3SS_JCEC <- as.data.frame(Cellcycle_A3SS_JCEC)
# Cellcycle_A3SS_JCEC$Key <- paste(Cellcycle_A3SS_JCEC$ID,Cellcycle_A3SS_JCEC$GeneID,sep = "_")
# rownames(Cellcycle_A3SS_JCEC) <- Cellcycle_A3SS_JCEC$Key
# AS = as.data.frame(rep("A3SS",dim(Cellcycle_A3SS_JCEC)[1]))
# names(AS)[names(AS) == 'rep("A3SS", dim(Cellcycle_A3SS_JCEC)[1])'] <- "Event"
# Cellcycle_A3SS_JCEC = cbind(Cellcycle_A3SS_JCEC,AS)
# 
# Cellcycle_A3SS_JCEC <- separate(Cellcycle_A3SS_JCEC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_A3SS_JCEC <- separate(Cellcycle_A3SS_JCEC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_A3SS_JCEC <- separate(Cellcycle_A3SS_JCEC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_A3SS_JCEC <- separate(Cellcycle_A3SS_JCEC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_A3SS_JCEC[11:40] <- as.numeric(unlist(Cellcycle_A3SS_JCEC[11:40]))
# 
# Cellcycle_A3SS_JCEC <- separate(Cellcycle_A3SS_JCEC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_A3SS_JCEC <- separate(Cellcycle_A3SS_JCEC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_A3SS_JCEC[45:58] <- as.numeric(unlist(Cellcycle_A3SS_JCEC[45:58]))
# 
# Cellcycle_A3SS_JCEC = unite(Cellcycle_A3SS_JCEC,"Position",c(4:12),sep = ":",remove = T)
# 
# Cellcycle_A3SS_JCEC50 = subset(Cellcycle_A3SS_JCEC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
#                                IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
#                                IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)
# 
# 
# Cellcycle_A3SS_JCEC50$PSIDiff <- (apply(Cellcycle_A3SS_JCEC50[37:50], MARGIN =  1, f,))
# 
# 
# # Cellcycle_A3SSsigs_JCEC = subset(Cellcycle_A3SS_JCEC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
# #                                  IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
# 


Cellcycle_A3SS_JC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/A3SS.MATS.JC.txt",header = T)
Cellcycle_A3SS_JC <- as.data.frame(Cellcycle_A3SS_JC)
Cellcycle_A3SS_JC$Key <- paste(Cellcycle_A3SS_JC$ID,Cellcycle_A3SS_JC$GeneID,sep = "_")
rownames(Cellcycle_A3SS_JC) <- Cellcycle_A3SS_JC$Key
AS = as.data.frame(rep("A3SS",dim(Cellcycle_A3SS_JC)[1]))
names(AS)[names(AS) == 'rep("A3SS", dim(Cellcycle_A3SS_JC)[1])'] <- "Event"
Cellcycle_A3SS_JC = cbind(Cellcycle_A3SS_JC,AS)

Cellcycle_A3SS_JC <- separate(Cellcycle_A3SS_JC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_A3SS_JC <- separate(Cellcycle_A3SS_JC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_A3SS_JC <- separate(Cellcycle_A3SS_JC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_A3SS_JC <- separate(Cellcycle_A3SS_JC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_A3SS_JC[11:40] <- as.numeric(unlist(Cellcycle_A3SS_JC[11:40]))

Cellcycle_A3SS_JC <- separate(Cellcycle_A3SS_JC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
Cellcycle_A3SS_JC <- separate(Cellcycle_A3SS_JC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
Cellcycle_A3SS_JC[45:58] <- as.numeric(unlist(Cellcycle_A3SS_JC[45:58]))

Cellcycle_A3SS_JC = unite(Cellcycle_A3SS_JC,"Position",c(4:11),sep = ":",remove = T)

Cellcycle_A3SS_JC50 = subset(Cellcycle_A3SS_JC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
                             IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
                             IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)


Cellcycle_A3SS_JC50$PSIDiff <- (apply(Cellcycle_A3SS_JC50[38:51], MARGIN =  1, f,))




# Cellcycle_A3SSsigs_JC = subset(Cellcycle_A3SS_JC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
#                                IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cellcycle_A5SS_JCEC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/A5SS.MATS.JCEC.txt",header = T)
# Cellcycle_A5SS_JCEC <- as.data.frame(Cellcycle_A5SS_JCEC)
# Cellcycle_A5SS_JCEC$Key <- paste(Cellcycle_A5SS_JCEC$ID,Cellcycle_A5SS_JCEC$GeneID,sep = "_")
# rownames(Cellcycle_A5SS_JCEC) <- Cellcycle_A5SS_JCEC$Key
# AS = as.data.frame(rep("A5SS",dim(Cellcycle_A5SS_JCEC)[1]))
# names(AS)[names(AS) == 'rep("A5SS", dim(Cellcycle_A5SS_JCEC)[1])'] <- "Event"
# Cellcycle_A5SS_JCEC = cbind(Cellcycle_A5SS_JCEC,AS)
# 
# Cellcycle_A5SS_JCEC <- separate(Cellcycle_A5SS_JCEC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_A5SS_JCEC <- separate(Cellcycle_A5SS_JCEC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_A5SS_JCEC <- separate(Cellcycle_A5SS_JCEC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_A5SS_JCEC <- separate(Cellcycle_A5SS_JCEC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_A5SS_JCEC[11:40] <- as.numeric(unlist(Cellcycle_A5SS_JCEC[11:40]))
# 
# Cellcycle_A5SS_JCEC <- separate(Cellcycle_A5SS_JCEC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_A5SS_JCEC <- separate(Cellcycle_A5SS_JCEC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_A5SS_JCEC[45:58] <- as.numeric(unlist(Cellcycle_A5SS_JCEC[45:58]))
# 
# Cellcycle_A5SS_JCEC = unite(Cellcycle_A5SS_JCEC,"Position",c(4:12),sep = ":",remove = T)
# 
# Cellcycle_A5SS_JCEC50 = subset(Cellcycle_A5SS_JCEC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
#                                IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
#                                IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)
# 
# 
# Cellcycle_A5SS_JCEC50$PSIDiff <- (apply(Cellcycle_A5SS_JCEC50[37:50], MARGIN =  1, f,))
# 
# 
# 
# # Cellcycle_A5SSsigs_JCEC = subset(Cellcycle_A5SS_JCEC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
# #                                  IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
# 


Cellcycle_A5SS_JC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/A5SS.MATS.JC.txt",header = T)
Cellcycle_A5SS_JC <- as.data.frame(Cellcycle_A5SS_JC)
Cellcycle_A5SS_JC$Key <- paste(Cellcycle_A5SS_JC$ID,Cellcycle_A5SS_JC$GeneID,sep = "_")
rownames(Cellcycle_A5SS_JC) <- Cellcycle_A5SS_JC$Key
AS = as.data.frame(rep("A5SS",dim(Cellcycle_A5SS_JC)[1]))
names(AS)[names(AS) == 'rep("A5SS", dim(Cellcycle_A5SS_JC)[1])'] <- "Event"
Cellcycle_A5SS_JC = cbind(Cellcycle_A5SS_JC,AS)

Cellcycle_A5SS_JC <- separate(Cellcycle_A5SS_JC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_A5SS_JC <- separate(Cellcycle_A5SS_JC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_A5SS_JC <- separate(Cellcycle_A5SS_JC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_A5SS_JC <- separate(Cellcycle_A5SS_JC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_A5SS_JC[11:40] <- as.numeric(unlist(Cellcycle_A5SS_JC[11:40]))

Cellcycle_A5SS_JC <- separate(Cellcycle_A5SS_JC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
Cellcycle_A5SS_JC <- separate(Cellcycle_A5SS_JC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
Cellcycle_A5SS_JC[45:58] <- as.numeric(unlist(Cellcycle_A5SS_JC[45:58]))

Cellcycle_A5SS_JC = unite(Cellcycle_A5SS_JC,"Position",c(4:11),sep = ":",remove = T)

Cellcycle_A5SS_JC50 = subset(Cellcycle_A5SS_JC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
                             IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
                             IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)


Cellcycle_A5SS_JC50$PSIDiff <- (apply(Cellcycle_A5SS_JC50[38:51], MARGIN =  1, f,))




# Cellcycle_A5SSsigs_JC = subset(Cellcycle_A5SS_JC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
#                                IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cellcycle_MXE_JCEC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/MXE.MATS.JCEC.txt",header = T)
# Cellcycle_MXE_JCEC <- as.data.frame(Cellcycle_MXE_JCEC)
# Cellcycle_MXE_JCEC$Key <- paste(Cellcycle_MXE_JCEC$ID,Cellcycle_MXE_JCEC$GeneID,sep = "_")
# rownames(Cellcycle_MXE_JCEC) <- Cellcycle_MXE_JCEC$Key
# AS = as.data.frame(rep("MXE",dim(Cellcycle_MXE_JCEC)[1]))
# names(AS)[names(AS) == 'rep("MXE", dim(Cellcycle_MXE_JCEC)[1])'] <- "Event"
# Cellcycle_MXE_JCEC = cbind(Cellcycle_MXE_JCEC,AS)
# 
# Cellcycle_MXE_JCEC <- separate(Cellcycle_MXE_JCEC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_MXE_JCEC <- separate(Cellcycle_MXE_JCEC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_MXE_JCEC <- separate(Cellcycle_MXE_JCEC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_MXE_JCEC <- separate(Cellcycle_MXE_JCEC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_MXE_JCEC[15:42] <- as.numeric(unlist(Cellcycle_MXE_JCEC[15:42]))
# 
# Cellcycle_MXE_JCEC <- separate(Cellcycle_MXE_JCEC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
# Cellcycle_MXE_JCEC <- separate(Cellcycle_MXE_JCEC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
# Cellcycle_MXE_JCEC[47:60] <- as.numeric(unlist(Cellcycle_MXE_JCEC[47:60]))
# 
# Cellcycle_MXE_JCEC = unite(Cellcycle_MXE_JCEC,"Position",c(4:14),sep = ":",remove = T)
# 
# Cellcycle_MXE_JCEC50 = subset(Cellcycle_MXE_JCEC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
#                                IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
#                                IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)
# 
# 
# Cellcycle_MXE_JCEC50$PSIDiff <- (apply(Cellcycle_MXE_JCEC50[37:50], MARGIN =  1, f,))
# 
# 
# 
# # Cellcycle_MXEsigs_JCEC = subset(Cellcycle_MXE_JCEC, IncLevel_TimePoint_1 & IncLevel_TimePoint_2 & IncLevel_TimePoint_3 & IncLevel_TimePoint_4 & IncLevel_TimePoint_5 & IncLevel_TimePoint_6 & IncLevel_TimePoint_7 & IncLevel_TimePoint_8 &
# #                                  IncLevel_TimePoint_9 & IncLevel_TimePoint_10 & IncLevel_TimePoint_11 & IncLevel_TimePoint_12 & IncLevel_TimePoint_13 & IncLevel_TimePoint_14 > IncLevelvar)
# # 
# 

Cellcycle_MXE_JC <- read.table("/pine/scr/k/w/kwamek/pengda_collab/RMATS-Cellcycle/MXE.MATS.JC.txt",header = T)
Cellcycle_MXE_JC <- as.data.frame(Cellcycle_MXE_JC)
Cellcycle_MXE_JC$Key <- paste(Cellcycle_MXE_JC$ID,Cellcycle_MXE_JC$GeneID,sep = "_")
rownames(Cellcycle_MXE_JC) <- Cellcycle_MXE_JC$Key
AS = as.data.frame(rep("MXE",dim(Cellcycle_MXE_JC)[1]))
names(AS)[names(AS) == 'rep("MXE", dim(Cellcycle_MXE_JC)[1])'] <- "Event"
Cellcycle_MXE_JC = cbind(Cellcycle_MXE_JC,AS)


Cellcycle_MXE_JC <- separate(Cellcycle_MXE_JC, col = IJC_SAMPLE_1, into = c("IJC_TimePoint_1","IJC_TimePoint_2","IJC_TimePoint_3","IJC_TimePoint_6","IJC_TimePoint_7","IJC_TimePoint_8","IJC_TimePoint_12","IJC_TimePoint_13","IJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_MXE_JC <- separate(Cellcycle_MXE_JC, col = SJC_SAMPLE_1, into = c("SJC_TimePoint_1","SJC_TimePoint_2","SJC_TimePoint_3","SJC_TimePoint_6","SJC_TimePoint_7","SJC_TimePoint_8","SJC_TimePoint_12","SJC_TimePoint_13","SJC_TimePoint_14"),sep = ",",remove = T)
Cellcycle_MXE_JC <- separate(Cellcycle_MXE_JC, col = IJC_SAMPLE_2, into = c("IJC_TimePoint_4","IJC_TimePoint_5","IJC_TimePoint_9","IJC_TimePoint_10","IJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_MXE_JC <- separate(Cellcycle_MXE_JC, col = SJC_SAMPLE_2, into = c("SJC_TimePoint_4","SJC_TimePoint_5","SJC_TimePoint_9","SJC_TimePoint_10","SJC_TimePoint_11"),sep = ",",remove = T)
Cellcycle_MXE_JC[15:42] <- as.numeric(unlist(Cellcycle_MXE_JC[15:42]))

Cellcycle_MXE_JC <- separate(Cellcycle_MXE_JC, col = IncLevel1, into = c("IncLevel_TimePoint_1","IncLevel_TimePoint_2","IncLevel_TimePoint_3","IncLevel_TimePoint_6","IncLevel_TimePoint_7","IncLevel_TimePoint_8","IncLevel_TimePoint_12","IncLevel_TimePoint_13","IncLevel_TimePoint_14"),sep = ",",remove = T)
Cellcycle_MXE_JC <- separate(Cellcycle_MXE_JC, col = IncLevel2, into = c("IncLevel_TimePoint_4","IncLevel_TimePoint_5","IncLevel_TimePoint_9","IncLevel_TimePoint_10","IncLevel_TimePoint_11"),sep = ",",remove = T)
Cellcycle_MXE_JC[47:60] <- as.numeric(unlist(Cellcycle_MXE_JC[47:60]))

Cellcycle_MXE_JC = unite(Cellcycle_MXE_JC,"Position",c(4:13),sep = ":",remove = T)

Cellcycle_MXE_JC50 = subset(Cellcycle_MXE_JC,IJC_TimePoint_1+SJC_TimePoint_1&IJC_TimePoint_2+SJC_TimePoint_2&IJC_TimePoint_3+SJC_TimePoint_3&IJC_TimePoint_4+SJC_TimePoint_4&IJC_TimePoint_5+SJC_TimePoint_5&IJC_TimePoint_6+SJC_TimePoint_6&
                             IJC_TimePoint_7+SJC_TimePoint_7&IJC_TimePoint_8+SJC_TimePoint_8&IJC_TimePoint_9+SJC_TimePoint_9&IJC_TimePoint_10+SJC_TimePoint_10&IJC_TimePoint_11+SJC_TimePoint_11&IJC_TimePoint_12+SJC_TimePoint_12&
                             IJC_TimePoint_13+SJC_TimePoint_13&IJC_TimePoint_14+SJC_TimePoint_14>=counts.variable)


Cellcycle_MXE_JC50$PSIDiff <- (apply(Cellcycle_MXE_JC50[38:51], MARGIN =  1, f,))



ASevents50 <- rbind(Cellcycle_RI_JC50,Cellcycle_SE_JC50,Cellcycle_A3SS_JC50,Cellcycle_A5SS_JC50,Cellcycle_MXE_JC50)
ASevents50 <- na.omit(ASevents50)
ASevents50Diff = subset(ASevents50,PSIDiff>=0.1)

ASevents50Diff = unite(ASevents50Diff,"EVENT",c(3,4),sep = "_", remove = F)
ASevents50Diffclean = subset(ASevents50Diff,select = c(55,3,2,39:41,48:49,42:44,50:52,45:47))
rownames(ASevents50Diffclean) <- NULL
#write.csv(ASevents50Diffclean,"AS-nonorm-50counts-10perdiffPSI.csv")



