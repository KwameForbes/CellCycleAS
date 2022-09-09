library(pheatmap)

AS = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/AS-nonorm-50counts-10perdiffPSI.csv",header = T,sep = ",")
AS = subset(AS,select = c(-1))

Score1 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score1.csv",sep = ",",header = T)
Score1 = subset(Score1,select = c(-1))
names(Score1)[names(Score1) == 'SCORE'] <- 'SEED1'
names(Score1)[names(Score1) == 'FDR'] <- 'FDR1'

Score2 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score2.csv",sep = ",",header = T)
Score2 = subset(Score2,select = c(-1))
names(Score2)[names(Score2) == 'SCORE'] <- 'SEED2'
names(Score2)[names(Score2) == 'FDR'] <- 'FDR2'

Score3 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score3.csv",sep = ",",header = T)
Score3 = subset(Score3,select = c(-1))
names(Score3)[names(Score3) == 'SCORE'] <- 'SEED3'
names(Score3)[names(Score3) == 'FDR'] <- 'FDR3'

Score4 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score4.csv",sep = ",",header = T)
Score4 = subset(Score4,select = c(-1))
names(Score4)[names(Score4) == 'SCORE'] <- 'SEED4'
names(Score4)[names(Score4) == 'FDR'] <- 'FDR4'

Score5 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score5.csv",sep = ",",header = T)
Score5 = subset(Score5,select = c(-1))
names(Score5)[names(Score5) == 'SCORE'] <- 'SEED5'
names(Score5)[names(Score5) == 'FDR'] <- 'FDR5'

Score6 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score6.csv",sep = ",",header = T)
Score6 = subset(Score6,select = c(-1))
names(Score6)[names(Score6) == 'SCORE'] <- 'SEED6'
names(Score6)[names(Score6) == 'FDR'] <- 'FDR6'

Score7 = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/ASscores/score7.csv",sep = ",",header = T)
Score7 = subset(Score7,select = c(-1))
names(Score7)[names(Score7) == 'SCORE'] <- 'SEED7'
names(Score7)[names(Score7) == 'FDR'] <- 'FDR7'


BigSeed = merge(Score1,Score2,by="EVENT")
BigSeed = merge(BigSeed,Score3,by="EVENT")
BigSeed = merge(BigSeed,Score4,by="EVENT")
BigSeed = merge(BigSeed,Score5,by="EVENT")
BigSeed = merge(BigSeed,Score6,by="EVENT")
BigSeed = merge(BigSeed,Score7,by="EVENT")

ASandSeed = merge(AS,BigSeed,by="EVENT")

pheatmap(ASandSeed[5:18], show_rownames=FALSE,cluster_cols= FALSE,cluster_rows = FALSE,annotation_row = GeneID)

Score1 <- subset(Score1, FDR1 <=0.25)
Score2 <- subset(Score2, FDR2 <=0.25)
Score3 <- subset(Score3, FDR3 <=0.25)
Score4 <- subset(Score4, FDR4 <=0.25)
Score5 <- subset(Score5, FDR5 <=0.25)
Score6 <- subset(Score6, FDR6 <=0.25)
Score7 <- subset(Score7, FDR7 <=0.25)

BiggerSeed <- merge(Score1,Score2, by="EVENT", all.x = T,all.y = T)
BiggerSeed <- merge(BiggerSeed,Score3,by="EVENT", all.x = T,all.y = T)
BiggerSeed <- merge(BiggerSeed,Score4,by="EVENT", all.x = T,all.y = T)
BiggerSeed <- merge(BiggerSeed,Score5,by="EVENT", all.x = T,all.y = T)
BiggerSeed <- merge(BiggerSeed,Score6,by="EVENT", all.x = T,all.y = T)
BiggerSeed <- merge(BiggerSeed,Score7,by="EVENT", all.x = T,all.y = T)
BiggerSeed <- subset(BiggerSeed, select = c(1,2,4,6,8,10,12,14,3,5,7,9,11,13,15))

BiggerSeed[is.na(BiggerSeed)] = 1000
BiggerSeed$SCORE <- apply(BiggerSeed[,c(2:8)], 1, FUN=min)
BiggerSeed$seed <- colnames(BiggerSeed[,c(2:8)])[apply(BiggerSeed[,c(2:8)], 1, FUN=which.min)] 

rownames(BiggerSeed) = NULL
BiggerSeed = subset(BiggerSeed,select=c(1,16,17))
BiggerSeed <- BiggerSeed[order(BiggerSeed$SCORE),]
BiggerSeed = subset(BiggerSeed,SCORE < 999)
# write.csv(BiggerSeed,"BiggerSeed.csv")

ASandSeed = merge(AS,BiggerSeed,by="EVENT")
write.csv(ASandSeed,"norm-PSI-and-best-Seed.csv")


anno <- data.frame(SEED=ASandSeed$seed,SCORE=ASandSeed$SCORE,EVENT =ASandSeed$Event)

rownames(ASandSeed) <- ASandSeed$EVENT
rownames(anno) <-  ASandSeed$EVENT

cat_df = data.frame("CellStage" = c("G1-S","S","S-G2","G2-M","M-G1","M-G1","G1","G1-S","S-G2","G2-M",
                                    "M-G1","M-G1","G1","G1-S"))
rownames(cat_df) = colnames(ASandSeed[5:18])

ASandSeed <- ASandSeed[order(ASandSeed$seed,ASandSeed$SCORE),]

pheatmap(ASandSeed[5:18], show_rownames=F,cluster_cols= F,cluster_rows = F,scale = 'row',
         annotation_row = anno,annotation_col = cat_df, cellwidth = 10,color=colorRampPalette(c("red", "white", "blue"))(150),
         border_color = NA)

