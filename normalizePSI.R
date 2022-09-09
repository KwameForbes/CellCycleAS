pheatmaptester = read.csv("/nas/longleaf/home/kwamek/CellCycleAS/AS-nonorm-50counts-10perdiffPSI.csv",header = T,sep = ",")
pheatmaptester  = subset(pheatmaptester,select = c(-1))


ASeventPSI = subset(ASeventPSI,select = c(-1))


TAGS = as.data.frame(ASeventPSI[1:4])
VALS = as.data.frame(ASeventPSI[5:18])



normie <- function(x, ...) {
  (x - min(x))/max(x)
}

VALS = as.data.frame(t(VALS))
VALS = apply(VALS,2,normie)
VALS = as.data.frame(t(VALS))
ASeventPSI = cbind(TAGS,VALS)

ASeventPSI2 = as.data.frame(ASeventPSI)
ASeventPSI2[5:18] = apply(ASeventPSI2[5:18], 1, normie)
#write.csv(ASeventPSI,"AS-norm-50counts-10perdiffPSI.csv")

pheatmaptester = merge(pheatmaptester,BiggerSeed,by ="EVENT")

pheatmaptester <- pheatmaptester[order(pheatmaptester$seed,pheatmaptester$SCORE),]

# Rename column where names
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_1"] <- "TP_1"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_2"] <- "TP_2"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_3"] <- "TP_3"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_4"] <- "TP_4"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_5"] <- "TP_5"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_6"] <- "TP_6"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_7"] <- "TP_7"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_8"] <- "TP_8"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_9"] <- "TP_9"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_10"] <- "TP_10"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_11"] <- "TP_11"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_12"] <- "TP_12"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_13"] <- "TP_13"
names(pheatmaptester)[names(pheatmaptester) == "IncLevel_TimePoint_14"] <- "TP_14"




anno <- data.frame(SEED=pheatmaptester$seed,EVENT =pheatmaptester$Event,SCORE =pheatmaptester$SCORE)

rownames(pheatmaptester) <- pheatmaptester$EVENT
rownames(anno) <- pheatmaptester$EVENT

cat_df = data.frame("CellStage" = c("G1-S","S","S-G2","G2-M","M-G1","M-G1","G1","G1-S","S-G2","G2-M",
                                     "M-G1","M-G1","G1","G1-S"))
rownames(cat_df) = colnames(pheatmaptester[5:18])

pheatmaptester2 = subset(pheatmaptester,SCORE<2.5)
write.csv(pheatmaptester2,"Periodic_ALT_Splice_under3.csv")

pheatmap(pheatmaptester2[5:18], show_rownames=FALSE,cluster_cols= F,cluster_rows = F,scale = "row",
         annotation_row = anno,annotation_col = cat_df, cellwidth = 10,color=colorRampPalette(c("red", "white", "blue"))(150),
         border_color = NA)


library(ggplot2)
# Basic dot plot
p<-ggplot(pheatmaptester, aes(x=seed, y=SCORE)) + 
  geom_boxplot(binaxis='y',stackdir='center', dotsize=.1)
p




  


