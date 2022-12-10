#plot of splicesite

table1 <- read.csv("sscount_low_TC_m.csv", sep = " ")
head(table1)
table1$lncRNA_count5 <- as.numeric(table1$lncRNA_count5)
table1$CDS_count5<- as.numeric(table1$CDS_count5)
table1$lncRNA_count3 <- as.numeric(table1$lncRNA_count3)
table1$CDS_count3 <- as.numeric(table1$CDS_count3)
table1

table2 <- read.csv("sscount_high_TC_m.csv", sep = " ")
head(table2)
table2$lncRNA_count5 <- as.numeric(table2$lncRNA_count5)
table2$CDS_count5<- as.numeric(table2$CDS_count5)
table2$lncRNA_count3 <- as.numeric(table2$lncRNA_count3)
table2$CDS_count3 <- as.numeric(table2$CDS_count3)
table2

lncRNA_t <- sum(table1$lncRNA_count5,table2$lncRNA_count5)
CDS_t <- sum(table1$CDS_count5,table2$CDS_count5)

per11 <- (table1[4,2]/lncRNA_t)*100
per12 <- (table2[4,2]/lncRNA_t)*100
per13 <- (table1[4,3]/CDS_t)*100
per14 <- (table2[4,3]/CDS_t)*100
per21 <- (table1[10,2]/lncRNA_t)*100
per22 <- (table2[10,2]/lncRNA_t)*100
per23 <- (table1[10,3]/CDS_t)*100
per24 <- (table2[10,3]/CDS_t)*100
per31 <- (table1[12,2]/lncRNA_t)*100
per32 <- (table2[12,2]/lncRNA_t)*100
per33 <- (table1[12,3]/CDS_t)*100
per34 <- (table2[12,3]/CDS_t)*100

ss <- c(rep("AT",4),rep("GC",4),rep("GT",4))
Type_of_gene <- rep(c("lncRNA (less TPE)","lncRNA (more TPE)","pcg (less TPE)","pcg (more TPE)"),3)
Frequency_of_introns <- c(per11,per12,per13,per14,per21,per22,per23,per24,per31,per32,per33,per34)
df <- data.frame(ss,Type_of_gene,Frequency_of_introns)
df
level_order <- c("AT","GC","GT")
tiff("/home/dell/Documents/koushiki/final_fig/fig_2b7.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(df, aes(factor(ss, level=level_order), Frequency_of_introns, fill=Type_of_gene)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = c("#31A2AC","#3A5199","#EB8A3E","#CB0000"))
p + xlab("") + ylab("Fraction of introns")+ 
  theme(axis.title = element_text(size = 19, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face="bold")) + theme(legend.position = "none")
dev.off()


per11 <- (table1[2,4]/lncRNA_t)*100
per12 <- (table2[2,4]/lncRNA_t)*100
per13 <- (table1[2,5]/CDS_t)*100
per14 <- (table2[2,5]/CDS_t)*100
per21 <- (table1[3,4]/lncRNA_t)*100
per22 <- (table2[3,4]/lncRNA_t)*100
per23 <- (table1[3,5]/CDS_t)*100
per24 <- (table2[3,5]/CDS_t)*100

ss <- c(rep("AC",4),rep("AG",4))
Type_of_gene <- rep(c("lncRNA (less TPE)","lncRNA (more TPE)","pcg (less TPE)","pcg (more TPE)"),2)
Frequency_of_introns <- c(per11,per12,per13,per14,per21,per22,per23,per24)
df <- data.frame(ss,Type_of_gene,Frequency_of_introns)
df
level_order <- c("AC","AG")
tiff("/home/dell/Documents/koushiki/final_fig/fig_2b8.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(df, aes(factor(ss, level=level_order), Frequency_of_introns, fill=Type_of_gene)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = c("#31A2AC","#3A5199","#EB8A3E","#CB0000"))
p + xlab("") + ylab("Fraction of introns")+ 
  theme(axis.title = element_text(size = 19, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face="bold")) + theme(legend.position = "none")
dev.off()

