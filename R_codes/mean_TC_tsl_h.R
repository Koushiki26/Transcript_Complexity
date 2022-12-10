#mean transcript complexity for different tsl for human

#tsl1
data_lncRNA <- read.csv("lncrna_tsl1_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS <- read.csv("mrna_tsl1_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_2a1.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.5,4,2,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA$TPE, data_CDS$TPE, col = c("#4D648D", "#E8A735"),ylab="TC", boxwex=2.0,
        xaxs="i",at=c(1,4), xaxt='n',
        outline=F, xlim=c(0,5), ylim=c(0.0,1), main="TSL1")
dev.off()

#tsl1+tsl2
data_lncRNA <- read.csv("lncrna_tsl2_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS <- read.csv("mrna_tsl2_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_2a2.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.5,4,2,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA$TPE, data_CDS$TPE, col = c("#4D648D", "#E8A735"),ylab="TC", boxwex=2.0,
        xaxs="i",at=c(1,4), xaxt='n',
        outline=F, xlim=c(0,5), ylim=c(0.0,1), main="TSL1+TSL2")
dev.off()

#tsl1+tsl2+tsl3
data_lncRNA <- read.csv("lncrna_tsl3_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS <- read.csv("mrna_tsl3_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_2a3.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.5,4,2,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA$TPE, data_CDS$TPE, col = c("#4D648D", "#E8A735"),ylab="TC", boxwex=2.0,
        xaxs="i",at=c(1,4), xaxt='n',
        outline=F, xlim=c(0,5), ylim=c(0.0,1), main="TSL1+TSL2+TSL3")
dev.off()

#tsl1+tsl2+tsl3+tsl4
data_lncRNA <- read.csv("lncrna_tsl4_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS <- read.csv("mrna_tsl4_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_2a4.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.5,4,2,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA$TPE, data_CDS$TPE, col = c("#4D648D", "#E8A735"),ylab="TC", boxwex=2.0,
        xaxs="i",at=c(1,4), xaxt='n',
        outline=F, xlim=c(0,5), ylim=c(0.0,1), main="TSL1+TSL2+TSL3+TSL4")
dev.off()

#tsl1+tsl2+tsl3+tsl4+tsl5
data_lncRNA <- read.csv("lncrna_tsl5_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS <- read.csv("mrna_tsl5_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_2a5.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.5,4,2,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA$TPE, data_CDS$TPE, col = c("#4D648D", "#E8A735"),ylab="TC", boxwex=2.0,
        xaxs="i",at=c(1,4), xaxt='n',
        outline=F, xlim=c(0,5), ylim=c(0.0,1), main="TSL1+TSL2+TSL3+TSL4+TSL5")
dev.off()

#tsl1+tsl2+tsl3+tsl4+tsl5+tslNA
data_lncRNA <- read.csv("lncrna_tsl6_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS <- read.csv("mrna_tsl6_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_2a6.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.5,4,2,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA$TPE, data_CDS$TPE, col = c("#4D648D", "#E8A735"),ylab="TC", boxwex=2.0,
        xaxs="i",at=c(1,4), xaxt='n',cex.main=1.05,font.main=2,
        outline=F, xlim=c(0,5), ylim=c(0.0,1), main="TSL1+TSL2+TSL3+TSL4+TSL5+TSLNA")
dev.off()