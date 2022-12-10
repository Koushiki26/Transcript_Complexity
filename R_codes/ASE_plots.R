#alternative splicing events plot

#human
data_lncRNA <- read.csv("lncrna_ASE_h.csv")
data_CDS <- read.csv("mrna_ASE_h.csv")

library(dplyr)
a <- sum(data_lncRNA$A3)
b <- sum(data_lncRNA$A5)
c <- sum(data_lncRNA$AF)
d <- sum(data_lncRNA$AL)
e <- sum(data_lncRNA$MX)
f <- sum(data_lncRNA$RI)
g <- sum(data_lncRNA$SE)
h <- a+b+c+d+e+f+g
a1 <- sum(data_CDS$A3)
b1 <- sum(data_CDS$A5)
c1 <- sum(data_CDS$AF)
d1 <- sum(data_CDS$AL)
e1 <- sum(data_CDS$MX)
f1 <- sum(data_CDS$RI)
g1 <- sum(data_CDS$SE)
h1 <- a1+b1+c1+d1+e1+f1+g1

library(ggplot2)
Splicing_events <- c(rep("A3",2),rep("A5",2),rep("AF",2),rep("AL",2),rep("MX",2)
                     ,rep("RI",2),rep("SE",2))
Type_of_gene <- rep(c("lncRNA","pcg"),7)
Fraction_of_transcripts <- as.numeric(c(a/h,a1/h1,b/h,b1/h1,c/h,c1/h1,d/h,d1/h1,
                                e/h,e1/h1,f/h,f1/h1,g/h,g1/h1))
df <- data.frame(Splicing_events,Type_of_gene,Fraction_of_transcripts)
df
level_order <- c("A3","A5","AF","AL","MX","RI","SE")
tiff("/home/dell/Documents/koushiki/final_fig/fig_c.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(df, aes(factor(Splicing_events, level=level_order), 
                    Fraction_of_transcripts, fill=Type_of_gene)) + 
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values = c("#4D648D", "#E8A735"))
p + xlab("Splicing events") + ylab("Fraction of transcripts")+ 
  theme(axis.title = element_text(size = 15, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black", face="bold")) + theme(legend.position = "none") 
dev.off()

#mouse
data_lncRNA <- read.csv("lncrna_ASE_m.csv")
data_CDS <- read.csv("mrna_ASE_m.csv")

library(dplyr)
a <- sum(data_lncRNA$A3)
b <- sum(data_lncRNA$A5)
c <- sum(data_lncRNA$AF)
d <- sum(data_lncRNA$AL)
e <- sum(data_lncRNA$MX)
f <- sum(data_lncRNA$RI)
g <- sum(data_lncRNA$SE)
h <- a+b+c+d+e+f+g
a1 <- sum(data_CDS$A3)
b1 <- sum(data_CDS$A5)
c1 <- sum(data_CDS$AF)
d1 <- sum(data_CDS$AL)
e1 <- sum(data_CDS$MX)
f1 <- sum(data_CDS$RI)
g1 <- sum(data_CDS$SE)
h1 <- a1+b1+c1+d1+e1+f1+g1

library(ggplot2)
Splicing_events <- c(rep("A3",2),rep("A5",2),rep("AF",2),rep("AL",2),rep("MX",2)
                     ,rep("RI",2),rep("SE",2))
Type_of_gene <- rep(c("lncRNA","pcg"),7)
Fraction_of_transcripts <- as.numeric(c(a/h,a1/h1,b/h,b1/h1,c/h,c1/h1,d/h,d1/h1,
                                        e/h,e1/h1,f/h,f1/h1,g/h,g1/h1))
df <- data.frame(Splicing_events,Type_of_gene,Fraction_of_transcripts)
df
level_order <- c("A3","A5","AF","AL","MX","RI","SE")
tiff("/home/dell/Documents/koushiki/final_fig/fig_d.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(df, aes(factor(Splicing_events, level=level_order), 
                    Fraction_of_transcripts, fill=Type_of_gene)) + 
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values = c("#4D648D", "#E8A735"))
p + xlab("Splicing events") + ylab("Fraction of transcripts")+ 
  theme(axis.title = element_text(size = 15, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black", face="bold")) + theme(legend.position = "none") 
dev.off()