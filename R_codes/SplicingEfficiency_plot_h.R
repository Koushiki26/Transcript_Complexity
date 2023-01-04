#splicing efficiency plot

library(dplyr)

x <- 0.333333
lncRNA_data1 <- read.csv("lncRNA_K562.csv", header = T)
head(lncRNA_data1)
cor.test(lncRNA_data1$score, lncRNA_data1$TC, method = "spearman")
a <- lncRNA_data1 %>% filter(TC <= x)
count(a)
a$type <- "lncRNA (low TC)"
b <- lncRNA_data1 %>% filter(TC > x)
count(b)
b$type <- "lncRNA (high TC)"
summary(b)
mRNA_data1 <- read.csv("mRNA_K562.csv", header = T)
head(mRNA_data1)
cor.test(mRNA_data1$score, mRNA_data1$TC, method = "spearman")
c <- mRNA_data1 %>% filter(TC <= x)
count(c)
c$type <- "mRNA (low TC)"
d <- mRNA_data1 %>% filter(TC > x)
count(d)
d$type <- "mRNA (high TC)"
total_data <- data.frame(rbind(a,b,c,d))
median(lncRNA_data1$score)
median(mRNA_data1$score)
ks.test(lncRNA_data1$score, mRNA_data1$score, alternative = "two.sided")

install.packages("ggridges")
library(ggridges)
library(ggplot2)

tiff("/home/dell/Documents/koushiki/final_fig/SE_ridge1.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(total_data, aes(x = score, y = type, fill = type)) +
  geom_density_ridges() + theme_ridges() + theme(legend.position = "none") + 
  scale_fill_manual(values = c("#3A5199","#31A2AC","#CB0000","#EB8A3E"))
p + xlab("Splicing Efficiency") + ylab("")+ xlim(c(-0.4,1.2)) +
  theme(plot.margin = margin(0.25, -0.5, 0.25, -2.75, "cm"),
        axis.title.x = element_text(size = 19, color = "black", face = "bold", hjust=0.6),
        axis.text = element_text(size = 16, color = "black", face="bold"),
        plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.6),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ggtitle("K562 Cell Line")
dev.off()


x <- 0.333333
lncRNA_data2 <- read.csv("lncRNA_A549.csv", header = T)
head(lncRNA_data2)
cor.test(lncRNA_data2$score, lncRNA_data2$TC, method = "spearman")
a1 <- lncRNA_data2 %>% filter(TC <= x)
count(a1)
a1$type <- "lncRNA (low TC)"
b1 <- lncRNA_data2 %>% filter(TC > x)
count(b1)
b1$type <- "lncRNA (high TC)"
mRNA_data2 <- read.csv("mRNA_A549.csv", header = T)
head(mRNA_data2)
cor.test(mRNA_data2$score, mRNA_data2$TC, method = "spearman")
c1 <- mRNA_data2 %>% filter(TC <= x)
count(c1)
c1$type <- "mRNA (low TC)"
d1 <- mRNA_data2 %>% filter(TC > x)
count(d1)
d1$type <- "mRNA (high TC)"
total_data1 <- data.frame(rbind(a1,b1,c1,d1))
median(lncRNA_data2$score)
median(mRNA_data2$score)
ks.test(lncRNA_data2$score, mRNA_data2$score, alternative = "two.sided")

tiff("/home/dell/Documents/koushiki/final_fig/SE_ridge2.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(total_data1, aes(x = score, y = type, fill = type)) +
  geom_density_ridges() + theme_ridges() + theme(legend.position = "none") + 
  scale_fill_manual(values = c("#3A5199","#31A2AC","#CB0000","#EB8A3E"))
p + xlab("Splicing Efficiency") + ylab("")+ xlim(c(-0.4,1.2)) +
  theme(plot.margin = margin(0.25, -0.5, 0.25, -2.75, "cm"),
        axis.title.x = element_text(size = 19, color = "black", face = "bold", hjust=0.6),
        axis.text = element_text(size = 16, color = "black", face="bold"),
        plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.6),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ggtitle("A549 Cell Line")
dev.off()
