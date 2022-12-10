#number of transcripts involved in different Alternative Splicing Events

data1 <- read.table("as_m_A3_strict.ioe", header=T)
data2 <- read.table("as_m_A5_strict.ioe", header=T)
data3 <- read.table("as_m_AF_strict.ioe", header=T)
data4 <- read.table("as_m_AL_strict.ioe", header=T)
data5 <- read.table("as_m_MX_strict.ioe", header=T)
data6 <- read.table("as_m_RI_strict.ioe", header=T)
data7 <- read.table("as_m_SE_strict.ioe", header=T)

data_lncRNA <- read.csv("lncrna_main_m.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
data_lncRNA$A3 <- ""
data_lncRNA$A5 <- ""
data_lncRNA$AF <- ""
data_lncRNA$AL <- ""
data_lncRNA$MX <- ""
data_lncRNA$RI <- ""
data_lncRNA$SE <- ""

for(i in 1:nrow(data_lncRNA)){
  test1 = grepl(data_lncRNA$Gene[i], data1$gene_id[1:nrow(data1)])
  if(sum(test1)>=1){
    data_lncRNA$A3[i] = sum(test1)
  }
  else{
    data_lncRNA$A3[i] = 0
  }
}
for(i in 1:nrow(data_lncRNA)){
  test2 = grepl(data_lncRNA$Gene[i], data2$gene_id[1:nrow(data2)])
  if(sum(test2)>=1){
    data_lncRNA$A5[i] = sum(test2)
  }
  else{
    data_lncRNA$A5[i] = 0
  }
}
for(i in 1:nrow(data_lncRNA)){
  test3 = grepl(data_lncRNA$Gene[i], data3$gene_id[1:nrow(data3)])
  if(sum(test3)>=1){
    data_lncRNA$AF[i] = sum(test3)
  }
  else{
    data_lncRNA$AF[i] = 0
  }
}
for(i in 1:nrow(data_lncRNA)){
  test4 = grepl(data_lncRNA$Gene[i], data4$gene_id[1:nrow(data4)])
  if(sum(test4)>=1){
    data_lncRNA$AL[i] = sum(test4)
  }
  else{
    data_lncRNA$AL[i] = 0
  }
}
for(i in 1:nrow(data_lncRNA)){
  test5 = grepl(data_lncRNA$Gene[i], data5$gene_id[1:nrow(data5)])
  if(sum(test5)>=1){
    data_lncRNA$MX[i] = sum(test5)
  }
  else{
    data_lncRNA$MX[i] = 0
  }
}
for(i in 1:nrow(data_lncRNA)){
  test6 = grepl(data_lncRNA$Gene[i], data6$gene_id[1:nrow(data6)])
  if(sum(test6)>=1){
    data_lncRNA$RI[i] = sum(test6)
  }
  else{
    data_lncRNA$RI[i] = 0
  }
}
for(i in 1:nrow(data_lncRNA)){
  test7 = grepl(data_lncRNA$Gene[i], data7$gene_id[1:nrow(data7)])
  if(sum(test7)>=1){
    data_lncRNA$SE[i] = sum(test7)
  }
  else{
    data_lncRNA$SE[i] = 0
  }
}
write.csv(data_lncRNA, file = "lncrna_ASE_m.csv", col.names = TRUE, row.names = F)

data_CDS <- read.csv("mrna_main_m.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS$A3 <- ""
data_CDS$A5 <- ""
data_CDS$AF <- ""
data_CDS$AL <- ""
data_CDS$MX <- ""
data_CDS$RI <- ""
data_CDS$SE <- ""

for(i in 1:nrow(data_CDS)){
  test1 = grepl(data_CDS$Gene[i], data1$gene_id[1:nrow(data1)])
  if(sum(test1)>=1){
    data_CDS$A3[i] = sum(test1)
  }
  else{
    data_CDS$A3[i] = 0
  }
}
for(i in 1:nrow(data_CDS)){
  test2 = grepl(data_CDS$Gene[i], data2$gene_id[1:nrow(data2)])
  if(sum(test2)>=1){
    data_CDS$A5[i] = sum(test2)
  }
  else{
    data_CDS$A5[i] = 0
  }
}
for(i in 1:nrow(data_CDS)){
  test3 = grepl(data_CDS$Gene[i], data3$gene_id[1:nrow(data3)])
  if(sum(test3)>=1){
    data_CDS$AF[i] = sum(test3)
  }
  else{
    data_CDS$AF[i] = 0
  }
}
for(i in 1:nrow(data_CDS)){
  test4 = grepl(data_CDS$Gene[i], data4$gene_id[1:nrow(data4)])
  if(sum(test4)>=1){
    data_CDS$AL[i] = sum(test4)
  }
  else{
    data_CDS$AL[i] = 0
  }
}
for(i in 1:nrow(data_CDS)){
  test5 = grepl(data_CDS$Gene[i], data5$gene_id[1:nrow(data5)])
  if(sum(test5)>=1){
    data_CDS$MX[i] = sum(test5)
  }
  else{
    data_CDS$MX[i] = 0
  }
}
for(i in 1:nrow(data_CDS)){
  test6 = grepl(data_CDS$Gene[i], data6$gene_id[1:nrow(data6)])
  if(sum(test6)>=1){
    data_CDS$RI[i] = sum(test6)
  }
  else{
    data_CDS$RI[i] = 0
  }
}
for(i in 1:nrow(data_CDS)){
  test7 = grepl(data_CDS$Gene[i], data7$gene_id[1:nrow(data7)])
  if(sum(test7)>=1){
    data_CDS$SE[i] = sum(test7)
  }
  else{
    data_CDS$SE[i] = 0
  }
}
write.csv(data_CDS, file = "mrna_ASE_m.csv", col.names = TRUE, row.names = F)

