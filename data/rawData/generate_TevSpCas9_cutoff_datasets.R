
library("ALDEx2")
library("reshape2")

Citro_guides <- read.csv("C_rodentium_sgRNA_target_sites.csv",sep="\t",header=TRUE)

Citro_tevspcas9 <- read.csv("C_rodentium_TevSpCas9_counts.csv",sep=",",header=TRUE)
row.names(Citro_tevspcas9) <- Citro_tevspcas9$X
Citro_tevspcas9 <- Citro_tevspcas9[,c(11,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10)]

outputLabel <- "Citro_TevSpCas9_Curated_Min_Epi300"
replicates <- 9
seedval <- 123
scaled <- FALSE

max_epi300_cutoff <- 100

metadata <- data.frame(cutoffVal = numeric(max_epi300_cutoff), sgRNAs = numeric(max_epi300_cutoff), proportionVsOne = numeric(max_epi300_cutoff), stringsAsFactors = FALSE)
cutoffOneVal <- 0

for (cutoff in 1:max_epi300_cutoff){
  print(nrow(Citro_tevspcas9))
  Citro_tevspcas9 <- Citro_tevspcas9[which(Citro_tevspcas9$e2 >= cutoff),]
  print(nrow(Citro_tevspcas9))
  conds <- c(rep("S",replicates),rep("NS",replicates))
  
  if(scaled){
    set.seed(seedval)
    Citro_tevspcas9_aldex.diff <- ALDEx2::aldex.clr(Citro_tevspcas9,conds, gamma=0.5)
  } else {
    set.seed(seedval)
    Citro_tevspcas9_aldex.diff <- ALDEx2::aldex.clr(Citro_tevspcas9,conds)
  }
  print(nrow(Citro_tevspcas9_aldex.diff))
  Citro_tevspcas9_aldex.diff <- ALDEx2::aldex.effect(Citro_tevspcas9_aldex.diff)
  Citro_tevspcas9_aldex.diff <- merge(Citro_guides, Citro_tevspcas9_aldex.diff,by.x=1,by.y=0)
  Citro_tevspcas9_aldex.diff$sgRNAs <- substr(Citro_tevspcas9_aldex.diff$sgRNAs_long,21,48)
  
  row.names(Citro_tevspcas9_aldex.diff) <- Citro_tevspcas9_aldex.diff$sgRNAs
  
  if(cutoff == 1){
    cutoffOneVal <- nrow(Citro_tevspcas9_aldex.diff)
  }
  metadata[cutoff,] <- c(cutoff, nrow(Citro_tevspcas9_aldex.diff), (nrow(Citro_tevspcas9_aldex.diff)/cutoffOneVal))
  
  # Label
  if(cutoff == 0){
    dataName <- outputLabel
  } else {
    dataName <- paste(outputLabel,"_",cutoff,sep="")
  }
  if(scaled){
    scaling <- "Scaled"
  } else {
    scaling <- "Unscaled"
  }
  dir.create(paste(dataName,"_long_",scaling,sep=''))
  write.csv(Citro_tevspcas9_aldex.diff[,c(2,6)],paste(dataName,"_long_",scaling,"/",dataName,"_long_",scaling,".csv", sep=''),quote=FALSE, row.names=FALSE)
}