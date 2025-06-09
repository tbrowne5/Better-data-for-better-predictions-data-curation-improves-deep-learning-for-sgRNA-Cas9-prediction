
library("ALDEx2")
library("reshape2")

all_ecoli_sites_found <- bind_rows(read.csv("E_coli_sgRNA_target_sites_part_1.csv",header=FALSE),read.csv("E_coli_sgRNA_target_sites_part_2.csv",header=FALSE),read.csv("E_coli_sgRNA_target_sites_part_3.csv",header=FALSE))
all_ecoli_sites_found$sgRNAs <- all_ecoli_sites_found$V1
all_ecoli_sites_found$V1 <- NULL
all_ecoli_sites_found$sgRNA20 <- substring(all_ecoli_sites_found[,1],541,560)

eSpCas9 <- read.csv("E_coli_eSpCas9_counts.tsv",sep="\t",header=TRUE,row.names=1)

eSpCas9 <- eSpCas9[,c(3,4,1,2)]

outputLabel <- "Guo_eSpCas9_Curated_Min_dCas"
replicates <- 2
seedval <- 123
scaled <- FALSE

min_dCas_cutoff <- 1
max_dCas_cutoff <- 2
original_nrow <- nrow(eSpCas9)

metadata <- data.frame(cutoffVal = numeric(max_dCas_cutoff), sgRNAs = numeric(max_dCas_cutoff), proportionVsOne = numeric(max_dCas_cutoff), stringsAsFactors = FALSE)
cutoffOneVal <- 0

for (cutoff in min_dCas_cutoff:max_dCas_cutoff){
  print(original_nrow)
  eSpCas9 <- eSpCas9[which(((eSpCas9$eSpd1_NS + eSpCas9$eSpd2_NS)/2) >= cutoff),]
  print(paste(nrow(eSpCas9),":",cutoff,sep=""))
  conds <- c(rep("S",replicates),rep("NS",replicates))
  
  if(scaled){
    set.seed(seedval)
    eSpCas9_aldex.diff <- ALDEx2::aldex.clr(eSpCas9,conds, gamma=0.5)
  } else {
    set.seed(seedval)
    eSpCas9_aldex.diff <- ALDEx2::aldex.clr(eSpCas9,conds)
  }
  eSpCas9_aldex.diff <- ALDEx2::aldex.effect(eSpCas9_aldex.diff)
  eSpCas9_aldex.diff <- merge(all_ecoli_sites_found, eSpCas9_aldex.diff,by.x=2,by.y=0)
  eSpCas9_aldex.diff <- eSpCas9_aldex.diff[!duplicated(eSpCas9_aldex.diff$sgRNAs),]
  print(nrow(eSpCas9_aldex.diff))

  row.names(eSpCas9_aldex.diff) <- eSpCas9_aldex.diff$sgRNAs
  
  if(cutoff == min_dCas_cutoff){
    cutoffOneVal <- nrow(eSpCas9_aldex.diff)
  }
  metadata[cutoff,] <- c(cutoff, nrow(eSpCas9_aldex.diff), (nrow(eSpCas9_aldex.diff)/cutoffOneVal))
  
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
  dir.create(paste(dataName,"_1118NT_",scaling,sep=''))
  write.csv(eSpCas9_aldex.diff[,c(2,6)],paste(dataName,"_1118NT_",scaling,"/",dataName,"_1118NT_",scaling,".csv", sep=''),quote=FALSE, row.names=FALSE)
}
