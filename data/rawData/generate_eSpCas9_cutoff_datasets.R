
library("ALDEx2")
library("reshape2")

all_ecoli_sites_found <- bind_rows(read.csv("E_coli_sgRNA_target_sites_part_1.csv",header=FALSE),read.csv("E_coli_sgRNA_target_sites_part_2.csv",header=FALSE),read.csv("E_coli_sgRNA_target_sites_part_3.csv",header=FALSE))
all_ecoli_sites_found$sgRNAs <- all_ecoli_sites_found$V1
all_ecoli_sites_found$V1 <- NULL
all_ecoli_sites_found$sgRNA20 <- substring(all_ecoli_sites_found[,1],541,560)

Citro_tevspcas9 <- read.csv("E_coli_eSpCas9_counts.tsv",sep="\t",header=TRUE,row.names=1)

Citro_tevspcas9 <- Citro_tevspcas9[,c(3,4,1,2)]

outputLabel <- "Guo_eSpCas9_Curated_Min_dCas"
replicates <- 2
seedval <- 123
scaled <- FALSE

min_epi300_cutoff <- 1
max_epi300_cutoff <- 100
original_nrow <- nrow(Citro_tevspcas9)

metadata <- data.frame(cutoffVal = numeric(max_epi300_cutoff), sgRNAs = numeric(max_epi300_cutoff), proportionVsOne = numeric(max_epi300_cutoff), stringsAsFactors = FALSE)
cutoffOneVal <- 0

for (cutoff in min_epi300_cutoff:max_epi300_cutoff){
  print(original_nrow)
  Citro_tevspcas9 <- Citro_tevspcas9[which(((Citro_tevspcas9$eSpd1_NS + Citro_tevspcas9$eSpd2_NS)/2) >= cutoff),]
  print(paste(nrow(Citro_tevspcas9),":",cutoff,sep=""))
  conds <- c(rep("S",replicates),rep("NS",replicates))
  
  if(scaled){
    set.seed(seedval)
    Citro_tevspcas9_aldex.diff <- ALDEx2::aldex.clr(Citro_tevspcas9,conds, gamma=0.5)
  } else {
    set.seed(seedval)
    Citro_tevspcas9_aldex.diff <- ALDEx2::aldex.clr(Citro_tevspcas9,conds)
  }
  Citro_tevspcas9_aldex.diff <- ALDEx2::aldex.effect(Citro_tevspcas9_aldex.diff)
  Citro_tevspcas9_aldex.diff <- merge(all_ecoli_sites_found, Citro_tevspcas9_aldex.diff,by.x=2,by.y=0)
  Citro_tevspcas9_aldex.diff <- Citro_tevspcas9_aldex.diff[!duplicated(Citro_tevspcas9_aldex.diff$sgRNAs),]
  print(nrow(Citro_tevspcas9_aldex.diff))

  row.names(Citro_tevspcas9_aldex.diff) <- Citro_tevspcas9_aldex.diff$sgRNAs
  
  if(cutoff == min_epi300_cutoff){
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
  dir.create(paste(dataName,"_1118NT_",scaling,sep=''))
  write.csv(Citro_tevspcas9_aldex.diff[,c(2,6)],paste(dataName,"_1118NT_",scaling,"/",dataName,"_318NT_",scaling,".csv", sep=''),quote=FALSE, row.names=FALSE)
}
