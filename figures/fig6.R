
# Required Packages and Density Function

library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(grid)
library(gridExtra)
library(ggExtra)
library(cowplot)
library(scales)

get_density <- function(x, y, n = 100, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Figure 6A

data <- read.csv(paste("data/Fig6A_pTox_TevSpCas9.txt",sep=""),sep="\t")

non0xmax <- 3
non0xmin <- -3.5
non0ymax <- 5
non0ymin <- -5

data$density <- get_density(data$prediction,data$actualscore)

dataplot <- NULL
dataplot <- ggplot(data[,c(2,1)])
dataplot <- dataplot + geom_point(aes(data$actualscore, data$prediction, color=abs(data$density)), size=1.5) + labs(tag = "A", color='Density')
dataplot <- dataplot + scale_color_viridis(option="magma",begin=0.1,end=0.9,direction=-1) + theme(panel.background = element_rect(fill = "#F0F0F8", color = NA), plot.tag = element_text(face = 'bold'), plot.margin = unit(c(5, 5, 5, 5), "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.5,'cm'), legend.title = element_text(size=7),legend.text = element_text(size=6)) + scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-4,-2, 0, 2, 4)) + scale_y_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-4,-2, 0, 2, 4)) + xlab("Activity Score") + ylab("Prediction")
dataplot = dataplot + coord_cartesian(ylim=c(non0ymin,non0ymax),xlim=c(non0xmin,non0xmax),clip="off")

print(dataplot)

cor(data[,2],data[,3],method="spearman")
cor(data[,2],data[,3],method="pearson")


# Figure 6B

data <- read.csv(paste("data/Fig6B_pTox_SpCas9.txt",sep=""),sep="\t")

non0xmax <- 5
non0xmin <- -5
non0ymax <- 5
non0ymin <- -5

data$density <- get_density(data$prediction,data$actualscore)

dataplot <- NULL
dataplot <- ggplot(data[,c(2,1)])
dataplot <- dataplot + geom_point(aes(data$actualscore, data$prediction, color=abs(data$density)), size=1.5) + labs(tag = "B", color='Density')
dataplot <- dataplot + scale_color_viridis(option="mako",begin=0,end=0.95,direction=-1) + theme(panel.background = element_rect(fill = "#F0F0F8", color = NA), plot.tag = element_text(face = 'bold'), plot.margin = unit(c(5, 5, 5, 5), "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.5,'cm'), legend.title = element_text(size=7),legend.text = element_text(size=6)) + scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-4,-2, 0, 2, 4)) + scale_y_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-4,-2, 0, 2, 4)) + xlab("Activity Score") + ylab("Prediction")
dataplot = dataplot + coord_cartesian(ylim=c(non0ymin,non0ymax),xlim=c(non0xmin,non0xmax),clip="off")

print(dataplot)

cor(data[,2],data[,3],method="spearman")
cor(data[,2],data[,3],method="pearson")


# Figure 6C

data <- read.csv(paste("data/Fig6C_KatG_TevSpCas9.txt",sep=""),sep="\t")

non0xmax <- 2.5
non0xmin <- -5.5
non0ymax <- 5
non0ymin <- -6

data$density <- get_density(data$prediction,data$actualscore)

dataplot <- NULL
dataplot <- ggplot(data[,c(2,1)])
dataplot <- dataplot + geom_point(aes(data$actualscore, data$prediction, color=abs(data$density)), size=1.5) + labs(tag = "C", color='Density')
dataplot <- dataplot + scale_color_viridis(option="magma",begin=0.1,end=0.9,direction=-1) + theme(panel.background = element_rect(fill = "#F0F0F8", color = NA), plot.tag = element_text(face = 'bold'), plot.margin = unit(c(5, 5, 5, 5), "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.5,'cm'), legend.title = element_text(size=7),legend.text = element_text(size=6)) + scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-4,-2, 0, 2, 4)) + scale_y_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-4,-2, 0, 2, 4)) + xlab("Activity Score") + ylab("Prediction")
dataplot = dataplot + coord_cartesian(ylim=c(non0ymin,non0ymax),xlim=c(non0xmin,non0xmax),clip="off")

print(dataplot)

cor(data[,2],data[,3],method="spearman")
cor(data[,2],data[,3],method="pearson")


# Fig 6D

pTox_Tev_performance <- read.csv("data/Fig6D_pTox_TevSpCas9.csv", row.names=1)
pTox_Tev_performance$model <- factor(pTox_Tev_performance$model,levels=pTox_Tev_performance$model)
ggplot(data=pTox_Tev_performance,aes(x=model,y=correlation,fill=colorVals)) + geom_col(fill=c("purple","#1E90FF","#737CA1","#737CA1")) + coord_cartesian(ylim = c(0.45,0.7)) + theme_classic() + theme(plot.tag = element_text(face = 'bold'), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + ylab("Spearman correlation") + labs(tag = "D")


#Fig 6E

pTox_Sp_performance <- read.csv("data/Fig6E_pTox_SpCas9.csv", row.names=1)
pTox_Sp_performance$model <- factor(pTox_Sp_performance$model,levels=pTox_Sp_performance$model)
ggplot(data=pTox_Sp_performance,aes(x=model,y=correlation,fill=colorVals)) + geom_col(fill=c("purple","#1E90FF","#737CA1","#737CA1")) + coord_cartesian(ylim = c(0.45,0.7)) + theme_classic() + theme(plot.tag = element_text(face = 'bold'), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + ylab("Spearman correlation") + labs(tag = "E")


# Fig 6F

KatG_Tev_performance <- read.csv("data/Fig6F_KatG_TevSpCas9.csv", row.names=1)
KatG_Tev_performance$model <- factor(KatG_Tev_performance$model,levels=KatG_Tev_performance$model)
ggplot(data=KatG_Tev_performance,aes(x=model,y=correlation,fill=colorVals)) + geom_col(fill=c("purple","gold","gold","#1E90FF","#737CA1","#737CA1")) + coord_cartesian(ylim = c(0.55,0.7)) + theme_classic() + theme(plot.tag = element_text(face = 'bold'), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + ylab("Spearman correlation") + labs(tag = "F")
