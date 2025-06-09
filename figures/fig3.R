
# Required Packages and Density Function

library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)

get_density <- function(x, y, n = 100, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Figure 3A

data <- read.csv(paste("data/Fig3A_TevSpCas9.txt",sep=""),sep="\t")
non0xmax <- max(data$actualscore)
non0xmin <- min(data$actualscore)
non0ymax <- max(data$prediction)
non0ymin <- min(data$prediction)

data$density <- get_density(data$prediction,data$actualscore)

dataplot <- NULL
dataplot <- ggplot(data[,c(2,1)])
dataplot <- dataplot + geom_point(aes(data$actualscore, data$prediction, color=abs(data$density)), size=0.5) + labs(tag = "A", color='Density')
dataplot <- dataplot + scale_color_viridis(option="magma",begin=0.1,end=0.9,direction=-1) + theme(panel.background = element_rect(fill = "white", color = NA), plot.tag = element_text(face = 'bold'), plot.margin = unit(c(5, 5, 5, 5), "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.5,'cm'), legend.title = element_text(size=7),legend.text = element_text(size=6)) + scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-8,-4, 0, 4)) + scale_y_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-8,-4, 0, 4)) + xlab("Activity Score") + ylab("Prediction")
dataplot = dataplot + coord_cartesian(ylim=c(non0ymin,non0ymax),xlim=c(non0xmin,non0xmax),clip="on")

print(dataplot)

cor(data[,2],data[,3],method="spearman")
cor(data[,2],data[,3],method="pearson")


# Figure 3B

data <- read.csv(paste("data/Fig3B_eSpCas9.txt",sep=""),sep="\t")
non0xmax <- max(data$actualscore)
non0xmin <- min(data$actualscore) +1
non0ymax <- max(data$prediction)
non0ymin <- min(data$prediction) +0.25

data$density <- get_density(data$prediction,data$actualscore)

dataplot <- NULL
dataplot <- ggplot(data[,c(2,1)])
dataplot <- dataplot + geom_point(aes(data$actualscore, data$prediction, color=abs(data$density)), size=0.5) + labs(tag = "B", color='Density')
dataplot <- dataplot + scale_color_viridis(begin=0,end=0.95,direction=-1) + theme(panel.background = element_rect(fill = "white", color = NA), plot.tag = element_text(face = 'bold'), plot.margin = unit(c(5, 5, 5, 5), "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.5,'cm'), legend.title = element_text(size=7),legend.text = element_text(size=6)) + scale_y_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-8,-4, 0, 4)) + xlab("Activity Score") + ylab("Prediction")
dataplot <- dataplot + coord_cartesian(ylim=c(non0ymin,non0ymax),xlim=c(non0xmin,non0xmax),clip="on")

print(dataplot)

cor(data[,2],data[,3],method="spearman")
cor(data[,2],data[,3],method="pearson")


# Fig 3C

data <- read.csv(paste("data/Fig3C_WTSpCas9.txt",sep=""),sep="\t")
non0xmax <- max(data$actualscore)
non0xmin <- min(data$actualscore)
non0ymax <- max(data$prediction)
non0ymin <- min(data$prediction)

data$density <- get_density(data$prediction,data$actualscore)

dataplot <- NULL
dataplot <- ggplot(data[,c(2,1)])
dataplot <- dataplot + geom_point(aes(data$actualscore, data$prediction, color=abs(data$density)), size=0.5) + labs(tag = "C", color='Density')
dataplot <- dataplot + scale_color_viridis(option="mako",begin=0,end=0.95,direction=-1) + theme(panel.background = element_rect(fill = "white", color = NA), plot.tag = element_text(face = 'bold'), plot.margin = unit(c(5, 5, 5, 5), "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.5,'cm'), legend.title = element_text(size=7),legend.text = element_text(size=6)) + scale_y_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-8,-4, 0, 4)) + scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = c(-8,-4, 0, 4)) + xlab("Activity Score") + ylab("Prediction")
dataplot <- dataplot + coord_cartesian(ylim=c(non0ymin,non0ymax),xlim=c(non0xmin,non0xmax),clip="on")

print(dataplot)

cor(data[,2],data[,3],method="spearman")
cor(data[,2],data[,3],method="pearson")

