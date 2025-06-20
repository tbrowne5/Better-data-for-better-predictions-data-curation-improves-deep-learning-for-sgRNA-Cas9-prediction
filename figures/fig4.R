
# Required Packages

library(scales)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(grid)
library(gridExtra)

get_density <- function(x, y, n = 100, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Fig 4A-C:

citro_aldex_all <- read.csv("data/Fig4A_TevSpCas9.csv")
citro_aldex_all$density <- get_density(citro_aldex_all$e2,citro_aldex_all$diff.btw)
citro_aldex_all_250 <- citro_aldex_all[which(citro_aldex_all$e2 <= 100),]
citro_aldex_all_250$density <- get_density(citro_aldex_all_250$e2,citro_aldex_all_250$diff.btw)
citro_score_trunc <- ggplot(data=citro_aldex_all_250,aes(x=e2,y=diff.btw,color=density)) + geom_point(size=1) + labs(tag="A", color='Density') +
  scale_y_continuous(labels = label_number(accuracy = 1), limits=c(-7,5), breaks=seq(-6,3,3)) +
  theme_bw() + theme(panel.grid = element_blank(), plot.tag = element_text(face = 'bold'), legend.key.width = unit(0.3, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 6),plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm")) +
  scale_color_viridis(option="magma",direction = -1,begin=0,end=0.9,labels = label_number(accuracy = 0.001)) + xlab("Control condition read count") + ylab("Activity score") + xlim(c(1,100)) + geom_vline(xintercept = c(55.5), linetype="dashed")

espcas9_aldex_all = read.csv("data/Fig4B_eSpCas9.csv",row.names=1)
espcas9_aldex_all$e2 <- (espcas9_aldex_all$eSpd1_NS + espcas9_aldex_all$eSpd2_NS) /2 
citro_aldex_all_250 <- espcas9_aldex_all[which(espcas9_aldex_all$e2 <= 100),]
citro_aldex_all_250$density <- get_density(citro_aldex_all_250$e2,citro_aldex_all_250$diff.btw)
esp_score_trunc <- ggplot(data=citro_aldex_all_250,aes(x=e2,y=diff.btw,color=density)) + geom_point(size=1) + labs(tag="B", color='Density') +
  scale_y_continuous(labels = label_number(accuracy = 1), limits=c(-7,7), breaks=seq(-6,6,3)) +
  theme_bw() + theme(panel.grid = element_blank(), plot.tag = element_text(face = 'bold'), legend.key.width = unit(0.3, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 6),plot.margin = unit(c(0.25, 0, 0.25, 0), "cm")) +
  scale_color_viridis(direction = -1,begin=0.2,end=1,labels = label_number()) + xlab("Control condition read count") + ylab("Activity score") + xlim(c(1,100)) + geom_vline(xintercept = c(30.75), linetype="dashed")


spcas9_aldex_all = read.csv("data/Fig4C_WTSpCas9.csv",row.names=1)
spcas9_aldex_all$e2 <- (spcas9_aldex_all$dCas1_NS + spcas9_aldex_all$dCas2_NS) /2 
citro_aldex_all_250 <- spcas9_aldex_all[which(spcas9_aldex_all$e2 <= 100),]
citro_aldex_all_250$density <- get_density(citro_aldex_all_250$e2,citro_aldex_all_250$diff.btw)
sp_score_trunc <- ggplot(data=citro_aldex_all_250,aes(x=e2,y=diff.btw,color=density)) + geom_point(size=1) + labs(tag="C", color='Density') +
  scale_y_continuous(labels = label_number(accuracy = 1), limits=c(-11,5), breaks=seq(-10,5,5)) +
  theme_bw() + theme(panel.grid = element_blank(), plot.tag = element_text(face = 'bold'), legend.key.width = unit(0.3, "cm"),legend.title = element_text(size = 8), legend.text = element_text(size = 6),plot.margin = unit(c(0.25, 0.25, 0.25,0), "cm")) +
  scale_color_viridis(option="mako",direction = -1,begin=0.1,end=0.9,labels = label_number(accuracy = 0.001)) + xlab("Control condition read count") + ylab("Activity score") + xlim(c(1,100)) + geom_vline(xintercept = c(72.75), linetype="dashed")

grid.arrange(citro_score_trunc,esp_score_trunc,sp_score_trunc, widths=c(1,1,1))


# Fig 4D-E:

citro_stats <- read.csv("data/Fig4D_TevSpCas9.csv",sep=",",header=TRUE,row.names=1)
citro_stats$reduction <- citro_stats$cutoffVal - 1
for(i in 2:nrow(citro_stats)){
  citro_stats[i,"reduction"] <- citro_stats[i-1,"sgRNAs"] - citro_stats[i,"sgRNAs"]
}
citro_stats$reductionProportion <- citro_stats$reduction / max(citro_stats$sgRNAs)

spcas9_stats <- read.csv("data/Fig4D_WTSpCas9.csv",sep=",",header=TRUE,row.names=1)
spcas9_stats$reduction <- spcas9_stats$cutoffVal - 1
for(i in 2:nrow(spcas9_stats)){
  spcas9_stats[i,"reduction"] <- spcas9_stats[i-1,"sgRNAs"] - spcas9_stats[i,"sgRNAs"]
}
spcas9_stats$reductionProportion <- spcas9_stats$reduction / max(spcas9_stats$sgRNAs)

espcas9_stats <- read.csv("data/Fig4D_eSpCas9.csv",sep=",",header=TRUE,row.names=1)
espcas9_stats$reduction <- espcas9_stats$cutoffVal - 1
for(i in 2:nrow(espcas9_stats)){
  espcas9_stats[i,"reduction"] <- espcas9_stats[i-1,"sgRNAs"] - espcas9_stats[i,"sgRNAs"]
}
espcas9_stats$reductionProportion <- espcas9_stats$reduction / max(espcas9_stats$sgRNAs)

citro_stats$proportionVsOne <- citro_stats$sgRNAs
spcas9_stats$proportionVsOne <- spcas9_stats$sgRNAs
espcas9_stats$proportionVsOne <- espcas9_stats$sgRNAs

citro_cutoff_proportion_plot = ggplot() +
  geom_line(data=spcas9_stats,aes(x=cutoffVal,y=proportionVsOne),color="#1E90FF",size=1) + 
  geom_line(data=citro_stats,aes(x=cutoffVal,y=proportionVsOne),color="purple",size=1) +
  geom_line(data=espcas9_stats,aes(x=cutoffVal,y=proportionVsOne),color="#02BA0F",size=1) +
  xlim(c(1,100)) + theme_bw() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 7), panel.grid = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank(), plot.tag = element_text(face = 'bold'),plot.margin = unit(c(0.5, 0.5, 0, 0.25), "cm")) + ylab("Data Points") + labs(tag="D")

delta_cit_top5_df <- read.csv("data/Fig4E_TevSpCas9.csv",row.names=1)
delta_eSp_top5_df <- read.csv("data/Fig4E_eSpCas9.csv",row.names=1)
delta_Sp_top5_df <- read.csv("data/Fig4E_WTSpCas9.csv",row.names=1)

citro_cutoff_plot = ggplot() +
  geom_vline(xintercept = c(20), size=0.5, linetype="dashed", color="darkgrey") +
  geom_smooth(data=delta_cit_top5_df, aes(x=Cutoff, y=mean_spearman), method = "loess", color = "purple", se = TRUE, span = 0.5) +
  geom_smooth(data=delta_eSp_top5_df, aes(x=Cutoff, y=mean_spearman), method = "loess", color = "#02BA0F", se = TRUE, span = 0.5) +
  geom_smooth(data=delta_Sp_top5_df, aes(x=Cutoff, y=mean_spearman), method = "loess", color = "#1E90FF", se = TRUE, span = 0.5) +
  geom_point(data = subset(delta_cit_top5_df, Cutoff == 56), aes(x = Cutoff, y = mean_spearman), 
             color = "purple", size = 3, shape = 21, stroke = 1.2, fill = "white") +
  geom_point(data = subset(delta_eSp_top5_df, Cutoff == 31), aes(x = Cutoff, y = mean_spearman), 
             color = "#02BA0F", size = 3, shape = 21, stroke = 1.2, fill = "white") +
  geom_point(data = subset(delta_Sp_top5_df, Cutoff == 73), aes(x = Cutoff, y = mean_spearman), 
             color = "blue", size = 3, shape = 21, stroke = 1.2, fill = "white") +
  geom_point(data=delta_cit_top5_df, aes(x=Cutoff, y=mean_spearman), alpha = 0.3, color="purple") + 
  geom_point(data=delta_eSp_top5_df, aes(x=Cutoff, y=mean_spearman), alpha = 0.4, color="#02BA0F") + 
  geom_point(data=delta_Sp_top5_df, aes(x=Cutoff, y=mean_spearman), alpha = 0.3, color="blue") + 
  xlab("Control condition minimum read count cutoff") + ylab(expression(Delta ~ "Spearman correlation")) + theme_bw() +
  theme(panel.grid = element_blank(), plot.tag = element_text(face = 'bold'),plot.margin = unit(c(-0.2, 0.5, 0.5, 0.25), "cm")) + xlim(c(1,100)) + labs(tag="E")

grid.arrange(citro_cutoff_proportion_plot,citro_cutoff_plot, heights=c(1,3))


# Plot Fig 4

cutoff_layout <- rbind(
  c(1, 2, 3),
  c(4, 4, 4),
  c(5, 5, 5)
)

grid.arrange(citro_score_trunc, esp_score_trunc, sp_score_trunc, citro_cutoff_proportion_plot,citro_cutoff_plot, layout_matrix = cutoff_layout, heights=c(1,0.6,1))
grid.arrange(citro_score_trunc, esp_score_trunc, sp_score_trunc, citro_cutoff_proportion_plot,citro_cutoff_plot, layout_matrix = cutoff_layout, heights=c(1,0.6,1.4))
