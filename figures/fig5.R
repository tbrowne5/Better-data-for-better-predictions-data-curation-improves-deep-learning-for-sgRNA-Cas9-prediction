
# Requirements...
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggExtra)
library(cowplot)
library(scales)


# Fig 5A

citro_ribbon <- read.csv("data/Fig5A_TevSpCas9.csv",row.names=1)
eSpCas9_ribbon <- read.csv("data/Fig5A_eSpCas9.csv",row.names=1)
SpCas9_ribbon <- read.csv("data/Fig5A_WTSpCas9.csv",row.names=1)

custom_labels_ribbon <- c(-500, -400, -300, seq(-200,-20,by=20),"sgRNA",seq(20,200,by=20), 300, 400, 500)

SLO_plot <- ggplot(citro_ribbon[4:404,], aes(x = Cutoff, y = mean_spearman), fill = "purple") +
  theme_classic() +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "purple", alpha = 0.15) +
  geom_line(color = "purple", size = 1) +
  geom_ribbon(data=eSpCas9_ribbon[54:454,], aes(ymin = ymin, ymax = ymax), fill = "green", alpha = 0.25) +
  geom_line(data=eSpCas9_ribbon[54:454,], color = "darkgreen", size = 1) +
  geom_ribbon(data=SpCas9_ribbon[54:454,], aes(ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.1) +
  geom_line(data=SpCas9_ribbon[54:454,], color = "#1E90FF", size = 1) +
  ylab(expression(Delta ~ "Spearman correlation")) + 
  geom_hline(yintercept = c(0), alpha = 1) +
  scale_x_continuous(breaks = seq(-260, 260, by = 20), labels = custom_labels_ribbon, limits = c(-260, 260)) + 
  scale_y_continuous(breaks = seq(0.0, 0.12, by = 0.04), labels = seq(0.00, 0.12, by = 0.04), limits = c(-0.01, 0.14)) +
  xlab("Included nucleotides relative to the sgRNA target") +
  labs(color = element_blank(), tag = "A") +
  theme(plot.tag = element_text(face = 'bold')) +
  geom_point(data = citro_ribbon[c("1","2","3","405","406","407"),], aes(x= c(-260,-240,-220,220,240,260), y=mean_spearman), color = "purple", size=2) +
  geom_point(data = SpCas9_ribbon[c("1","2","3","505","506","507"),], aes(x= c(-260,-240,-220,220,240,260), y=mean_spearman), color = "#1E90FF", size=2) +
  geom_point(data = eSpCas9_ribbon[c("1","2","3","394","395","396"),], aes(x= c(-260,-240,-220,220,240,260), y=mean_spearman), color = "darkgreen", size=2) +
  geom_vline(xintercept = c(-250,-230,-210,210,230,250), linetype="dashed", color="grey")
print(SLO_plot)


# Fig 5B-D

for(i in 0:2){
  nuclease_colouring = TRUE
  dataNum <- i
  if(dataNum == 0){
    data <- read.csv("data/Fig5C_eSpCas9.tsv", header = TRUE, row.names=1,sep="\t") #pnsi_tevsacas9_SLO.txt
    pTag <- "C"
  }
  if(dataNum == 1){
    data <- read.csv("data/Fig5D_WTSpCas9.tsv", header = TRUE, row.names=1,sep="\t") #pnsi_tevsacas9_SLO.txt
    pTag <- "D"
  }
  if(dataNum == 2){
    data <- read.csv("data/Fig5B_TevSpCas9.tsv", header = TRUE, row.names=1,sep="\t") #pnsi_tevsacas9_SLO.txt
    pTag <- "B"
  }
  
  if(nuclease_colouring){
    colour_start = "#3f3f3f" #C49102"
    colour_mid = "#ECECEC"
    if(dataNum == 0){
      colour_end = "darkgreen"
    }
    if(dataNum == 1){
      colour_end = "blue"
      colour_start = "darkorange"
    }
    if(dataNum == 2){
      colour_end = "#6C3BAA"
      colour_start = "darkorange"
    }
  } else {
    colour_start = "blue"
    colour_mid = "#ECECEC"
    colour_end = "#ff4b33"
  }
  
  data <- data * -1
  data[data != 0] <- data[data != 0] - data[162,11]
  data[162,11] <- -10 # USED TO DISTINGUISH THE G FROM NGG
  data <- data[126:187,]
  data2 <- data
  
  data_long <- data %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "Value")
  
  # Add row numbers
  data_long <- data_long %>%
    mutate(Row = rep(c(-35:-1,1:27), each = 16))
  
  # Separate zero and non-zero values
  data_non_zero <- data_long %>% filter(Value != 0)
  data_zero <- data_long %>% filter(Value == 0)
  data_G <- data_long %>% filter(Value == -10)
  
  # Custom labels
  custom_labels <- c(sapply(-35:-2, function(x) paste(x,"\n", x+1, sep="")),"-1\nN","N\nG","G\nG","G\n1",sapply(1:24, function(x) paste(x,"\n", x+1, sep="")))
  
  # Create the heatmap
  heatmap_plot <- ggplot() +
    geom_tile(data = data_non_zero, aes(y = Feature, x = factor(Row), fill = Value), color = "darkgrey") +
    geom_tile(data = data_G, aes(y = Feature, x = factor(Row)), fill = "white", color = "darkgrey") +
    scale_fill_gradientn(colors = c(colour_start, colour_mid, colour_end), breaks = c(-1, 0, 1), limits=c(-1.5,1.5), oob=squish) +
    theme_minimal() +
    scale_x_discrete(labels = custom_labels) +
    labs(tag=pTag,x = "Dinucleotide position relative to the NGG PAM", y = "Dinucleotides", fill = "Mean\nscore\n") + #, title = title_text) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6), legend.key.height = unit(0.4, "cm"), plot.tag = element_text(face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank())
  
  if(dataNum == 0){
    heatmap_plot_eSp <- heatmap_plot
  }
  if(dataNum == 1){
    heatmap_plot_Sp <- heatmap_plot
  }
  if(dataNum == 2){
    heatmap_plot_Cit <- heatmap_plot
  }
}

grid.arrange(heatmap_plot_Cit, heatmap_plot_eSp, heatmap_plot_Sp, heights=c(1,1,1))


# Fig 4A-D
grid.arrange(SLO_plot, heatmap_plot_Cit, heatmap_plot_eSp, heatmap_plot_Sp, heights=c(1.2,1,1,1))
