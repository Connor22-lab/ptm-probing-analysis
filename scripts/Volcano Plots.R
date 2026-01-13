#R script for plotting volcano plots using ggplot2

#using Candiates lists exported from spectronaut

rm(list=ls())#clears global environment

#set working directory to your project file (includes all relevent data tables etc)
setwd("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw")


#load relevent packages

library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(tidyverse)
library(data.table)
library(plyr)
library(ggpubr)

#set theme- for aesthetics- can customise this
theme_set(theme_bw(base_size = 15) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              plot.title = element_text(hjust = 0.5)
            ))

#load in data- Cands refers to table of results from spectronaut differntial abundance testing (t test)

#

Cands1<-as.data.frame(fread("PTM1.tsv"))

#subset out specific comparisons required
Cands1_LPS_act <- subset(Cands1, Cands1$`Comparison (group1/group2)` == "LPS / LPS+Nigericin")
Cands1_P3C4_act <- subset(Cands1, Cands1$`Comparison (group1/group2)` == "P3C4 / P3C4+Nigericin")

# Flip log2FC so x-axis represents (LPS+Nigericin / LPS) and (P3C4+Nigericin / P3C4)
Cands1_LPS_act$log2FC_state <- -Cands1_LPS_act$`AVG Log2 Ratio`
Cands1_P3C4_act$log2FC_state <- -Cands1_P3C4_act$`AVG Log2 Ratio`

#annotate up or down regulated
#Log2FC 0.58 = FC 1.5, Q<0.05 - spectronaut default
Cands1_LPS_act$Enrichment <- "NO"
Cands1_LPS_act$Enrichment[Cands1_LPS_act$log2FC_state >= 0.58 & Cands1_LPS_act$Qvalue < 0.05] <- "UP"
Cands1_LPS_act$Enrichment[Cands1_LPS_act$log2FC_state <= -0.58 & Cands1_LPS_act$Qvalue < 0.05] <- "DOWN"


Cands1_P3C4_act$Enrichment <- "NO"
Cands1_P3C4_act$Enrichment[Cands1_P3C4_act$log2FC_state >= 0.58 & Cands1_P3C4_act$Qvalue < 0.05] <- "UP"
Cands1_P3C4_act$Enrichment[Cands1_P3C4_act$log2FC_state <= -0.58 & Cands1_P3C4_act$Qvalue < 0.05] <- "DOWN"
#save for reference
write.csv(Cands1_LPS_act, "Candidates_LPS_act_tr1.csv")
write.csv(Cands1_P3C4_act , "Candidates_P3C4_act_tr1.csv")


#want to change label colours to refelct if enriched or not


### all proteins- unlabelled
# 6h
p1<-ggplot(data = Cands1_LPS_act, aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Enrichment, alpha = Enrichment, size = Enrichment)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = c("DOWN" = "darkcyan", "UP" = "deeppink", "NO" = "grey"),
                     labels = c("DOWN","NO", "UP")) +#assign colour and label to points
  scale_alpha_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 0.5)) +#adjust transparency of points 0= transparent 1= solid
  scale_size_manual(values = c("DOWN" = 2, "UP" = 2, "NO" = 2)) +#size of points
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +# can set y and x axis limits if desired
  ggtitle("Activation: LPS+Nigericin vs LPS") +#main title
  guides(colour = guide_legend(override.aes = list(size = 6))) #fix legend

print(p1)


#save plot to working directory
#you can adjust the sizing as required
ggsave( "VolcanoPlot1.png", plot = p1, width = 32, height = 20, units = c("cm"))

# 24h
p2<-ggplot(data = Cands1_P3C4_act , aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Enrichment, alpha = Enrichment, size = Enrichment)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = c("DOWN" = "darkcyan", "UP" = "deeppink", "NO" = "grey"),
                     labels = c("DOWN","NO", "UP")) +#assign colour and label to points
  scale_alpha_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 0.5)) +#adjust transparency of points 0= transparent 1= solid
  scale_size_manual(values = c("DOWN" = 2, "UP" = 2, "NO" = 2)) +#size of points
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +# can set y and x axis limits if desired
  ggtitle("Activation: P3C4+Nigericin vs P3C4") +#main title
  guides(colour = guide_legend(override.aes = list(size = 6))) #fix legend
print(p2)

ggsave( "VolcanoPlot2.png", plot = p2, width = 32, height = 20, units = c("cm"))

#############################################################################
#label all changed proteins

#make labels - differentially expressed proteins
lab_LPS_act_up <- Cands1_LPS_act %>% filter(Cands1_LPS_act$Enrichment %in% "UP")
lab_LPS_act_down <- Cands1_LPS_act %>% filter(Cands1_LPS_act$Enrichment %in% "DOWN")
#make plots

p3<-ggplot(data = Cands1_LPS_act, aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Enrichment, alpha = Enrichment, size = Enrichment)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = c("DOWN" = "darkcyan", "UP" = "deeppink", "NO" = "grey"),
                     labels = c("DOWN","NO", "UP")) +#assign colour and label to points
  scale_alpha_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 0.5)) +#adjust transparency of points 0= transparent 1= solid
  scale_size_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 1)) +#size of points
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: LPS+Nigericin vs LPS") +#main title
  geom_text_repel(data = lab_LPS_act_up, aes(label = Genes), force = 2, col = "deeppink3", size = 3, max.overlaps = 15, fontface= "bold")+ #labelling
  geom_text_repel(data = lab_LPS_act_down, aes(label = Genes), force = 2, col = "darkcyan", size = 3, max.overlaps = 15, fontface= "bold")+
  guides(colour = guide_legend(override.aes = list(size = 6))) #fix legend
print(p3)

ggsave( "VolcanoPlot3.png", plot = p3, width = 32, height = 20, units = c("cm"))


lab_P3C4_act_up <- Cands1_P3C4_act  %>% filter(Cands1_P3C4_act$Enrichment %in% "UP")
lab_P3C4_act_down <- Cands1_P3C4_act  %>% filter(Cands1_P3C4_act$Enrichment %in% "DOWN")

p4<-ggplot(data = Cands1_P3C4_act , aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Enrichment, alpha = Enrichment, size = Enrichment)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = c("DOWN" = "darkcyan", "UP" = "deeppink", "NO" = "grey"),
                     labels = c("DOWN","NO", "UP")) +#assign colour and label to points
  scale_alpha_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 0.5)) +#adjust transparency of points 0= transparent 1= solid
  scale_size_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 1)) +#size of points
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: P3C4+Nigericin vs P3C4") +#main title
  geom_text_repel(data = lab_P3C4_act_up, aes(label = Genes), force = 2, col = "deeppink3", size = 3, max.overlaps = 15, fontface= "bold")+
  geom_text_repel(data = lab_P3C4_act_down, aes(label = Genes), force = 2, col = "darkcyan", size = 3, max.overlaps = 15, fontface= "bold")+
  guides(colour = guide_legend(override.aes = list(size = 6))) #fix legend

print(p4)

ggsave( "VolcanoPlot4.png", plot = p4, width = 32, height = 20, units = c("cm"))


###############################################################################
#want to label a subset of proteo=ins of interest

#load in some csv format file with uniprot ids matched to gene names
#can download these from uniprot- this example shows inflammasome related prots
#establish proteins of interest (poi)
poi <- as.data.frame(fread("POI.csv"))
poi_list<-poi$Entry
poi$ProteinGroups<-poi$Entry

poi_pattern <- paste(poi_list, collapse="|")

Cands1_LPS_act <- Cands1_LPS_act %>%
  mutate(Label = ifelse(grepl(poi_pattern, UniProtIds), "poi", "Other"))

# Create a new variable that combines Label and Enrichment
Cands1_LPS_act$Combined <- paste(Cands1_LPS_act$Label, Cands1_LPS_act$Enrichment, sep = "_")

#make labels - inflammasome poi
lab_LPS_act <- Cands1_LPS_act %>% filter(Cands1_LPS_act$Combined %in% "poi_NO")
lab_LPS_act_up <- Cands1_LPS_act %>% filter(Cands1_LPS_act$Combined %in% "poi_UP")
lab_LPS_act_down <- Cands1_LPS_act %>% filter(Cands1_LPS_act$Combined %in% "poi_DOWN")
#make plots

p5<-ggplot(data = Cands1_LPS_act, aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Enrichment, alpha = Enrichment, size = Enrichment)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = c("DOWN" = "darkcyan", "UP" = "deeppink", "NO" = "grey"),
                     labels = c("DOWN","NO", "UP")) +#assign colour and label to points
  scale_alpha_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 0.5)) +#adjust transparency of points 0= transparent 1= solid
  scale_size_manual(values = c("DOWN" = 2, "UP" = 2, "NO" = 2)) +#size of points
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: LPS+Nigericin vs LPS") +#main title
  geom_text_repel(data = lab_LPS_act, aes(label = Genes), force = 5, col = "grey43", size = 4, max.overlaps = 10, fontface= "bold")+ #label poi
  geom_text_repel(data = lab_LPS_act_up, aes(label = Genes), force = 5, col = "deeppink3", size = 4, max.overlaps = 10, fontface= "bold")+
  geom_text_repel(data = lab_LPS_act_down, aes(label = Genes), force = 5, col = "darkcyan", size = 4, max.overlaps = 10, fontface= "bold")+
  guides(colour = guide_legend(override.aes = list(size = 6))) #fix legend
print(p5)

ggsave( "VolcanoPlot5.png", plot = p5, width = 32, height = 20, units = c("cm"))


################################################################################

#alternative volcano plot- slightly fancier colour labelling

# Define the color, alpha, and size mappings for the combined variable
color_mapping <- c("poi_NO" = "black", "poi_UP" = "deeppink3","poi_DOWN" = "darkcyan",
                   "Other_NO" = "grey", "Other_UP" = "deeppink3","Other_DOWN" = "darkcyan")
alpha_mapping <- c("poi_NO" = 1, "poi_UP" = 1,"poi_DOWN" = 1,
                   "Other_NO" = 0.3, "Other_UP" = 0.3, "Other_DOWN" = 0.3)
size_mapping <- c("poi_NO" = 2, "poi_UP" = 2,"poi_DOWN" = 2,
                  "Other_NO" = 2, "Other_UP" = 2, "Other_DOWN" = 2)

p6<-ggplot(data = Cands1_LPS_act, aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Combined, alpha = Combined, size = Combined)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = color_mapping) +
  scale_alpha_manual(values = alpha_mapping) +
  scale_size_manual(values = size_mapping) + labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: LPS+Nigericin vs LPS") +#main title
  geom_text_repel(data = lab_LPS_act, aes(label = Genes), force = 5, col = "black", size = 4, max.overlaps = 10, fontface= "bold")+ #label poi
  geom_text_repel(data = lab_LPS_act_up, aes(label = Genes), force = 5, col = "deeppink", size = 4, max.overlaps = 10, fontface= "bold")+
  geom_text_repel(data = lab_LPS_act_down, aes(label = Genes), force = 5, col = "darkcyan", size = 4, max.overlaps = 10, fontface= "bold")+
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),  # Background of the entire plot area
        legend.background = element_rect(fill = "transparent", color = NA), # Background of the legend
        legend.box.background = element_rect(fill = "transparent", color = NA))

print(p6)

ggsave( "VolcanoPlot6.png", plot = p6, width = 32, height = 20, units = c("cm"))

#remove legend if desired

p7<-ggplot(data = Cands1_LPS_act, aes(x = log2FC_state, y = -log10(Qvalue))) + # assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') + # log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') + # log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') + # 0 intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + # qvalue intercept
  geom_point(aes(color = Combined, alpha = Combined, size = Combined)) + # tell ggplot what labels to use for plotting
  scale_color_manual(values = color_mapping, guide = "none") + # remove color legend
  scale_alpha_manual(values = alpha_mapping, guide = "none") + # remove alpha legend
  scale_size_manual(values = size_mapping, guide = "none") + # remove size legend
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) + # axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: LPS+Nigericin vs LPS") + # main title
  geom_text_repel(data = lab_LPS_act, aes(label = Genes), force = 5, col = "black", size = 4, max.overlaps = 10, fontface = "bold") + # label poi
  geom_text_repel(data = lab_LPS_act_up, aes(label = Genes), force = 5, col = "deeppink", size = 4, max.overlaps = 10, fontface = "bold") +
  geom_text_repel(data = lab_LPS_act_down, aes(label = Genes), force = 5, col = "darkcyan", size = 4, max.overlaps = 10, fontface = "bold") +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA), # Background of the entire plot area
        legend.background = element_rect(fill = "transparent", color = NA), # Background of the legend
        legend.box.background = element_rect(fill = "transparent", color = NA))

print(p7)

ggsave( "VolcanoPlot7.png", plot = p7, width = 32, height = 20, units = c("cm"))


Cands1_P3C4_act <- Cands1_P3C4_act %>%
  mutate(Label = ifelse(grepl(poi_pattern, UniProtIds), "poi", "Other"))

# Create a new variable that combines Label and Enrichment
Cands1_P3C4_act$Combined <- paste(Cands1_P3C4_act$Label, Cands1_P3C4_act$Enrichment, sep = "_")

#make labels - inflammasome poi
lab_P3C4_act <- Cands1_P3C4_act  %>% filter(Cands1_P3C4_act$Combined %in% "poi_NO")
lab_P3C4_act_up <- Cands1_P3C4_act  %>% filter(Cands1_P3C4_act$Combined %in% "poi_UP")
lab_P3C4_act_down <- Cands1_P3C4_act  %>% filter(Cands1_P3C4_act$Combined %in% "poi_DOWN")
#make plots

p8<-ggplot(data = Cands1_P3C4_act , aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Enrichment, alpha = Enrichment, size = Enrichment)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = c("DOWN" = "darkcyan", "UP" = "deeppink", "NO" = "grey"),
                     labels = c("DOWN","NO", "UP")) +#assign colour and label to points
  scale_alpha_manual(values = c("DOWN" = 1, "UP" = 1, "NO" = 0.5)) +#adjust transparency of points 0= transparent 1= solid
  scale_size_manual(values = c("DOWN" = 2, "UP" = 2, "NO" = 2)) +#size of points
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: P3C4+Nigericin vs P3C4") +#main title
  geom_text_repel(data = lab_P3C4_act, aes(label = Genes), force = 5, col = "grey43", size = 4, max.overlaps = 10, fontface= "bold")+ #label poi
  geom_text_repel(data = lab_P3C4_act_up, aes(label = Genes), force = 5, col = "deeppink3", size = 4, max.overlaps = 10, fontface= "bold")+
  geom_text_repel(data = lab_P3C4_act_down, aes(label = Genes), force = 5, col = "darkcyan", size = 4, max.overlaps = 10, fontface= "bold")+
  guides(colour = guide_legend(override.aes = list(size = 6))) #fix legend

print(p8)


ggsave( "VolcanoPlot8.png", plot = p8, width = 32, height = 20, units = c("cm"))


################################################################################


#alternativ volcano plot- slightly fancier colour labelling


# Define the color, alpha, and size mappings for the combined variable
color_mapping <- c("poi_NO" = "black", "poi_UP" = "deeppink3","poi_DOWN" = "darkcyan",
                   "Other_NO" = "grey", "Other_UP" = "deeppink3","Other_DOWN" = "darkcyan")
alpha_mapping <- c("poi_NO" = 1, "poi_UP" = 1,"poi_DOWN" = 1,
                   "Other_NO" = 0.3, "Other_UP" = 0.3, "Other_DOWN" = 0.3)
size_mapping <- c("poi_NO" = 2, "poi_UP" = 2,"poi_DOWN" = 2,
                  "Other_NO" = 2, "Other_UP" = 2, "Other_DOWN" = 2)

p9<-ggplot(data = Cands1_P3C4_act , aes(x = log2FC_state, y = -log10(Qvalue))) +#assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') +#log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') +# 0  intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + #qvalue intercept
  geom_point(aes(color = Combined, alpha = Combined, size = Combined)) + #tell ggplot what labels to use for plotting
  scale_color_manual(values = color_mapping) +
  scale_alpha_manual(values = alpha_mapping) +
  scale_size_manual(values = size_mapping) + labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) +#axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: P3C4+Nigericin vs P3C4") +#main title
  geom_text_repel(data = lab_P3C4_act, aes(label = Genes), force = 5, col = "black", size = 4, max.overlaps = 10, fontface= "bold")+ #label poi
  geom_text_repel(data = lab_P3C4_act_up, aes(label = Genes), force = 5, col = "deeppink", size = 4, max.overlaps = 10, fontface= "bold")+
  geom_text_repel(data = lab_P3C4_act_down, aes(label = Genes), force = 5, col = "darkcyan", size = 4, max.overlaps = 10, fontface= "bold")+
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),  # Background of the entire plot area
        legend.background = element_rect(fill = "transparent", color = NA), # Background of the legend
        legend.box.background = element_rect(fill = "transparent", color = NA))
print(p9)

ggsave( "VolcanoPlot9.png", plot = p9, width = 32, height = 20, units = c("cm"))


#remove legend if desired

p10<-ggplot(data = Cands1_P3C4_act , aes(x = log2FC_state, y = -log10(Qvalue))) + # assign data set and parameters
  geom_vline(xintercept = 0.58, col = "black", linetype = 'dashed') + # log2 FC intercept
  geom_vline(xintercept = -0.58, col = "black", linetype = 'dashed') + # log2 FC intercept
  geom_vline(xintercept = 0, col = "darkgrey", linetype = 'dashed') + # 0 intercept
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + # qvalue intercept
  geom_point(aes(color = Combined, alpha = Combined, size = Combined)) + # tell ggplot what labels to use for plotting
  scale_color_manual(values = color_mapping, guide = "none") + # remove color legend
  scale_alpha_manual(values = alpha_mapping, guide = "none") + # remove alpha legend
  scale_size_manual(values = size_mapping, guide = "none") + # remove size legend
  labs(x = expression("log"[2] * "FC (Active / Primed)"), y = expression("-log"[10] * "Q-value")) + # axis labels
  xlim(-3,6) +
  ylim(0,8) +
  ggtitle("Activation: P3C4+Nigericin vs P3C4") + # main title
  geom_text_repel(data = lab_P3C4_act, aes(label = Genes), force = 7, col = "black", size = 4, max.overlaps = 10, fontface = "bold") + # label poi
  geom_text_repel(data = lab_P3C4_act_up, aes(label = Genes), force = 7, col = "deeppink", size = 4, max.overlaps = 10, fontface = "bold") +
  geom_text_repel(data = lab_P3C4_act_down, aes(label = Genes), force = 7, col = "darkcyan", size = 4, max.overlaps = 10, fontface = "bold") +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA), # Background of the entire plot area
        legend.background = element_rect(fill = "transparent", color = NA), # Background of the legend
        legend.box.background = element_rect(fill = "transparent", color = NA))

print(p10)

ggsave( "VolcanoPlot10.png", plot = p10, width = 32, height = 20, units = c("cm"))
