### final figures
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
#source('scripts/NCBA2_anosim.r')

### Start of code for figures

setkey(amr_melted_raw_analytic,ID) 
setkey(amr_melted_analytic,ID) 

setkey(microbiome_melted_analytic,ID)

# Set keys for both metadata files
setkey(metadata,ID)
setkey(microbiome_metadata,ID)

microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
amr_melted_analytic <- amr_melted_analytic[metadata]

#### AMR figures
###
##
# Heatmap at group level - by sample
AMR_group_samples <- amr_melted_raw_analytic[Level_ID=="Group"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Group)]#[order(-sum_class )]

group_annotations <- annotations[,-1]
group_annotations <- unique(group_annotations,by = c("group"))

setkey(group_annotations, group)
setkey(AMR_group_samples , Name)

AMR_group_samples  <- group_annotations[AMR_group_samples ]
setkey(AMR_group_samples , Group)
AMR_group_samples  <- metadata[AMR_group_samples]

ggplot(data = AMR_group_samples , aes(x = Group, y = group)) +
  geom_tile(aes(fill = log_Normalized_Count)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=2.5) 



#write.csv(AMR_group_samples, "16S_norm_group_melted_sample_counts_.csv")


## Heatmap at Mechanism level - by sample
AMR_mech_samples <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, ID)]#[order(-sum_class )]

mech_annotations <- annotations[,-1]
mech_annotations <- unique(mech_annotations,by = c("mechanism"))

setkey(mech_annotations, mechanism)
setkey(AMR_mech_samples , Name)

AMR_mech_samples  <- mech_annotations[AMR_mech_samples ]
setkey(AMR_mech_samples , ID)
AMR_mech_samples  <- metadata[AMR_mech_samples]
#write.csv(AMR_mech_samples, "16S_norm_melted_sample_counts.csv")



## Heatmap at Group level
AMR_mech_sample_groups <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Group)]#[order(-sum_class )]
AMR_mech_sample_groups$Group <- factor(AMR_mech_sample_groups$Group,levels = c("Arrival", "Exit", "Shipment"))

mech_annotations <- annotations[,-1]
mech_annotations <- unique(mech_annotations,by = c("mechanism"))

setkey(mech_annotations, mechanism)
setkey(AMR_mech_sample_groups, Name)

AMR_mech_sample_groups <- mech_annotations[AMR_mech_sample_groups]

write.csv(AMR_mech_sample_groups, "16S_norm_melted_Group_counts.csv")


ggplot(data = AMR_mech_sample_groups , aes(x = Group, y = mechanism)) +
  geom_tile(aes(fill = log_Normalized_Count)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=2.5) 


setkey(AMR_mech_samples, ID)
AMR_mech_samples <- AMR_mech_samples[metadata]



## Stacked bars - average sum of AMR class
AMR_class_mean <- amr_melted_analytic[Level_ID=="Class", .(sum_class= mean(Normalized_Count)),by=.(Group, Name)][order(-sum_class )]
ggplot(AMR_class_mean, aes(x = Group, y = sum_class, fill = Name)) + 
  geom_bar(stat = "identity")

##
## Stacked 100% bars
##
##
AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(DOF_categories, Name)][order(-sum_class )]
AMR_class_sum[,total:= sum(sum_class), by=.(DOF_categories)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(DOF_categories, Name) ]
AMR_class_sum$Class <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(DOF_categories, Name) ] removes some with low proportions
ggplot(AMR_class_sum, aes(x = DOF_categories, y = percentage, fill = Class)) + 
  geom_bar(stat = "identity")+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample DOF_categories') +
  ylab('Relative abundance\n') 




### Figure 1 - AMR Class stacked graph ####

AMR_class_sum <- amr_melted_raw_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(Group, Name)][order(-sum_class )]
AMR_class_sum$Name = droplevels(AMR_class_sum$Name)
AMR_class_sum <- AMR_class_sum[with(AMR_class_sum, order(-sum_class)), ]
#AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=sort(levels(AMR_class_sum$Name ), FALSE))
#AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=c("Fluoroquinolones","Phenicol","Bacitracin","Cationic antimicrobial peptides", "Multi-drug resistance" ,"Aminoglycosides",   "betalactams" ,"MLS" ,"Tetracyclines"))
AMR_class_sum[,total:= sum(sum_class), by=.(Group)]
AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(Group, Name) ]
AMR_class_sum[percentage < .025,percentage:=as.character('')]
#AMR_class_sum[percentage > 2.5,percentage:= paste0(percentage,"%")]

#myColors <- rev(brewer.pal(14,"Set1"))
#names(myColors) <- levels(AMR_class_sum$Name)
#colScale <- scale_colour_manual(name = "Name",values = myColors)
AMR_class_sum$Group <- factor(AMR_class_sum$Group,levels = c("Arrival", "Exit", "Shipment"))


fig1 <- AMR_class_sum %>% 
  ggplot(aes(x=Group,y=sum_class, fill=Name, label= percentage)) +
  geom_bar(stat='identity')+
  #facet_wrap(~Cow_ID) +
  #scale_fill_manual(name = "Name",values = myColors) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample group') +
  ylab('16S normalized AMR gene abundance\n') 
#ggtitle('            Total AMR Gene Abundance and proportion at the class level, by sample group\n')
fig1
ggsave("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure1_AMR_barplot.jpeg", width = 30, height = 20, units = "cm")
#dev.off()
### By average of AMR class counts

AMR_class_mean <- amr_melted_raw_analytic[Level_ID=="Class", .(mean_class= mean(Normalized_Count)),by=.(Group, Name)]
AMR_class_mean$Name = droplevels(AMR_class_mean$Name)
AMR_class_mean <- AMR_class_mean[with(AMR_class_mean, order(-mean_class)), ]
AMR_class_mean$Name = factor(AMR_class_mean$Name ,levels=c("Phenicol","Bacitracin","Cationic antimicrobial peptides", "Multi-drug resistance" ,"Aminoglycosides",   "betalactams" ,"MLS" ,"Tetracyclines"))
AMR_class_mean[,total:= sum(mean_class), by=.(Group)]
AMR_class_mean[,percentage:= round(mean_class/total, digits=2) ,by=.(Group, Name) ]
AMR_class_mean[percentage < .03,percentage:='']

fig_mean <- AMR_class_mean %>% 
  ggplot(aes(x=Group,y=mean_class,fill=Name,label= percentage)) +
  geom_bar(stat='identity')+
  #facet_wrap(~Cow_ID) +
  scale_fill_manual(name = "Name",values = myColors) +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank()
  ) +
  xlab('\n Sample group') +
  ylab('16S normalized AMR gene abundance\n') +
  ggtitle('Average AMR Gene Abundance and proportion at the class level, for each treatment group at Day 1 and day 11 \n')
fig_mean


## Figure 2 -AMR Group ordination #####
jpeg("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure2_AMR_Group.jpeg")

par(mfrow = c(2,1))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots 
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

##
plot(metaMDS_AMR_class, type="none", display=c("sites"))
points(metaMDS_AMR_class, display="sites", pch=20, col=as.numeric(pData(AMR_analytic_data[[1]])$Group_cat))
groupz <- sort(unique(pData(AMR_analytic_data[[1]])$Group_cat))
for(i in seq(groupz)) {ordispider(metaMDS_AMR_class, pData(AMR_analytic_data[[1]])$Group_cat,font=2, cex=1.5, col=i, show.groups=groupz[i])}

#mtext(side = 3, "AMR Class", line = 1, cex=2)
#mtext(side = 3, "Resistome Ordination", line = 1, cex=2)
mtext(side = 2, "AMR Class", line = 4, cex=2)

plot(metaMDS_AMR_mech, type="none", display=c("sites"))
points(metaMDS_AMR_mech, display="sites", pch=20, col=as.numeric(pData(AMR_analytic_data[[2]])$Group_cat))
groupz <- sort(unique(pData(AMR_analytic_data[[2]])$Group_cat))
for(i in seq(groupz)) {ordispider(metaMDS_AMR_mech, pData(AMR_analytic_data[[2]])$Group_cat,font=2, cex=1.5, col=i, show.groups=groupz[i])}

mtext(side = 2, "AMR Mechanism", line = 4, cex=2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Day 1-Treated","Day 1-Untreated","Day11-Treated","Day11-Untreated"), xpd = TRUE, horiz = TRUE, inset = c(0, 
                                                                                                                             0), bty = "n", pch = c(19), col = 1:4, cex = 1.5)
#identify(fig,"sites")
dev.off()

##Figure 3 Log Fold change in AMR abundance over time?


#### Figure 4 - AMR richness and diversity ##########
jpeg("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure4_AMR_richness_diversity_boxplots.jpeg")

par(mfrow = c(2, 2))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots 
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

boxplot(metadata$AMR_class_Richness ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
mtext(side = 3, "AMR Class", line = 1, cex=2)
mtext(side = 2, "Richness", line = 4, cex=2)

boxplot(metadata$AMR_mech_Richness ~ metadata$Group_cat,col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
#title("AMR mechanism richness by group")
mtext(side = 3, "AMR Mechanism", line = 1, cex=2)

boxplot(metadata$AMR_class_Shannon ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
mtext(side = 2, "Shannon\'s Diversity", line = 4, cex=2)

#title("AMR Class diversity by group")
boxplot(metadata$AMR_mech_Shannon ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
#title("AMR mechanism diversity by group")


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Day 1-Treated","Day 1-Untreated", "Day11-Treated","Day11-Untreated"),x.intersp = .5, xpd = TRUE, horiz = TRUE, inset = c(9, 
                                                                                                                                             0), bty = "n", pch = 15, fill = c("Red","White","deepskyblue3","Grey"), cex = 1.5)


dev.off()






### Figure 5 - Microbiome Phylum stacked graph ####

fig1 <- microbiome_melted_raw_analytic %>% 
  filter(Level_ID=="Class") %>%
  ggplot(aes(x=reorder(Pair_description))) +
  geom_bar(aes(fill=factor(Name), weight=Normalized_Count))+
  #facet_wrap(~Cow_ID) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample type') +
  ylab('Count of raw observations\n') +
  ggtitle('raw AMR Mechanism hits - by cow ID \n')
fig1


microbiome_phylum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(Group, Name)]
microbiome_phylum[,total:= sum(sum_phylum), by=.(Group)]
microbiome_phylum[,proportion:= sum_phylum/total ,by=.(Group, Name) ]
rare_phyla <- microbiome_phylum[proportion < .03,Name]

## Changing name of 
microbiome_phylum_edited <- microbiome_melted_analytic[Level_ID=="Phylum", .(Normalized_Count,Group),by=.(Name,ID)]
microbiome_phylum_edited[Name %in% rare_phyla ,Name := 'Low abundance phyla']
microbiome_phylum_edited[,total:= sum(Normalized_Count), by=.(Group)]

microbiome_phylum_mean <- microbiome_phylum_edited[, .(mean_phylum = mean(Normalized_Count)) ,by=.(Group,Name) ]
microbiome_phylum_mean[, total := sum(mean_phylum) ,by=.(Group) ]
microbiome_phylum_mean[, percentage := round(mean_phylum/total,digits=2) ,by=.(Group, Name) ]
microbiome_phylum_mean[, percentage:= as.character(percentage)]
microbiome_phylum_mean[percentage=='0', percentage := '<0.01']

microbiome_phylum_mean$Name = droplevels(microbiome_phylum_mean$Name)
#microbiome_phylum_mean$Name = rev(reorder(microbiome_phylum_mean$Name, X=as.numeric(microbiome_phylum_mean$phylum), FUN=sum))
microbiome_phylum_mean$Name = factor(microbiome_phylum_mean$Name ,levels=c("Low abundance phyla","Actinobacteria","Proteobacteria","Bacteroidetes","Firmicutes"))

myColors <- rev(brewer.pal(5,"Set3"))
names(myColors) <- levels(microbiome_phylum_mean$Name)
colScale <- scale_colour_manual(name = "Name",values = myColors)

microbiome_phylum_mean$Group <- factor(microbiome_phylum_mean$Group,levels = c("Day1-Treated", "Day1-Untreated", "Day11-Treated", "Day11-Untreated"))


fig1 <- microbiome_phylum_mean %>% 
  ggplot(aes(x=Group,y=mean_phylum, fill=Name, label= percentage)) +
  geom_bar(stat='identity')+
  #facet_wrap(~Cow_ID) +
  scale_fill_manual(name = "Name",values = myColors) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample group') +
  ylab('CSS normalized counts\n') 
#ggtitle('Average phyla counts and proportion, by sample group \n')
fig1

ggsave("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/Review/Review_submission_files/figures/Figure5_Microbiome_barplot_update.jpeg", width = 40, height = 30, units = "cm")





### Figure 6 - Microbiome group ordination ####

#pdf("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure6_Microbiome_Group.pdf")
jpeg("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure6_Microbiome_Group.jpeg")

par(mfrow = c(3,1))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots 
par(tcl = -0.5)
par(mgp = c(2, 0.6, 0))


plot(metaMDS_microbiome_phylum, type="none", display=c("sites"))
points(metaMDS_microbiome_phylum, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[2]])$Group_cat))
groupz <- sort(unique(pData(microbiome_analytic_data[[2]])$Group_cat))
for(i in seq(groupz)) {ordispider(metaMDS_microbiome_phylum, pData(microbiome_analytic_data[[2]])$Group_cat,font=2, cex=1.5, col=i, show.groups=groupz[i])}

#mtext(side = 3, "Microbiome Ordination", line = 1, cex=2)
mtext(side = 2, "Phylum", line = 4, cex=2)

plot(metaMDS_microbiome_class, type="none", display=c("sites"))
points(metaMDS_microbiome_class, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[3]])$Group_cat))
groupz <- sort(unique(pData(microbiome_analytic_data[[3]])$Group_cat))
for(i in seq(groupz)) {ordispider(metaMDS_microbiome_class, pData(microbiome_analytic_data[[3]])$Group_cat,font=2, cex=1.5, col=i, show.groups=groupz[i])}

mtext(side = 2, "Class", line = 4, cex=2)

plot(metaMDS_microbiome_order, type="none", display=c("sites"))
points(metaMDS_microbiome_order, display="sites", pch=20, col=as.numeric(pData(microbiome_analytic_data[[4]])$Group_cat))
groupz <- sort(unique(pData(microbiome_analytic_data[[4]])$Group_cat))
for(i in seq(groupz)) {ordispider(metaMDS_microbiome_order, pData(microbiome_analytic_data[[4]])$Group_cat,font=2, cex=1.5, col=i, show.groups=groupz[i])}

mtext(side = 2, "Order", line = 4, cex=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Day 1-Treated","Day 1-Untreated","Day11-Treated","Day11-Untreated"),xpd = TRUE, horiz = TRUE, inset = c(0, 
                                                                                                                            0), bty = "n", pch = c(19), col = 1:4, cex = 1.5)
#identify(fig,"sites")
dev.off()




##### Figure 7 , Micro richness and diversity ######
#pdf("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure7_Micro_richness_diversity.pdf")
jpeg("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure7_Micro_richness_diversity.jpeg")

par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(3, 4, 1, 1))
par(oma = c(4,3,4,1)) ## moves the plot area, not for individual plots 
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

boxplot(metadata$microbiome_phylum_Richness ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
mtext(side = 3, "Phylum", line = 1, cex=2)
mtext(side = 2, "Richness", line = 4, cex=2)
boxplot(metadata$microbiome_class_Richness ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
mtext(side = 3, "Class", line = 1, cex=2)
boxplot(metadata$microbiome_order_Richness ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
mtext(side = 3, "Order", line = 1, cex=2)
boxplot(metadata$microbiome_phylum_Shannon ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
mtext(side = 2, "Shannon's Diversity", line = 4, cex=2)
boxplot(metadata$microbiome_class_Shannon ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)
boxplot(metadata$microbiome_order_Shannon ~ metadata$Group_cat, col=c("Red","White","deepskyblue3","Grey"),names=FALSE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Day 1-Treated","Day 1-Untreated","Day11-Treated","Day11-Untreated"),x.intersp = .5, xpd = TRUE, horiz = TRUE, inset = c(9, 
                                                                                                                                            0), bty = "n", pch = 15, fill = c("Red","White","deepskyblue3","Grey"), cex = 1.5)
dev.off()
