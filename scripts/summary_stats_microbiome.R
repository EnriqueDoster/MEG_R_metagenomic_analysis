library(data.table)
library(tidyr)
library(ggthemes)

setkey(microbiome_melted_raw_analytic,ID) 
setkey(microbiome_melted_analytic,ID) 
setkey(microbiome_metadata,ID)

microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
microbiome_melted_raw_analytic <- microbiome_melted_raw_analytic[microbiome_metadata]


#
##
### MEGARes - gene level relative abundance - raw counts
## 
## With label for "low abundance classes"
Microbiome_total_Phylum <- microbiome_melted_raw_analytic[Level_ID=="Phylum", .(sum_class= sum(Normalized_Count)),by=.(Name)]
Microbiome_total_Phylum[,total:= sum(sum_class)]
Microbiome_total_Phylum[,proportion:= sum_class/total , by=.(Name) ]

length(unique(Microbiome_total_Phylum$Name))
rare_class <- Microbiome_total_Phylum[proportion < .001 ,Name]
length(rare_class)

Microbiome_phylum_by_ID_edited <- microbiome_melted_raw_analytic[Level_ID=="Phylum", .(sum_class = sum(Normalized_Count),Group,Unique_sample), by = .(Name,ID)]

Microbiome_phylum_by_ID_edited[(Name %in% rare_class), Name := 'Low abundance <0.1% Phyla']
#Microbiome_phylum_by_ID_edited[(Name %in% rare_class),]

Microbiome_phylum_by_ID_edited[,sample_total:= sum(sum_class), by=.(ID)]
Microbiome_phylum_by_ID_edited[, percentage := sum_class/sample_total *100 ,by=.( Name, ID) ]
Microbiome_phylum_by_ID_edited[, percentage_label:= as.character(percentage)]
Microbiome_phylum_by_ID_edited[percentage_label=='0', percentage_label := '<0.01']

# Drop extra factors 
Microbiome_phylum_by_ID_edited$Name = droplevels(Microbiome_phylum_by_ID_edited$Name)
unique(factor(Microbiome_phylum_by_ID_edited$Name))
# Check which factors you have and change the order here
Microbiome_phylum_by_ID_edited$Name = factor(Microbiome_phylum_by_ID_edited$Name ,levels=c("Low abundance <0.1% Phyla","GN02","", "TM7" ,"Tenericutes",
                                                                           "Verrucomicrobia","Planctomycetes","Actinobacteria","Phenicol",
                                                                           "Chloroflexi","Proteobacteria","Firmicutes"))

Microbiome_phylum_by_ID_edited$Phyla <- Microbiome_phylum_by_ID_edited$Name

ggplot(Microbiome_phylum_by_ID_edited[percentage >0], aes(x = ID, y = percentage, fill = Phyla)) + 
  geom_bar(stat = "identity") +
  facet_wrap( ~ Unique_sample, scales='free_x',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=12),
    legend.title=element_text(size=16),
    panel.background = element_rect(fill = "white")
  ) +
  #ggtitle("CSS counts -resistome composition") +
  xlab('Sample ID') +
  ylab('Relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1)

ggsave("graphs/Raw_counts_relative_abundance_Microbiome_Phylum_byGroup.jpeg", width = 60, height = 30, units = "cm")



ggplot(Microbiome_phylum_by_ID_edited[percentage >0], aes(x = ID, y = sum_class, fill = Class)) + 
  geom_bar(stat = "identity") +
  facet_wrap( ~ Group, scales='free_x',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
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
    legend.text=element_text(size=12),
    legend.title=element_text(size=16),
    panel.background = element_rect(fill = "white")
  ) +
  #ggtitle("CSS counts -resistome composition") +
  xlab('Sample ID') +
  ylab('counts') +
  scale_fill_tableau("Tableau 20", direction = -1)

ggsave("graphs/Raw_counts_Microbiome_Phylum_byGroup.jpeg", width = 60, height = 30, units = "cm")



