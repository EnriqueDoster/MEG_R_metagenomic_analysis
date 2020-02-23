## Script to summarize counts
### Code for dataset analysis 
library(data.table)
library(tidyr)
#setwd("~/Dropbox/WRITING/NCBA2_JAN2018/NCBA2_analysis/") 
#source('scripts/Frontiers_NCBA2_analysis.R')
options(scipen = 999) # to decrease the use of scientific notation



setkey(amr_melted_raw_analytic,ID) 
setkey(amr_melted_analytic,ID) 
setkey(microbiome_melted_analytic,ID)
setkey(metadata,ID)
setkey(microbiome_metadata,ID)

amr_melted_analytic <- amr_melted_analytic[metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]

### Looking for clinically important AMR genes

amr_group_check <- amr_group_raw[!group %in% snp_regex, ]
important_AMR_regex = c('OXA',
                        'SME',
                        'sme',
                        'IMI',
                        'NDM',
                        'GES',
                        'KPC',
                        'CPHA',
                        'TEM',
                        'SHV',
                        'CTX',
                        'CMY',
                        'VGA',
                        'VGAB',
                        'VGAD',
                        'VATA',
                        'VATB',
                        'VATC',
                        'VATD',
                        'VATE',
                        'CFRA')

amr_raw_important_AMR <- amr_group_check[group %in% important_AMR_regex, ]

View(amr_raw_important_AMR[group %in% important_AMR_regex,  ])

melted_important_AMR <- amr_melted_raw_analytic[ Level_ID =='Group' & Name %in% important_AMR_regex,  ][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, Packaging_samples)]#[order(-sum_class )]
melted_important_AMR[,.(sum_important_genes = sum(Normalized_Count))]

ggplot(data = melted_important_AMR , aes(x = Packaging_samples, y = Name)) +
  geom_tile(aes(fill = log(Normalized_Count))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=5) 
ggsave("~/Dropbox/WRITING/FC_meat_2019/FC_meat_manuscript/FCmeat_figures/Supp_AMR_important_genes.jpeg", width = 30, height = 20, units = "cm")



##
###
#### Comparison of Sequencing read counts
###
##

Treatment_read_count_comparison <- metadata[,.(total_raw_reads = sum(Raw_paired_reads),mean_raw_reads = mean(Raw_paired_reads), sd_raw_reads = sd(Raw_paired_reads)), by=Treatment]

Packaging_read_count_comparison <- metadata[,.(total_raw_reads = sum(Raw_paired_reads),mean_raw_reads = mean(Raw_paired_reads), sd_raw_reads = sd(Raw_paired_reads), min_raw_reads = min(Raw_paired_reads),max_raw_reads = max(Raw_paired_reads)), by=Packaging]

ggplot(microbiome_metadata, aes(x=Packaging, y=Raw_paired_reads, color=Packaging)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1) , dotsize=1)
ggsave("~/Dropbox/WRITING/FC_meat_2019/FC_meat_manuscript/FCmeat_figures/Supp_FC_meat_Microbiome_raw_reads-by_Packaging.jpeg", width = 30, height = 20, units = "cm")

ggplot(metadata, aes(x=Packaging, y=Raw_paired_reads, color=Packaging)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1) , dotsize=1)
ggsave("~/Dropbox/WRITING/FC_meat_2019/FC_meat_manuscript/FCmeat_figures/Supp_FC_meat_AMR_raw_reads-by_Packaging.jpeg", width = 30, height = 20, units = "cm")




## AMR 
ggplot(metadata, aes(x=Packaging, y=Raw_paired_reads, color=Packaging)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1) , dotsize=1)+
  theme(
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  )



ggsave("FC_meat_AMR_raw_reads-by_Packaging.jpeg", width = 30, height = 20, units = "cm")


ggplot(metadata, aes(x=Dilution, y=Raw_paired_reads, color=Treatment)) +
  #geom_boxplot() + 
  geom_point(binaxis='y', stackdir='center',
             position=position_dodge(1) , dotsize=1)
ggplot(metadata, aes(x=sample, y=Raw_paired_reads, color=Dilution)) +
  geom_point(stat = 'identity')

ggsave("FC_meat_AMR_raw_reads-by_Packaging.jpeg", width = 30, height = 20, units = "cm")



##
### AMR majority of counts?
##

gene_counts = amr_melted_analytic[Level_ID=='Gene',]
setkey(gene_counts,Name)
setkey(annotations,header)
gene_counts <- gene_counts[annotations][ID != "NA"]
setkey(gene_counts,ID)

## Total AMR percentage for study
total_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(class) ][,total:=sum(class_sum)][,class_percentage:= class_sum/total * 100]
mech_percentage_byclass <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(mechanism,class) ][,class_percentage:= class_sum/sum(class_sum) * 100 , by = .(class)]

## Total AMR abundance by treatment group
AMR_counts_by_Treatment <- gene_counts[,.(class_sum = sum(Normalized_Count), mean_class = mean(Normalized_Count), sd_class = sd(Normalized_Count)), by= .(Treatment) ]
AMR_counts_by_Packaging <- gene_counts[,.(class_sum = sum(Normalized_Count), mean_class = mean(Normalized_Count), sd_class = sd(Normalized_Count)), by= .(Packaging) ]

total_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(Group, class) ][,class_percentage:= class_sum/sum(class_sum) * 100, by = Group]
total_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(Treatment, class) ][,class_percentage:= class_sum/sum(class_sum) * 100, by = Treatment]

#gene_counts <- gene_counts[metadata]
total_amr_percentage <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(class,mechanism) ][,class_percentage:= class_sum/sum(class_sum) * 100]
total_amr_percentage[, .(sum = sum(class_percentage))]

## AMR percentage of each mechanism per class, in complete dataset
mech_percentage_byclass <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(mechanism,class) ][,class_percentage:= class_sum/sum(class_sum) * 100 , by = .(class)]
mech_percentage_byclass[,.(test=sum(class_sum)), by=.(class, mechanism)]
## AMR classes by sample group
AMR_classes_percent_bygroup <- gene_counts[,.(class_sum = sum(Normalized_Count)), by= .(class,Group) ][,class_percentage:= class_sum/sum(class_sum) * 100 , by = .(Group)]






########################
##############
########
###
# MICROBIOME
###
########
##############
########################

## Microbiome and resistome Total, mean, std, and median counts per sample

# Microbiome
micro_count_summary_per_sample = microbiome_melted_analytic[Level_ID=="Phylum", .(total_counts = sum(Normalized_Count), mean_counts = mean(Normalized_Count), std_dev_counts = sd(Normalized_Count)), by = ID]
micro_count_summary_project = count_summary_per_sample[,.(total_reads = sum(total_counts),mean_per_sample=mean(total_counts), std_dev_counts= sd(total_counts))]



order_counts = microbiome_melted_analytic[Level_ID=='Order',]
setkey(order_counts,Name)


order_annotations = microbiome_taxonomy[,id := NULL]
setkey(order_annotations, Order)  # Data tables are SQL objects with optional primary keys
order_annotations <- unique(order_annotations, by = key(order_annotations))

setkey(order_annotations,Order)
order_counts <- order_counts[order_annotations][ID != "NA"]
setkey(order_counts,ID)

## Total AMR percentage for study
total_microbiome_percentage <- order_counts[,.(Order_sum = sum(Normalized_Count)), by= .(Name) ][,total:=sum(Order_sum)][,Order_percentage:= signif(Order_sum/total * 100,digits=3)]

microbiome_percentage_perPhylum <- order_counts[,.(sum_hits = sum(Normalized_Count)), by= .(Phylum, Class) ][,Class_percentage:= sum_hits/sum(sum_hits) * 100, by = .(Phylum)]

total_phylum <- microbiome_melted_analytic[Level_ID=="Phylum"][,.(sum_phylum = sum(Normalized_Count)), by = Name][,total := sum(sum_phylum)][,proportion := signif(sum_phylum/total * 100,digits=3) , by = Name]
total_class <- microbiome_melted_analytic[Level_ID=="Class"][,.(sum_class = sum(Normalized_Count)), by = Name][,total := sum(sum_class)][,proportion := signif(sum_class/total * 100,digits=3), by = Name]
total_order <- microbiome_melted_analytic[Level_ID=="Order"][,.(sum_order = sum(Normalized_Count)), by = Name][,total := sum(sum_order)][,proportion := signif(sum_order/total * 100,digits=3), by = Name]



########################
##############
########
###
# METADATA
###
########
##############
########################






## Average and standard deviations  ######
metadata_DT <- data.table(metadata)
metadata_DT[,.(mean_phred = mean(phred), sd_phred = sd(phred), min_phred = min(phred), max_phred = max(phred),
               mean_human_reads = mean(human_contaminant_reads), max_human_reads = max(human_contaminant_reads),
               min_human_reads = min(human_contaminant_reads), sum_human_reads = sum(human_contaminant_reads),
               total_reads = sum(as.numeric(Non_host_paired_reads)))][,.(percent_hum_contamination = sum_human_reads/total_reads * 100)]

metadata_DT[,.(mean_unclassified = mean(Unclassified_read_percentage), max_unclassified = max(Unclassified_read_percentage), min_unclassified = min(Unclassified_read_percentage)) ]


summarize_metadata_DT <- metadata_DT[, .(mean_phred = mean(phred), sd_phred = sd(phred), 
                                         mean_raw_reads = mean(Paired_reads), sd_raw_reads = sd(Paired_reads), 
                                         mean_16S_counts = mean(GG_16S_counts),sd_16S_counts = sd(GG_16S_counts), 
                                         mean_trimmed_reads = mean(Trimmed_paired_reads), sd_trimmed_reads = sd(Trimmed_paired_reads),
                                         mean_non_host = mean(Non_host_paired_reads), sd_non_host = sd(Non_host_paired_reads),
                                         mean_human_reads = mean(human_contaminant_reads), sd_human_reads = sd(human_contaminant_reads)),
                                     by = Group]

setkey(metadata,ID)
setkey(amr_melted_analytic,ID)
combined_melted <- amr_melted_analytic[metadata]
combined_melted <- combined_melted[Level_ID=="Group"]
#write.csv(combined_melted[Normalized_Count > 0, .(unique_mech = .N) , by = .(Group,Name)][unique_mech > 4], "unique_mech")


##### AMR exploratory #####

# count of each class by sample group


#amr_melted_raw_analytic[Level_ID=='Class'][Normalized_Count > 0 , .(sample_count = .N), by = .(Name, Group) ]

ggplot(amr_melted_raw_analytic[Level_ID=='Group'][Normalized_Count > 0 , .(sample_count = .N), by = .(Name, Group) ],
       aes(Group, Name , label=sample_count)) +
  geom_tile(aes(fill = sample_count), color = "white") +
  geom_text() + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of genes ") +
  xlab("List of patients") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16), 
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")


setkey(annotations,mechanism)

amr_melted <- amr_melted_raw_analytic[Level_ID=='Mechanism']
setkey(amr_melted,Name)
amr_melted <- annotations[amr_melted, allow.cartesian=TRUE]
amr_melted[is.na(Level_ID)]
# https://stackoverflow.com/questions/23197787/grouped-frequency-distribution-plot-in-r?rq=1
ggplot(amr_melted[Level_ID=='Mechanism'][Normalized_Count > 0 , .(sample_count = .N), by = .(mechanism, Group,class) ],
       aes(Group, interaction(class, mechanism) , label=sample_count)) +
  geom_tile(aes(fill = sample_count), color = "white") +
  #coord_flip()+
  #geom_text() +
  #facet_grid(~class) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of genes ") +
  xlab("List of patients") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")

amr_melted_raw_analytic[Level_ID=='Class'][Normalized_Count > 0 , .(sample_count = .N), by = .(Name, Group,Time,Treatment) ] %>%
  ggplot(aes(x=Time,y=sample_count, color=Name, group=Group)) +
  #geom_line(aes(x=Time,y=sample_count,linetype=Treatment))+
  geom_line()+
  geom_point()+
  #facet_wrap(~Cow_ID) +
  scale_fill_manual(name = "Name",values = myColors) +
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
  ggtitle('Total AMR Gene Abundance at the class level, for each treatment group at arrival and day 11 \n')

#########
# increase in abundance
amr_melted_raw_analytic[Level_ID=='Class', .(AMR_sum = sum(Normalized_Count)), by = .(Group,Time,Treatment) ] 

########
#######
#length(unique(microbiome_melted_analytic[Level_ID=="Phylum"]$Name))






#### AMR important genes #########
#source('scripts/print_results.R')
amr_melted_analytic[[1]]
amr_melted_analytic[contains("OXA"),]


amr_raw_SNP <- amr_raw[!group %in% snp_regex, ]

important_AMR_regex = c('OXA',
                        'SME',
                        'sme',
                        'IMI',
                        'NDM',
                        'GES',
                        'KPC',
                        'CPHA',
                        'TEM',
                        'SHV',
                        'CTX',
                        'CMY',
                        'VGA',
                        'VGAB',
                        'VGAD',
                        'VATA',
                        'VATB',
                        'VATC',
                        'VATD',
                        'VATE',
                        'CFRA')

amr_raw_important_AMR <- amr_raw[group %in% important_AMR_regex, ]
metadata["X156p2_GTTTCG_L006"]

View(amr_raw[group %in% important_AMR_regex,  ])






##### A FEW MORE SUMMARY FIGURES #####
#### AMR figures
###
##

## By ID, faceted by packaging, Class level
AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(ID, Name) ]
AMR_class_sum$Class <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_class_sum, aes(x = ID, y = percentage, fill = Class)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("CSS counts -resistome composition") +
  xlab('Sample ID') +
  ylab('Relative abundance') 

ggsave("FC_meat_CSS_AMR_barplot-RGISNPconfirmed_deduped_amrplusplus.jpeg", width = 30, height = 20, units = "cm")




AMR_class_sum_raw <- amr_melted_raw_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
AMR_class_sum_raw[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum_raw[,percentage:= sum_class/total ,by=.(ID, Name) ]
AMR_class_sum_raw$Class <- AMR_class_sum_raw$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_class_sum_raw, aes(x = ID, y = percentage, fill = Class)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("Raw counts -resistome composition") +
  xlab('\n Sample ID') +
  ylab('Relative abundance\n') 

ggsave("FC_meat_raw_AMR_barplot-RGISNPconfirmed_deduped_amrplusplus.jpeg", width = 30, height = 20, units = "cm")



noSNP_CSS_amr_analytic_data <- amr_melted_analytic
noSNP_CSS_amr_analytic_data$Type <- "noSNP_samdedup"
#write.csv(noSNP_CSS_amr_analytic_data,"FC_meat_noSNP_samdedup_amrplusplus.csv", row.names = FALSE)


### Microbiome summary barplots

microbiome_phylum_sum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Packaging )]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]
# only keep taxa greater than 1%
microbiome_phylum_sum <- microbiome_phylum_sum[percentage > .01]
microbiome_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
microbiome_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]

microbiome_phylum_sum$Class <- microbiome_phylum_sum$Name
#microbiome_phylum_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(microbiome_phylum_sum, aes(x = ID, y = percentage, fill = Class)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging, scales='free',ncol = 2) +
  #scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("Microbiome composition in Fort Collins ground \n beef, by packaging type (only taxa > 1% of sample)") +
  xlab('Sample ID') +
  ylab('Relative abundance') 

ggsave("FC_meat_CSS_microbiome_barplot-noSNP_samtools_dedup_amrplusplus.jpeg", width = 30, height = 20, units = "cm")



