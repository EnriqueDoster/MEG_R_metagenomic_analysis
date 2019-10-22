## NCBA2 - Salmonella identification 
setwd("~/Dropbox/WRITING/NCBA2_JAN2018/NCBA2_analysis") 
library(dplyr)
library(ggthemes)
source('scripts/NCBA2_analysis.R')

### Looking at total hits for filtered standard and modified database
total_hits_std <- kraken_melted_raw_analytic[Level_ID=="Domain", .(Reads = sum(Normalized_Count)), by=.(ID)]

total_hits_std$Database = 'Standard database'
total_hits_custom$Database = 'Modified database'

kraken_v2_hits_results = rbind(total_hits_std,total_hits_custom)
write.csv(kraken_v2_hits_results,"Kraken_classification_totals.csv",row.names = FALSE)

kraken_v2_hits_results$Database = factor(kraken_v2_hits_results$Database ,levels=c("Standard database","Modified database"))

kraken_v2_hits_results %>%
ggplot(aes(x=reorder(ID,-Reads), y = Reads, group=rev(Database), color=Database, fill=Database)) +
  geom_bar(stat='identity', position= position_dodge(width=0))+
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
    legend.text=element_text(size=16),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  scale_colour_tableau() +
  scale_y_continuous(breaks= seq(0, 900000, by = 50000)) +
  xlab('Samples') +
  ylab('Number of reads') +
  ggtitle('Total count of classified reads in each sample, by database used')
ggsave("~/Dropbox/WRITING/NCBA2_JAN2018/Salmonella_manuscript/Figures/Fig1_total_classified_reads.jpeg", width = 40, height = 30, units = "cm")


total_hits_std$Database = 'Standarddatabase'
total_hits_custom$Database = 'Modifieddatabase'

kraken_v2_hits_results = rbind(total_hits_std,total_hits_custom)


wide_kraken2_results = dcast(kraken_v2_hits_results, ID ~ Database , value.var="Reads")
write.csv(wide_kraken2_results,"Kraken_classification_totals_wide.csv",row.names = FALSE)

percen_kraken = wide_kraken2_results[,.(percent_plasmid = 100 - (Modifieddatabase/Standarddatabase * 100)), by= .(ID)]

### Subsetting Salmonella Counts
Salmonella_custom_filter_v2 
Salmonella_standard_filter_v2 <- kraken_melted_raw_analytic[Name=="Salmonella enterica" & Normalized_Count > 0]


Salmonella_standard_v2$Type = "Salmonella_standard_v2"
Salmonella_standard_filter_v2$Type = "Salmonella_standard_filter_v2"
Salmonella_custom_v2$Type = "Salmonella_custom_v2"

Salmonella_custom_filter_v2$Type = "Salmonella_custom_filter_v2"


Salmonella_kraken_v2_results = rbind(Salmonella_custom_filter_v2,Salmonella_custom_v2,Salmonella_standard_filter_v2,Salmonella_standard_v2)
write.csv(Salmonella_kraken_v2_results,"Salmonella_kraken_v2_results.csv",row.names = FALSE)
#Salmonella_kraken_v2_results = read.csv("Salmonella_kraken_v2_results.csv", header = TRUE)



S_id <- metadata_play[,.(ID, PCR,Senterica_counts_custom_filtered,Senterica_counts_std_filter,Senterica_counts_std,Senterica_counts_custom)]

S_id_long_hits <-  melt(S_id, id.vars = c("ID","PCR"),variable.name="count", measure.vars = c("Senterica_counts_std","Senterica_counts_custom","Senterica_counts_custom_filtered","Senterica_counts_std_filter"))

S_id_long_hits %>% 
  ggplot(aes(x=reorder(ID,-value), y = value, group=count, color=count)) +
  geom_point(aes(fill=count))+
  geom_line() +
  #facet_wrap(~PCR) + #, scales='free'
   theme(
     panel.grid.major=element_blank(),
     panel.grid.minor=element_blank(),
     panel.background=element_blank(),
     axis.title.x=element_blank(),
     strip.text.x=element_blank(),
     axis.text.x=element_blank(),
     axis.ticks.x=element_blank(),
     strip.text.y=element_text(size=24, angle=0),
     axis.text.y=element_text(size=22),
     axis.title=element_text(size=26),
     legend.position="right",
     panel.spacing=unit(0.1, "lines"),
     plot.title=element_text(size=32, hjust=0.5),
     legend.text=element_text(size=18),
     legend.title=element_blank(),
     legend.key=element_blank()
   ) +
  scale_colour_tableau() +
  #theme_tufte(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  xlab('\nSample') +
  ylab('Count of raw observations\n') +
  ggtitle('Total counts of Salmonella enterica, by sample \n')


## Number positive
as.data.table(Salmonella_kraken_v2_results)[Type=="Salmonella_standard_filter_v2" | Type=="Salmonella_custom_filter_v2"][,.N,by=Type]

# Graph with sample on x axis
as.data.table(Salmonella_kraken_v2_results)[Type=="Salmonella_standard_filter_v2" | Type=="Salmonella_custom_filter_v2"] %>% 
  ggplot(aes(x=reorder(ID,-Normalized_Count), y = Normalized_Count, group=Type, color=Type)) +
  geom_point(aes(fill=Normalized_Count))+
  geom_line() +
  #facet_wrap(~PCR) + #, scales='free'
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    axis.title.x=element_blank(),
    strip.text.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank(),
    legend.key=element_blank()
  ) +
  scale_colour_tableau() +
  #theme_tufte(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  xlab('\nSample') +
  ylab('Count of raw observations\n') +
  ggtitle('Total counts of Salmonella enterica, by sample \n')
