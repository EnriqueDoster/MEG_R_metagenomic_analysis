##wilcox analytical tests
#source('scripts/NCBA2_anosim.r')
boxplot(pData(AMR_raw_analytic_data[[1]])$AMR_mech_Shannon ~ pData(AMR_raw_analytic_data[[1]])$Treatment)
boxplot(pData(AMR_raw_analytic_data[[1]])$kraken_phylum_Richness ~ pData(AMR_raw_analytic_data[[1]])$Treatment)

## GLM for paired reads, trimmed reads, nonhost. tested by batch, lane, run
summary(glm(metadata$Raw_paired_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$Raw_paired_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$Raw_paired_reads ~ as.factor(metadata$Packaging))) 

summary(glm(metadata$QC_filtered_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$QC_filtered_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$QC_filtered_reads ~ as.factor(metadata$Packaging)))

summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Packaging)))

summary(glm(metadata$deduped_SNP_confirmed_counts ~ as.factor(metadata$Treatment)))
summary(glm(metadata$deduped_SNP_confirmed_counts ~ as.factor(metadata$Dilution)))
summary(glm(metadata$deduped_SNP_confirmed_counts ~ as.factor(metadata$Packaging)))

## GLM microbome
summary(glm(microbiome_metadata$Raw_paired_reads ~ as.factor(microbiome_metadata$Treatment)))
summary(glm(microbiome_metadata$Raw_paired_reads ~ as.factor(microbiome_metadata$Packaging))) 

summary(glm(microbiome_metadata$QC_filtered_reads ~ as.factor(microbiome_metadata$Treatment)))
summary(glm(microbiome_metadata$QC_filtered_reads ~ as.factor(microbiome_metadata$Packaging)))

summary(glm(microbiome_metadata$nonhost_filtered_reads ~ as.factor(metadata$Treatment)))
summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Dilution)))
summary(glm(metadata$nonhost_filtered_reads ~ as.factor(metadata$Packaging)))



### Wilcox testing of diversity 

#### AMR diversity comparisons 
###
##
#
## Test by treatment, Class and Mech
wilcox.test(pData(AMR_analytic_data[[1]])$AMR_class_Richness ~ pData(AMR_analytic_data[[1]])$Treatment)
wilcox.test(pData(AMR_analytic_data[[1]])$AMR_class_Shannon ~ pData(AMR_analytic_data[[1]])$Treatment)
wilcox.test(pData(AMR_analytic_data[[2]])$AMR_class_Richness ~ pData(AMR_analytic_data[[2]])$Treatment)
wilcox.test(pData(AMR_analytic_data[[2]])$AMR_class_Shannon ~ pData(AMR_analytic_data[[2]])$Treatment)

## Test by Dilution, Class and Mech
wilcox.test(pData(AMR_analytic_data[[1]])$AMR_class_Richness ~ pData(AMR_analytic_data[[1]])$Dilution)
wilcox.test(pData(AMR_analytic_data[[1]])$AMR_class_Shannon ~ pData(AMR_analytic_data[[1]])$Dilution)
wilcox.test(pData(AMR_analytic_data[[2]])$AMR_class_Richness ~ pData(AMR_analytic_data[[2]])$Dilution)
wilcox.test(pData(AMR_analytic_data[[2]])$AMR_class_Shannon ~ pData(AMR_analytic_data[[2]])$Dilution)


summary(glm(pData(AMR_analytic_data[[2]])$AMR_class_Richness ~ pData(AMR_analytic_data[[2]])$Packaging))
summary(glm(pData(AMR_analytic_data[[2]])$AMR_class_Shannon ~ pData(AMR_analytic_data[[2]])$Packaging))

summary(glm(pData(AMR_analytic_data[[2]])$AMR_class_Richness ~ pData(AMR_analytic_data[[2]])$Blinded_Store))
summary(glm(pData(AMR_analytic_data[[2]])$AMR_class_Shannon ~ pData(AMR_analytic_data[[2]])$Blinded_Store))


#### Microbiome diversity comparisons 
###
##
#
## Test by treatment, Class and Mech
wilcox.test(pData(microbiome_analytic_data[[1]])$microbiome_phylum_Richness ~ pData(microbiome_analytic_data[[1]])$Treatment)
wilcox.test(pData(microbiome_analytic_data[[1]])$microbiome_phylum_Shannon ~ pData(microbiome_analytic_data[[1]])$Treatment)

wilcox.test(pData(microbiome_analytic_data[[1]])$microbiome_class_Richness ~ pData(microbiome_analytic_data[[1]])$Treatment)
wilcox.test(pData(microbiome_analytic_data[[1]])$microbiome_class_Shannon ~ pData(microbiome_analytic_data[[1]])$Treatment)

wilcox.test(pData(microbiome_analytic_data[[1]])$microbiome_order_Richness ~ pData(microbiome_analytic_data[[1]])$Treatment)
wilcox.test(pData(microbiome_analytic_data[[1]])$microbiome_order_Shannon ~ pData(microbiome_analytic_data[[1]])$Treatment)



summary(glm(pData(microbiome_analytic_data[[1]])$microbiome_phylum_Richness ~ pData(microbiome_analytic_data[[1]])$Packaging))
summary(glm(pData(microbiome_analytic_data[[1]])$microbiome_phylum_Shannon ~ pData(microbiome_analytic_data[[1]])$Packaging))

summary(glm(pData(microbiome_analytic_data[[2]])$microbiome_class_Richness ~ pData(microbiome_analytic_data[[2]])$Packaging))
summary(glm(pData(microbiome_analytic_data[[2]])$microbiome_class_Shannon ~ pData(microbiome_analytic_data[[2]])$Packaging))

summary(glm(pData(microbiome_analytic_data[[2]])$microbiome_class_Richness ~ pData(microbiome_analytic_data[[2]])$Packaging))
summary(glm(pData(microbiome_analytic_data[[2]])$microbiome_class_Shannon ~ pData(microbiome_analytic_data[[2]])$Packaging))




summary(glm(pData(microbiome_analytic_data[[1]])$microbiome_phylum_Richness ~ pData(microbiome_analytic_data[[1]])$Blinded_Store))
summary(glm(pData(microbiome_analytic_data[[1]])$microbiome_phylum_Shannon ~ pData(microbiome_analytic_data[[1]])$Blinded_Store))

summary(glm(pData(microbiome_analytic_data[[2]])$microbiome_class_Richness ~ pData(microbiome_analytic_data[[2]])$Blinded_Store))
summary(glm(pData(microbiome_analytic_data[[2]])$microbiome_class_Shannon ~ pData(microbiome_analytic_data[[2]])$Blinded_Store))

summary(glm(pData(microbiome_analytic_data[[3]])$microbiome_order_Richness ~ pData(microbiome_analytic_data[[3]])$Blinded_Store))
summary(glm(pData(microbiome_analytic_data[[3]])$microbiome_order_Shannon ~ pData(microbiome_analytic_data[[3]])$Blinded_Store))







wilcox.test(pData(AMR_class_arrival)$AMR_mech_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$AMR_mech_Shannon ~ pData(AMR_class_arrival)$Treatment)
## Test by treatment, day11
wilcox.test(pData(AMR_class_day11)$AMR_class_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$AMR_class_Shannon ~ pData(AMR_class_day11)$Treatment)

boxplot(pData(AMR_class_day11)$AMR_mech_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$AMR_mech_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$AMR_mech_Shannon ~ pData(AMR_class_day11)$Treatment)

## Test over Time, treated
wilcox.test(pData(AMR_class_arrival_treated)$AMR_class_Richness, pData(AMR_class_day11_treated)$AMR_class_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$AMR_class_Shannon, pData(AMR_class_day11_treated)$AMR_class_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$AMR_mech_Richness, pData(AMR_class_day11_treated)$AMR_mech_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$AMR_mech_Shannon, pData(AMR_class_day11_treated)$AMR_mech_Shannon , paired = TRUE) ## Sig

# Test over time, Untreated
wilcox.test(pData(AMR_class_arrival_untreated)$AMR_class_Richness, pData(AMR_class_day11_untreated)$AMR_class_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$AMR_class_Shannon, pData(AMR_class_day11_untreated)$AMR_class_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$AMR_mech_Richness, pData(AMR_class_day11_untreated)$AMR_mech_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$AMR_mech_Shannon, pData(AMR_class_day11_untreated)$AMR_mech_Shannon , paired = TRUE)

#### Kraken diversity comparisons 
###
##
#
## Test by treatment, Arrival
wilcox.test(pData(AMR_class_arrival)$kraken_phylum_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$kraken_phylum_Shannon ~ pData(AMR_class_arrival)$Treatment)

wilcox.test(pData(AMR_class_arrival)$kraken_class_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$kraken_class_Shannon ~ pData(AMR_class_arrival)$Treatment)

wilcox.test(pData(AMR_class_arrival)$kraken_order_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$kraken_order_Shannon ~ pData(AMR_class_arrival)$Treatment)

wilcox.test(pData(AMR_class_arrival)$kraken_family_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$kraken_family_Shannon ~ pData(AMR_class_arrival)$Treatment)

wilcox.test(pData(AMR_class_arrival)$kraken_genus_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$kraken_genus_Shannon ~ pData(AMR_class_arrival)$Treatment)

wilcox.test(pData(AMR_class_arrival)$kraken_species_Richness ~ pData(AMR_class_arrival)$Treatment)
wilcox.test(pData(AMR_class_arrival)$kraken_species_Shannon ~ pData(AMR_class_arrival)$Treatment)

## Test by treatment, day11
wilcox.test(pData(AMR_class_day11)$kraken_phylum_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$kraken_phylum_Shannon ~ pData(AMR_class_day11)$Treatment)

wilcox.test(pData(AMR_class_day11)$kraken_class_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$kraken_class_Shannon ~ pData(AMR_class_day11)$Treatment)

wilcox.test(pData(AMR_class_day11)$kraken_order_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$kraken_order_Shannon ~ pData(AMR_class_day11)$Treatment)

wilcox.test(pData(AMR_class_day11)$kraken_family_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$kraken_family_Shannon ~ pData(AMR_class_day11)$Treatment)

wilcox.test(pData(AMR_class_day11)$kraken_genus_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$kraken_genus_Shannon ~ pData(AMR_class_day11)$Treatment)

wilcox.test(pData(AMR_class_day11)$kraken_species_Richness ~ pData(AMR_class_day11)$Treatment)
wilcox.test(pData(AMR_class_day11)$kraken_species_Shannon ~ pData(AMR_class_day11)$Treatment)




## Test over Time, treated
wilcox.test(pData(AMR_class_arrival_treated)$kraken_domain_Richness, pData(AMR_class_day11_treated)$kraken_domain_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_domain_Shannon, pData(AMR_class_day11_treated)$kraken_domain_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$kraken_phylum_Richness, pData(AMR_class_day11_treated)$kraken_phylum_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_phylum_Shannon, pData(AMR_class_day11_treated)$kraken_phylum_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$kraken_class_Richness, pData(AMR_class_day11_treated)$kraken_class_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_class_Shannon, pData(AMR_class_day11_treated)$kraken_class_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$kraken_order_Richness, pData(AMR_class_day11_treated)$kraken_order_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_order_Shannon, pData(AMR_class_day11_treated)$kraken_order_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$kraken_family_Richness, pData(AMR_class_day11_treated)$kraken_family_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_family_Shannon, pData(AMR_class_day11_treated)$kraken_family_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$kraken_genus_Richness, pData(AMR_class_day11_treated)$kraken_genus_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_genus_Shannon, pData(AMR_class_day11_treated)$kraken_genus_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_treated)$kraken_species_Richness, pData(AMR_class_day11_treated)$kraken_species_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_treated)$kraken_species_Shannon, pData(AMR_class_day11_treated)$kraken_species_Shannon , paired = TRUE)


## Test over Time, untreated
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_domain_Richness, pData(AMR_class_day11_untreated)$kraken_domain_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_domain_Shannon, pData(AMR_class_day11_untreated)$kraken_domain_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$kraken_phylum_Richness, pData(AMR_class_day11_untreated)$kraken_phylum_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_phylum_Shannon, pData(AMR_class_day11_untreated)$kraken_phylum_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$kraken_class_Richness, pData(AMR_class_day11_untreated)$kraken_class_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_class_Shannon, pData(AMR_class_day11_untreated)$kraken_class_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$kraken_order_Richness, pData(AMR_class_day11_untreated)$kraken_order_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_order_Shannon, pData(AMR_class_day11_untreated)$kraken_order_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$kraken_family_Richness, pData(AMR_class_day11_untreated)$kraken_family_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_family_Shannon, pData(AMR_class_day11_untreated)$kraken_family_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$kraken_genus_Richness, pData(AMR_class_day11_untreated)$kraken_genus_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_genus_Shannon, pData(AMR_class_day11_untreated)$kraken_genus_Shannon , paired = TRUE)

wilcox.test(pData(AMR_class_arrival_untreated)$kraken_species_Richness, pData(AMR_class_day11_untreated)$kraken_species_Richness , paired = TRUE)
wilcox.test(pData(AMR_class_arrival_untreated)$kraken_species_Shannon, pData(AMR_class_day11_untreated)$kraken_species_Shannon , paired = TRUE)


