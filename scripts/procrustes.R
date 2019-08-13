####### Procrustes analysis ######

AMR_analytic_data <- meg_filter_data(AMR_raw_analytic_data, filter_min_threshold = 0.15)

## AMR subsets
AMR_Nursery_samples = which(pData(AMR_analytic_data[[1]])$Housing_facility == "Nursery")
AMR_Farrowing_weaning_samples = which(pData(AMR_analytic_data[[1]])$Housing_facility == "Farrowing_weaning")
AMR_Finisher_samples = which(pData(AMR_analytic_data[[1]])$Housing_facility == "Finisher")

#Class
AMR_class_Nursery <- AMR_analytic_data[[1]][, AMR_Nursery_samples]
AMR_class_Farrowing_weaning <- AMR_analytic_data[[1]][, AMR_Farrowing_weaning_samples]
AMR_class_Finisher <- AMR_analytic_data[[1]][, AMR_Finisher_samples]

#Mechanism
AMR_mech_Nursery <- AMR_analytic_data[[2]][, AMR_Nursery_samples]
AMR_mech_Farrowing_weaning <- AMR_analytic_data[[2]][, AMR_Farrowing_weaning_samples]
AMR_mech_Finisher <- AMR_analytic_data[[2]][, AMR_Finisher_samples]


## Microbiome subsets
microbiome_Nursery_samples = which(pData(microbiome_analytic_data[[1]])$Housing_facility == "Nursery")
microbiome_Farrowing_weaning_samples = which(pData(microbiome_analytic_data[[1]])$Housing_facility == "Farrowing_weaning")
microbiome_Finisher_samples = which(pData(microbiome_analytic_data[[1]])$Housing_facility == "Finisher")

#phylum
microbiome_phylum_Nursery <- microbiome_analytic_data[[2]][,microbiome_Nursery_samples]
microbiome_phylum_Farrowing_weaning <- microbiome_analytic_data[[2]][, microbiome_Farrowing_weaning_samples]
microbiome_phylum_Finisher <- microbiome_analytic_data[[2]][, microbiome_Finisher_samples]

#class
microbiome_class_Nursery <- microbiome_analytic_data[[3]][,microbiome_Nursery_samples]
microbiome_class_Farrowing_weaning <- microbiome_analytic_data[[3]][, microbiome_Farrowing_weaning_samples]
microbiome_class_Finisher <- microbiome_analytic_data[[3]][, microbiome_Finisher_samples]

#order
microbiome_order_Nursery <- microbiome_analytic_data[[4]][,microbiome_Nursery_samples]
microbiome_order_Farrowing_weaning <- microbiome_analytic_data[[4]][, microbiome_Farrowing_weaning_samples]
microbiome_order_Finisher <- microbiome_analytic_data[[4]][, microbiome_Finisher_samples]

#family
microbiome_family_Nursery <- microbiome_analytic_data[[5]][,microbiome_Nursery_samples]
microbiome_family_Farrowing_weaning <- microbiome_analytic_data[[5]][, microbiome_Farrowing_weaning_samples]
microbiome_family_Finisher <- microbiome_analytic_data[[5]][, microbiome_Finisher_samples]

#genus
microbiome_genus_Nursery <- microbiome_analytic_data[[6]][,microbiome_Nursery_samples]
microbiome_genus_Farrowing_weaning <- microbiome_analytic_data[[6]][, microbiome_Farrowing_weaning_samples]
microbiome_genus_Finisher <- microbiome_analytic_data[[6]][, microbiome_Finisher_samples]

#species
microbiome_species_Nursery <- microbiome_analytic_data[[7]][,microbiome_Nursery_samples]
microbiome_species_Farrowing_weaning <- microbiome_analytic_data[[7]][, microbiome_Farrowing_weaning_samples]
microbiome_species_Finisher <- microbiome_analytic_data[[7]][, microbiome_Finisher_samples]






#### PROCRUSTES #####
### Compare Resistomes and microbiome at each housing facility
dist_AMR_class_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(AMR_class_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_AMR_class_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(AMR_class_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_AMR_class_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(AMR_class_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_AMR_mech_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(AMR_mech_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_AMR_mech_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(AMR_mech_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_AMR_mech_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(AMR_mech_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)


## Microbiome
dist_microbiome_phylum_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_class_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_class_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_class_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_order_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_order_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_order_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_family_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_family_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_family_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_genus_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_genus_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_genus_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_species_Nursery <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_Nursery, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_species_Farrowing_weaning <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_Farrowing_weaning, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_species_Finisher <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_Finisher, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)




# Protest AMR vs microbiome phylum
proTest_AMR_class_Nursery_vs_phylum_Nursery <- protest(dist_AMR_class_Nursery,dist_microbiome_phylum_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Farrowing_weaning_vs_phylum_Farrowing_weaning <- protest(dist_AMR_class_Farrowing_weaning,dist_microbiome_phylum_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Finisher_vs_phylum_Finisher <- protest(dist_AMR_class_Finisher,dist_microbiome_phylum_Finisher , scores= "sites", symmetric = TRUE)

proTest_AMR_mech_Nursery_vs_phylum_Nursery <- protest(dist_AMR_mech_Nursery,dist_microbiome_phylum_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Farrowing_weaning_vs_phylum_Farrowing_weaning <- protest(dist_AMR_mech_Farrowing_weaning,dist_microbiome_phylum_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Finisher_vs_phylum_Finisher <- protest(dist_AMR_mech_Finisher,dist_microbiome_phylum_Finisher , scores= "sites", symmetric = TRUE)

# Protest AMR vs microbiome class
proTest_AMR_class_Nursery_vs_class_Nursery <- protest(dist_AMR_class_Nursery,dist_microbiome_class_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Farrowing_weaning_vs_class_Farrowing_weaning <- protest(dist_AMR_class_Farrowing_weaning,dist_microbiome_class_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Finisher_vs_class_Finisher <- protest(dist_AMR_class_Finisher,dist_microbiome_class_Finisher , scores= "sites", symmetric = TRUE)

proTest_AMR_mech_Nursery_vs_class_Nursery <- protest(dist_AMR_mech_Nursery,dist_microbiome_class_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Farrowing_weaning_vs_class_Farrowing_weaning <- protest(dist_AMR_mech_Farrowing_weaning,dist_microbiome_class_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Finisher_vs_class_Finisher <- protest(dist_AMR_mech_Finisher,dist_microbiome_class_Finisher , scores= "sites", symmetric = TRUE)

# Protest AMR vs microbiome order
proTest_AMR_class_Nursery_vs_order_Nursery <- protest(dist_AMR_class_Nursery,dist_microbiome_order_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Farrowing_weaning_vs_order_Farrowing_weaning <- protest(dist_AMR_class_Farrowing_weaning,dist_microbiome_order_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Finisher_vs_order_Finisher <- protest(dist_AMR_class_Finisher,dist_microbiome_order_Finisher , scores= "sites", symmetric = TRUE)

proTest_AMR_mech_Nursery_vs_order_Nursery <- protest(dist_AMR_mech_Nursery,dist_microbiome_order_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Farrowing_weaning_vs_order_Farrowing_weaning <- protest(dist_AMR_mech_Farrowing_weaning,dist_microbiome_order_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Finisher_vs_order_Finisher <- protest(dist_AMR_mech_Finisher,dist_microbiome_order_Finisher , scores= "sites", symmetric = TRUE)

# Protest AMR vs microbiome phylum
proTest_AMR_class_Nursery_vs_phylum_Nursery <- protest(dist_AMR_class_Nursery,dist_microbiome_phylum_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Farrowing_weaning_vs_phylum_Farrowing_weaning <- protest(dist_AMR_class_Farrowing_weaning,dist_microbiome_phylum_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Finisher_vs_phylum_Finisher <- protest(dist_AMR_class_Finisher,dist_microbiome_phylum_Finisher , scores= "sites", symmetric = TRUE)

proTest_AMR_mech_Nursery_vs_phylum_Nursery <- protest(dist_AMR_mech_Nursery,dist_microbiome_phylum_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Farrowing_weaning_vs_phylum_Farrowing_weaning <- protest(dist_AMR_mech_Farrowing_weaning,dist_microbiome_phylum_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Finisher_vs_phylum_Finisher <- protest(dist_AMR_mech_Finisher,dist_microbiome_phylum_Finisher , scores= "sites", symmetric = TRUE)

# Protest AMR vs microbiome phylum
proTest_AMR_class_Nursery_vs_phylum_Nursery <- protest(dist_AMR_class_Nursery,dist_microbiome_phylum_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Farrowing_weaning_vs_phylum_Farrowing_weaning <- protest(dist_AMR_class_Farrowing_weaning,dist_microbiome_phylum_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_class_Finisher_vs_phylum_Finisher <- protest(dist_AMR_class_Finisher,dist_microbiome_phylum_Finisher , scores= "sites", symmetric = TRUE)

proTest_AMR_mech_Nursery_vs_phylum_Nursery <- protest(dist_AMR_mech_Nursery,dist_microbiome_phylum_Nursery , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Farrowing_weaning_vs_phylum_Farrowing_weaning <- protest(dist_AMR_mech_Farrowing_weaning,dist_microbiome_phylum_Farrowing_weaning , scores= "sites", symmetric = TRUE)
proTest_AMR_mech_Finisher_vs_phylum_Finisher <- protest(dist_AMR_mech_Finisher,dist_microbiome_phylum_Finisher , scores= "sites", symmetric = TRUE)





## Plotting procrustes
# Finisher, Phylum
plot(proTest_AMR_class_Finisher_vs_phylum_Finisher, main="Finisher facility - AMR class vs Microbiome phylum Procrustes")
text(0, 0.4, "M^2= 0.62, p-value = 0.23",cex = 1)
plot(proTest_AMR_mech_Finisher_vs_phylum_Finisher, main="Finisher facility - AMR mech vs Microbiome phylum Procrustes")
text(0, 0.4, "M^2= 0.57, p-value = 0.17",cex = 1)

# Finisher, class
plot(proTest_AMR_class_Finisher_vs_class_Finisher, main="Finisher facility - AMR class vs Microbiome class Procrustes")
text(0, 0.4, "M^2= 0.46, p-value = 0.09",cex = 1)
plot(proTest_AMR_mech_Finisher_vs_class_Finisher, main="Finisher facility - AMR mech vs Microbiome class Procrustes")
text(0, 0.4, "M^2= 0.57, p-value = 0.23",cex = 1)

# Finisher, order
plot(proTest_AMR_class_Finisher_vs_order_Finisher, main="Finisher facility - AMR class vs Microbiome order Procrustes")
text(0, 0.3, "M^2= 0.48, p-value = 0.10",cex = 1.5)
plot(proTest_AMR_mech_Finisher_vs_order_Finisher, main="Finisher facility - AMR mech vs Microbiome order Procrustes")
text(0, 0.3, "M^2= 0.59, p-value = 0.24",cex = 1.5)

# Nursery, order
plot(proTest_AMR_class_Nursery_vs_order_Nursery, main="Nursery facility - AMR class vs Microbiome order Procrustes")
text(0, 0.3, "M^2= 0.93, p-value = 0.49",cex = 1.5)
plot(proTest_AMR_mech_Nursery_vs_order_Nursery, main="Nursery facility - AMR mech vs Microbiome order Procrustes")
text(0, 0.25, "M^2= 0.77, p-value = 0.28",cex = 1.5)

# Farrowing_weaning, order
plot(proTest_AMR_class_Farrowing_weaning_vs_order_Farrowing_weaning, main="Farrowing_weaning facility - AMR class vs Microbiome order Procrustes")
text(0, 0.3, "M^2= 0.77, p-value = 0.79",cex = 1.5)
plot(proTest_AMR_mech_Farrowing_weaning_vs_order_Farrowing_weaning, main="Farrowing_weaning facility - AMR mech vs Microbiome order Procrustes")
text(0, 0.3, "M^2= 0.69, p-value = 0.79",cex = 1.5)









dist_AMR_mech_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_AMR_class_arrival_comparison <- protest(dist_AMR_class_Nursery,dist_AMR_class_Farrowing_weaning, scores= "sites", symmetric = TRUE)
proTest_AMR_class_day11_comparison <- protest(dist_AMR_class_day11_treated,dist_AMR_class_day11_untreated, scores= "sites", symmetric = TRUE)
proTest_AMR_mech_arrival_comparison <- protest(dist_AMR_mech_arrival_treated,dist_AMR_mech_Farrowing_weaning, scores= "sites", symmetric = TRUE)
proTest_AMR_mech_day11_comparison <- protest(dist_AMR_mech_day11_treated,dist_AMR_mech_day11_untreated, scores= "sites", symmetric = TRUE)



## microbiome, at arrival then at day 11
dist_microbiome_phylum_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_day11_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_day11_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

dist_microbiome_phylum_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_microbiome_phylum_arrival_comparison <- protest(dist_microbiome_phylum_arrival_treated,dist_microbiome_phylum_arrival_untreated, scores= "sites", symmetric = TRUE)
proTest_microbiome_phylum_day11_comparison <- protest(dist_microbiome_phylum_day11_treated,dist_microbiome_phylum_day11_untreated, scores= "sites", symmetric = TRUE)

dist_microbiome_class_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_class_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_class_day11_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_day11_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_class_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_class_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_microbiome_class_arrival_comparison <- protest(dist_microbiome_class_arrival_treated,dist_microbiome_class_arrival_untreated, scores= "sites", symmetric = TRUE)
proTest_microbiome_class_day11_comparison <- protest(dist_microbiome_class_day11_treated,dist_microbiome_class_day11_untreated, scores= "sites", symmetric = TRUE)

dist_microbiome_order_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_order_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_order_day11_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_day11_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_order_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_order_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_microbiome_order_arrival_comparison <- protest(dist_microbiome_order_arrival_treated,dist_microbiome_order_arrival_untreated, scores= "sites", symmetric = TRUE)
proTest_microbiome_order_day11_comparison <- protest(dist_microbiome_order_day11_treated,dist_microbiome_order_day11_untreated, scores= "sites", symmetric = TRUE)


dist_microbiome_family_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_family_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_family_day11_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_day11_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_family_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_family_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_microbiome_family_arrival_comparison <- protest(dist_microbiome_family_arrival_treated,dist_microbiome_family_arrival_untreated, scores= "sites", symmetric = TRUE)
proTest_microbiome_family_day11_comparison <- protest(dist_microbiome_family_day11_treated,dist_microbiome_family_day11_untreated, scores= "sites", symmetric = TRUE)


dist_microbiome_genus_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_genus_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_genus_day11_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_day11_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_genus_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_genus_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_microbiome_genus_arrival_comparison <- protest(dist_microbiome_genus_arrival_treated,dist_microbiome_genus_arrival_untreated, scores= "sites", symmetric = TRUE)
proTest_microbiome_genus_day11_comparison <- protest(dist_microbiome_genus_day11_treated,dist_microbiome_genus_day11_untreated, scores= "sites", symmetric = TRUE)


dist_microbiome_species_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_species_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_species_day11_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_day11_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_species_day11_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_species_day11_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)

proTest_microbiome_species_arrival_comparison <- protest(dist_microbiome_species_arrival_treated,dist_microbiome_species_arrival_untreated, scores= "sites", symmetric = TRUE)
proTest_microbiome_species_day11_comparison <- protest(dist_microbiome_species_day11_treated,dist_microbiome_species_day11_untreated, scores= "sites", symmetric = TRUE)




### AMR class and mech vs microbiome phylum ##
dist_AMR_class_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(AMR_class_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
proTest_AMR_class_vs_phylum_arrival_treated <- protest(dist_AMR_class_arrival_treated,dist_microbiome_phylum_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(AMR_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
proTest_AMR_class_vs_phylum_arrival_untreated <- protest(dist_AMR_class_arrival_untreated,dist_microbiome_phylum_arrival_untreated, scores= "sites", symmetric = TRUE)

dist_AMR_mech_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(AMR_mech_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_arrival_treated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_arrival_treated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
proTest_AMR_mech_vs_phylum_arrival_treated <- protest(dist_AMR_mech_arrival_treated,dist_microbiome_phylum_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(AMR_mech_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
dist_microbiome_phylum_arrival_untreated <- metaMDS(vegdist(decostand(t(MRcounts(microbiome_phylum_arrival_untreated, norm=FALSE)), "hell"), "euclidean"),distance="none",symmetric=TRUE, trymax=100)
proTest_AMR_mech_vs_phylum_arrival_untreated <- protest(dist_AMR_mech_arrival_untreated,dist_microbiome_phylum_arrival_untreated, scores= "sites", symmetric = TRUE)

dist_AMR_class_day11_treated <- vegdist(decostand(t(MRcounts(AMR_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_phylum_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_phylum_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_phylum_day11_treated <- protest(dist_AMR_class_day11_treated,dist_microbiome_phylum_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_phylum_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_phylum_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_phylum_day11_untreated <- protest(dist_AMR_class_day11_untreated,dist_microbiome_phylum_day11_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_treated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_phylum_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_phylum_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_phylum_day11_treated <- protest(dist_AMR_mech_day11_treated,dist_microbiome_phylum_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_phylum_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_phylum_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_phylum_day11_untreated <- protest(dist_AMR_mech_day11_untreated,dist_microbiome_phylum_day11_untreated, scores= "sites", symmetric = TRUE)

### AMR class and mech vs microbiome class ##
dist_AMR_class_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_class_arrival_treated <- protest(dist_AMR_class_arrival_treated,dist_microbiome_class_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_class_arrival_untreated <- protest(dist_AMR_class_arrival_untreated,dist_microbiome_class_arrival_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_class_arrival_treated <- protest(dist_AMR_mech_arrival_treated,dist_microbiome_class_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_class_arrival_untreated <- protest(dist_AMR_mech_arrival_untreated,dist_microbiome_class_arrival_untreated, scores= "sites", symmetric = TRUE)

dist_AMR_class_day11_treated <- vegdist(decostand(t(MRcounts(AMR_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_class_day11_treated <- protest(dist_AMR_class_day11_treated,dist_microbiome_class_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_class_day11_untreated <- protest(dist_AMR_class_day11_untreated,dist_microbiome_class_day11_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_treated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_class_day11_treated <- protest(dist_AMR_mech_day11_treated,dist_microbiome_class_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_class_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_class_day11_untreated <- protest(dist_AMR_mech_day11_untreated,dist_microbiome_class_day11_untreated, scores= "sites", symmetric = TRUE)

### AMR class and mech vs microbiome order ##
dist_AMR_class_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_order_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_order_arrival_treated <- protest(dist_AMR_class_arrival_treated,dist_microbiome_order_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_order_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_order_arrival_untreated <- protest(dist_AMR_class_arrival_untreated,dist_microbiome_order_arrival_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_order_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_order_arrival_treated <- protest(dist_AMR_mech_arrival_treated,dist_microbiome_order_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_order_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_order_arrival_untreated <- protest(dist_AMR_mech_arrival_untreated,dist_microbiome_order_arrival_untreated, scores= "sites", symmetric = TRUE)

dist_AMR_class_day11_treated <- vegdist(decostand(t(MRcounts(AMR_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_order_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_order_day11_treated <- protest(dist_AMR_class_day11_treated,dist_microbiome_order_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_order_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_order_day11_untreated <- protest(dist_AMR_class_day11_untreated,dist_microbiome_order_day11_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_treated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_order_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_order_day11_treated <- protest(dist_AMR_mech_day11_treated,dist_microbiome_order_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_order_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_order_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_order_day11_untreated <- protest(dist_AMR_mech_day11_untreated,dist_microbiome_order_day11_untreated, scores= "sites", symmetric = TRUE)

### AMR class and mech vs microbiome family ##
dist_AMR_class_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_family_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_family_arrival_treated <- protest(dist_AMR_class_arrival_treated,dist_microbiome_family_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_family_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_family_arrival_untreated <- protest(dist_AMR_class_arrival_untreated,dist_microbiome_family_arrival_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_family_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_family_arrival_treated <- protest(dist_AMR_mech_arrival_treated,dist_microbiome_family_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_family_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_family_arrival_untreated <- protest(dist_AMR_mech_arrival_untreated,dist_microbiome_family_arrival_untreated, scores= "sites", symmetric = TRUE)







dist_AMR_class_day11_treated <- vegdist(decostand(t(MRcounts(AMR_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_family_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_family_day11_treated <- protest(dist_AMR_class_day11_treated,dist_microbiome_family_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_family_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_family_day11_untreated <- protest(dist_AMR_class_day11_untreated,dist_microbiome_family_day11_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_treated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_family_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_family_day11_treated <- protest(dist_AMR_mech_day11_treated,dist_microbiome_family_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_family_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_family_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_family_day11_untreated <- protest(dist_AMR_mech_day11_untreated,dist_microbiome_family_day11_untreated, scores= "sites", symmetric = TRUE)

### AMR class and mech vs microbiome genus ##
dist_AMR_class_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_genus_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_genus_arrival_treated <- protest(dist_AMR_class_arrival_treated,dist_microbiome_genus_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_genus_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_genus_arrival_untreated <- protest(dist_AMR_class_arrival_untreated,dist_microbiome_genus_arrival_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_genus_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_genus_arrival_treated <- protest(dist_AMR_mech_arrival_treated,dist_microbiome_genus_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_genus_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_genus_arrival_untreated <- protest(dist_AMR_mech_arrival_untreated,dist_microbiome_genus_arrival_untreated, scores= "sites", symmetric = TRUE)

dist_AMR_class_day11_treated <- vegdist(decostand(t(MRcounts(AMR_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_genus_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_genus_day11_treated <- protest(dist_AMR_class_day11_treated,dist_microbiome_genus_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_genus_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_genus_day11_untreated <- protest(dist_AMR_class_day11_untreated,dist_microbiome_genus_day11_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_treated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_genus_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_genus_day11_treated <- protest(dist_AMR_mech_day11_treated,dist_microbiome_genus_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_genus_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_genus_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_genus_day11_untreated <- protest(dist_AMR_mech_day11_untreated,dist_microbiome_genus_day11_untreated, scores= "sites", symmetric = TRUE)


### AMR class and mech vs SPECIES

dist_AMR_class_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_species_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_species_arrival_treated <- protest(dist_AMR_class_arrival_treated,dist_microbiome_species_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_class_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_species_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_species_arrival_untreated <- protest(dist_AMR_class_arrival_untreated,dist_microbiome_species_arrival_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_treated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_arrival_treated <- vegdist(decostand(t(MRcounts(microbiome_species_arrival_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_species_arrival_treated <- protest(dist_AMR_mech_arrival_treated,dist_microbiome_species_arrival_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_arrival_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_arrival_untreated <- vegdist(decostand(t(MRcounts(microbiome_species_arrival_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_species_arrival_untreated <- protest(dist_AMR_mech_arrival_untreated,dist_microbiome_species_arrival_untreated, scores= "sites", symmetric = TRUE)


dist_AMR_class_day11_treated <- vegdist(decostand(t(MRcounts(AMR_class_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_species_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_species_day11_treated <- protest(dist_AMR_class_day11_treated,dist_microbiome_species_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_class_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_class_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_species_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_class_vs_species_day11_untreated <- protest(dist_AMR_class_day11_untreated,dist_microbiome_species_day11_untreated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_treated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_treated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_day11_treated <- vegdist(decostand(t(MRcounts(microbiome_species_day11_treated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_species_day11_treated <- protest(dist_AMR_mech_day11_treated,dist_microbiome_species_day11_treated, scores= "sites", symmetric = TRUE)
dist_AMR_mech_day11_untreated <- vegdist(decostand(t(MRcounts(AMR_mech_day11_untreated, norm=FALSE)), "hell"), "euclidean")
dist_microbiome_species_day11_untreated <- vegdist(decostand(t(MRcounts(microbiome_species_day11_untreated, norm=FALSE)), "hell"), "euclidean")
proTest_AMR_mech_vs_species_day11_untreated <- protest(dist_AMR_mech_day11_untreated,dist_microbiome_species_day11_untreated, scores= "sites", symmetric = TRUE)


### Go over the results of procrustes ###
## Arrival comparisons of microbiome and resistome
proTest_AMR_class_arrival_comparison
proTest_AMR_mech_arrival_comparison

proTest_AMR_class_day11_comparison
proTest_AMR_mech_day11_comparison

proTest_microbiome_phylum_arrival_comparison
proTest_microbiome_class_arrival_comparison
proTest_microbiome_order_arrival_comparison
proTest_microbiome_family_arrival_comparison
proTest_microbiome_genus_arrival_comparison
proTest_microbiome_species_arrival_comparison

proTest_microbiome_phylum_day11_comparison
proTest_microbiome_class_day11_comparison
proTest_microbiome_order_day11_comparison
proTest_microbiome_family_day11_comparison
proTest_microbiome_genus_day11_comparison
proTest_microbiome_species_day11_comparison


## comparison of resistome vs microbiome
proTest_AMR_class_vs_phylum_arrival_treated
proTest_AMR_mech_vs_phylum_arrival_treated
proTest_AMR_class_vs_phylum_arrival_untreated
proTest_AMR_mech_vs_phylum_arrival_untreated

proTest_AMR_class_vs_class_arrival_treated
proTest_AMR_mech_vs_class_arrival_treated
proTest_AMR_class_vs_class_arrival_untreated
proTest_AMR_mech_vs_class_arrival_untreated

proTest_AMR_class_vs_order_arrival_treated
proTest_AMR_mech_vs_order_arrival_treated
proTest_AMR_class_vs_order_arrival_untreated
proTest_AMR_mech_vs_order_arrival_untreated

proTest_AMR_class_vs_family_arrival_treated
proTest_AMR_mech_vs_family_arrival_treated
proTest_AMR_class_vs_family_arrival_untreated
proTest_AMR_mech_vs_family_arrival_untreated

proTest_AMR_class_vs_genus_arrival_treated
proTest_AMR_mech_vs_genus_arrival_treated
proTest_AMR_class_vs_genus_arrival_untreated
proTest_AMR_mech_vs_genus_arrival_untreated

proTest_AMR_class_vs_species_arrival_treated
proTest_AMR_mech_vs_species_arrival_treated
proTest_AMR_class_vs_species_arrival_untreated
proTest_AMR_mech_vs_species_arrival_untreated


## Day11
proTest_AMR_class_vs_phylum_day11_treated
proTest_AMR_mech_vs_phylum_day11_treated
proTest_AMR_class_vs_phylum_day11_untreated
proTest_AMR_mech_vs_phylum_day11_untreated

proTest_AMR_class_vs_class_day11_treated
proTest_AMR_mech_vs_class_day11_treated
proTest_AMR_class_vs_class_day11_untreated
proTest_AMR_mech_vs_class_day11_untreated

proTest_AMR_class_vs_order_day11_treated
proTest_AMR_mech_vs_order_day11_treated
proTest_AMR_class_vs_order_day11_untreated
proTest_AMR_mech_vs_order_day11_untreated

proTest_AMR_class_vs_family_day11_treated
proTest_AMR_mech_vs_family_day11_treated
proTest_AMR_class_vs_family_day11_untreated
proTest_AMR_mech_vs_family_day11_untreated

proTest_AMR_class_vs_genus_day11_treated
proTest_AMR_mech_vs_genus_day11_treated
proTest_AMR_class_vs_genus_day11_untreated
proTest_AMR_mech_vs_genus_day11_untreated

proTest_AMR_class_vs_species_day11_treated
proTest_AMR_mech_vs_species_day11_treated
proTest_AMR_class_vs_species_day11_untreated
proTest_AMR_mech_vs_species_day11_untreated
