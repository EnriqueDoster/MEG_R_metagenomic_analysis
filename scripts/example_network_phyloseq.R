#phyloseq network modeling attempt
#source("https://bioconductor.org/biocLite.R")
#biocLite("phyloseq")

library("phyloseq")
library("ggplot2")
library("biomformat")

# Count information
taxaCounts = as.matrix(MRcounts(AMR_analytic_data[[1]]))
OTU = otu_table(taxaCounts , taxa_are_rows = TRUE)

# Taxonomic information
taxaAMR = as.data.table(fData(AMR_analytic_data[[1]]))
# colnames(taxaAMR) <- "header"
# setkey(taxaAMR, header)
# setkey(annotations,header)
# taxaAMR <- annotations[taxaAMR]
TAX = tax_table(as.matrix(taxaAMR))
# taxa_names(TAX) <- taxaAMR$header
taxa_names(TAX) <- taxaAMR$Feature


# metadata
sampledata = sample_data(metadata)
sample_names(sampledata) <- metadata$ID
#make phyloseq object 
AMR_ps <- phyloseq(OTU,TAX,sampledata)

ntaxa()
nsamples()
sample_sums(Bushman)[1:10]
# rank_names(AMR_ps)
# sample_variables(AMR_ps)
# get_variable(Bushman, sample_variables(Bushman)[5])[1:10]
# get_taxa(Bushman, sample_names(Bushman)[5])[1:10]
# 
# taxa_names(AMR_ps)
# colnames(tax_table(Bushman)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")

#Note how a sample name is required by get_taxa, and vice versa. This might seem confusing at first, but get_taxa is returning 
#all the OTU abundances from one sample, while get_sample is returning the abundances from all samples for one OTU.
sample_sums(Bushman)[1:10]
taxa_sums(Bushman)[1:10]
get_taxa(Bushman, sample_names(Bushman)[5])[1:10]
get_sample(Bushman, taxa_names(Bushman)[5])[1:10]

#can subset
sub_physeq = subset_taxa(AMR_ps, class == "Tetracyclines")


#identical(sum(ent.logi), ntaxa(ent.trim))
plot_bar(AMR_ps, x="ID",fill = "class") #facet_grid=~Seq_time

library("ape")
### building trees ##
random_tree = rtree(ntaxa(AMR_ps), rooted=TRUE, tip.label=taxa_names(AMR_ps))
plot(random_tree)
AMR_ps=merge_phyloseq(AMR_ps,random_tree)

# Some plots
plot_net(AMR_ps , maxdist=0.4, color="DOF_categories", shape="Seq_time")
plot_tree(AMR_ps, color="DOF_categories", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(AMR_ps, color="Group", shape="Time", label.tips="taxa_names", ladderize="right", plot.margin=0.3)

ig = make_network(AMR_ps, type = "samples", distance = "bray", max.dist = 0.2)
plot_network(ig, AMR_ps, color = "DOF_categories", shape = "DDD_cat", line_weight = 0.4, label = NULL)
plot_network(ig, AMR_ps, color = "DDD_cat", shape = "DOF_categories", line_weight = 0.4, label = NULL,layout.method=layout.bipartite)

# layout.fruchterman.reingold
# layout.reingold.tilford
# layout.bipartite

## Top 10 features
AMR_ps_100 = prune_taxa(names(sort(taxa_sums(AMR_ps), TRUE))[1:6], AMR_ps)
jg = make_network(AMR_ps_100, "taxa", "jaccard", 0.3)
plot_network(jg, AMR_ps_100, "taxa",line_weight = 0.4, label = NULL)



data(enterotype)





### PreProcessing #########
## rarefy to even depth
eso = rarefy_even_depth(esophagus)
# plot(as(otu_table(eso), 'vector'), as(otu_table(esophagus), 'vector'))

data(GlobalPatterns)
GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
# remove the samples that have less than 20 total reads from Chlamydiae
GP.chl = prune_samples(names(which(sample_sums(GP.chl) >= 20)), GP.chl)
# (p = plot_tree(GP.chl, color='SampleType', shape='Family',
# label.tips='Genus', size='abundance'))
GP.chl.r = rarefy_even_depth(GP.chl)
# plot(as(otu_table(GP.chl.r), 'vector'), as(otu_table(GP.chl), 'vector'))
plot_ordination(GP.chl, ordinate(GP.chl, "MDS"), color = "SampleType") + geom_point(size = 5)

GPr  = transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

GP.chl = subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
GP.chl = prune_samples(sample_sums(GP.chl)>=20, GP.chl)

gpsfbg = tax_glom(gpsfb, "Family")
plot_tree(gpsfbg, color="SampleType", shape="Class", size="abundance")
transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU) )

#remove sparse
GP = filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Add new variable that matches list of factors
sample_data(GP)$human = factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue") )

#Standardize abundances to the median sequencing depth
total = median(sample_sums(GP))
standf = function(x, t=total) round(t * (x / sum(x)))
gps = transform_sample_counts(GP, standf)
#Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)


# differences in prune_taxa() vs. subset_taxa()
#These are two very different methods for subsetting OTUs in a dataset.
#That was prune_taxa, applicable when we know (or can provide) the OTU IDs of the OTUs we want to retain in the dataset,
#or a logical vector with the same length as ntaxa reseulting from a test (useful for genefilter orgenefilter_sample results (see next section)).
#Alternatively, you can subset based on taxonomy expressions, using subset_taxa.

topN = 20
most_abundant_taxa = sort(taxa_sums(GP), TRUE)[1:topN]
print(most_abundant_taxa)
GP20 = prune_taxa(names(most_abundant_taxa), GP)
# Exploratory tree #1 (Note, too ma) plot_tree(ex2, color='SampleType',
# label.tips='Family', size='abundance')
ntaxa(GP20)
length(get_taxa_unique(GP20, "Phylum"))

#####   Filtering low-variance OTUs
#Suppose we wanted to use the variance of OTUs across samples as a condition for filtering. For example, to remove OTUs that do not change much across all (or most) samples. Note that the value of the variance is highly-dependent on the sequencing effort of each sample (the total number of reads sequenced from a particular sample). Thus, we must first transform the sample counts to relative abundance, which is shown in more detail in the next section. The following code will create a version of the GP dataset in which the abundance values have been transformed to relative abundance within each sample, and then OTUs have been filtered to keep only those with variance greater than 0.00001 (assuming we wanted to pick an arbitrary threshold in this way).

GPr = transform_sample_counts(GP, function(x) x/sum(x))
GPf = filter_taxa(GPr, function(x) var(x) > 1e-05, TRUE)
#Get the names of the most-abundant phyla, and use for subsetting.
top.phyla = sort(tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum), TRUE)
top.phyla = top.phyla[1:5]
# Prune to just the most-abundant 5 phyla
GP2 = subset_taxa(GP, Phylum %in% names(top.phyla))
get_taxa_unique(GP2, "Phylum")
# Prune further, to top 200 most abundant taxa of top 5 phyla
GP2 = prune_taxa(names(sort(taxa_sums(GP2), TRUE)[1:200]), GP2)


### Network #########
ig = make_network(GP, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, GP, color = "SampleType", shape = "human", line_weight = 0.4, 
             label = NULL)
data(enterotype)
ig = make_network(enterotype, max.dist = 0.3)
plot_network(ig, enterotype, color = "SeqTech", shape = "Enterotype", line_weight = 0.4, 
             label = NULL)
data(GlobalPatterns)
# prune to just the top 100 most abundant OTUs across all samples (crude).
GP100 = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE))[1:100], GlobalPatterns)
jg = make_network(GP100, "taxa", "jaccard", 0.3)
plot_network(jg, GP100, "taxa", color = "Phylum", line_weight = 0.4, label = NULL)


enterotype = subset_samples(enterotype, !is.na(Enterotype))
plot_net(enterotype, maxdist = 0.4, point_label = "Sample_ID")
plot_net(enterotype, maxdist = 0.3, color = "SeqTech", shape="Enterotype")
ig <- make_network(enterotype, max.dist=0.3)
plot_network(ig, enterotype)
plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.4, label=NULL)

ig <- make_network(enterotype, max.dist=0.2)
plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.4, label=NULL)

ig <- make_network(enterotype, dist.fun="bray", max.dist=0.3)
plot_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.4, label=NULL)


######## Ordination ##########
#Here are some other examples. There are some 45 or so methods.
distance(esophagus)  # unweighted UniFrac
distance(esophagus, weighted = TRUE)  # weighted UniFrac
distance(esophagus, "jaccard")  # vegdist jaccard
distance(esophagus, "bray")  # vegdist bray-curtis
distance(esophagus, "gower")  # vegdist option 'gower'
distance(esophagus, "g")  # designdist method option 'g'
distance(esophagus, "minkowski")  # invokes a method from the base dist() function.
distance(esophagus, "(A+B-2*J)/(A+B)")  # designdist custom distance
distance("help")
distance("list")

GP.MDS = ordinate(GP100, method = "MDS", distance = "unifrac")
Here are just a few examples of other supported combinations.

GP.NMDS = ordinate(GP, "NMDS", "gower")
GP.NMDS = ordinate(GP, "NMDS", "bray")  # perform NMDS on bray-curtis distance
GP.NMDS.UF.ord = ordinate(GP, "NMDS")  # UniFrac. Takes a while.
GP.NMDS.wUF.ord = ordinate(GP, "NMDS", "unifrac", weighted = TRUE)  # weighted-UniFrac
GP.NMDS.gower = ordinate(GP, "NMDS", "gower")
require("ggplot2")
ptitle = "GP PCoA of UniFrac distance, GP most abundant 100 OTUs only"
p = plot_ordination(GP100, GP.MDS, type = "samples", color = "SampleType", title = ptitle)
p + geom_point(size = 5) + geom_polygon(aes(fill = SampleType))

## CCA ####
p2 = plot_ordination(GP2, ordinate(GP2, "CCA"), type = "samples", color = "SampleType")
p2 + geom_point(size = 5) + geom_polygon(aes(fill = SampleType))
plot_ordination(GP2, ordinate(GP2, "CCA"), type = "taxa", color = "Phylum") + 
  geom_point(size = 4)
plot_ordination(GP2, ordinate(GP2, "CCA"), type = "split")
plot_ordination(GP2, ordinate(GP2, "CCA"), type = "split", color = "SampleType")
plot_ordination(GP2, ordinate(GP2, "CCA"), type = "biplot", shape = "Phylum")
plot_ordination(GP2, ordinate(GP2, "CCA"), type = "split", color = "Phylum", 
                label = "SampleType")
plot_ordination(GP2, ordinate(GP2, "CCA"), type = "split", color = "SampleType", 
                shape = "Phylum", label = "SampleType")


#Examples of ordination plots
#(Note: can't map continuous variables to shape with plot_ordination function)
p4title = "Bushman dataset, PCoA/MDS ordination on Bray-Curtis distance"
Bushman.ord = ordinate(Bushman, method = "MDS", distance = "bray")
plot_ordination(Bushman, Bushman.ord, "samples", color = "OMEGA3_FATTY_ACIDS_G_AVE", 
                title = p4title) + geom_point(size = 4)



GPUF <- UniFrac(AMR_ps)
GloPa.pcoa = ordinate(AMR_ps, method="PCoA", distance=GPUF)
plot_scree(GloPa.pcoa, "Scree plot for Global Patterns, UniFrac/PCoA")

(p12 <- plot_ordination(AMR_ps, GloPa.pcoa, "samples", color="Group") + 
   geom_point(size=5) + geom_path() + scale_colour_hue(guide = FALSE) )

(p13 <- plot_ordination(AMR_ps, GloPa.pcoa, "samples", axes=c(1, 3),
                        color="Group") + geom_line() + geom_point(size=5) )

GP.NMDS <- ordinate(AMR_ps, "NMDS", GPUF)
(p <- plot_ordination(AMR_ps, GP.NMDS, "samples", color="Group") +
   geom_line() + geom_point(size=5) )

#Exploratory heatmaps and other stuff ############

my.ord    <- ordinate(AMR_ps)
plot_ordination(AMR_ps, my.ord, color="Group")

plot_heatmap(AMR_ps, taxa.label="Phylum")
gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
(p <- plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family"))
p$scales$scales[[1]]$name <- "My X-Axis"
p$scales$scales[[2]]$name <- "My Y-Axis"
print(p)

### More bar plots

data("GlobalPatterns")
gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")

plot_bar(gp.ch, "Family", fill="Genus", facet_grid=~SampleType)
p = plot_bar(gp.ch, "Family", fill="Genus", facet_grid=~SampleType)
p + geom_point(aes(x=Family, y=Abundance), color="black", position="jitter", size=3)

data("enterotype")
TopNOTUs <- names(sort(taxa_sums(enterotype), TRUE)[1:10])
ent10   <- prune_species(TopNOTUs, enterotype)
p = plot_bar(ent10, "Genus", fill="Genus", facet_grid=SeqTech~Enterotype)
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

## Richness plots ####
plot_richness(Bushman, x = "SEX", color = "SOLUBLE_DIETARY_FIBER_G_AVE") + geom_boxplot()

(p = plot_richness(Bushman, x = "SEX"))
p + geom_boxplot(data = p$data, aes(x = SEX, y = value, color = NULL), alpha = 0.1)
### Heatmaps #####
plot_heatmap(GP2, "NMDS", "bray", "SampleType", "Family")
plot_heatmap(GP2, "NMDS", "bray", "SampleType", "Family", trans = log_trans(10))
plot_heatmap(GP2, "NMDS", "bray", "SampleType", "Family", trans = identity_trans())
plot_heatmap(GP2, "NMDS", "bray", "SampleType", "Family", trans = boxcox_trans(0.15))
### Multiple testing #####
#In this example we will perform testing on fractional abundances to remove effect of differences in total sequencing across samples for same taxa.
#We first filter the low-variance taxa, avoiding the noise and "penalty" we pay for testing taxa that don't vary at much across samples. 
#In practice, these OTU abundances probably need additional preprocessing prior to testing, and many good methods in microarray analysis and 
#differential expression sequencing could probably apply to this data (but are not directly implemented/supported in phyloseq, yet). 
#Otherwise, beware that the old motto "garbage-in, garbage-out" can definitely apply to your data if you are not careful.
GPr = transform_sample_counts(GP, function(x) x/sum(x))
GP3f = filter_taxa(GPr, function(x) var(x) > 1e-05, TRUE)
#We are going to use the multtest wrapper included in phyloseq, mt. Try ?mt for help on this function. The following code uses this wrapper to calculate the multiple-inference-adjusted P-values, using the "human" sample variable.
GP.fwer.table = mt(GP3f, "human")
subset(GP.fwer.table, adjp < alpha)
#What if we want FDR instead of FWER?The easiest thing to do is use the stats::p.adjust function.
mtm = mt(GP3f, "human")
#Add the FDR correction based on the unadjusted ("raw") P-values.
mtm$FDR <- p.adjust(c(mtm[, "rawp"]), method = c("fdr"))
subset(mtm, FDR < 0.05)


#### Convert phyloseq into vegan ######
library("vegan")
packageVersion("vegan")
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

keepvariables = which(sapply(sample_data(Bushman), is.numeric))
bushsd = data.frame(sample_data(Bushman))[keepvariables]

bioenv(veganotu(Bushman), bushsd)

bioenv(veganotu(Bushman) ~ DEPTH + AGE + TOTAL_FAT_G_AVE + INSOLUBLE_DIETARY_FIBER_G_AVE, 
       bushsd)

#Ordination Example on the Gap Statistic   ############
library("cluster")
# Load data
data(enterotype)
# ordination
ent.ca = ordinate(enterotype, method = "CCA", distance = NULL)
pam1 = function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
x = scores(ent.ca, display = "sites")
# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn = clusGap(x[, 1:2], FUN = pam1, K.max = 6, B = 50)
gskmn

#That's nice. Just in case it is useful, let's look at what the wrapper-function might look like.

gap_statistic_ordination = function(ord, FUNcluster, type = "sites", K.max = 6, 
                                    axes = c(1:2), B = 500, verbose = interactive(), ...) {
  require("cluster")
  # If 'pam1' was chosen, use this internally defined call to pam
  if (FUNcluster == "pam1") {
    FUNcluster = function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
  }
  # Use the scores function to get the ordination coordinates
  x = scores(ord, display = type)
  # If axes not explicitly defined (NULL), then use all of them
  if (is.null(axes)) {
    axes = 1:ncol(x)
  }
  # Finally, perform, and return, the gap statistic calculation using
  # cluster::clusGap
  clusGap(x[, axes], FUN = FUNcluster, K.max = K.max, B = B, verbose = verbose, 
          ...)
}

plot_clusgap = function(clusgap, title = "Gap Statistic calculation results") {
  require("ggplot2")
  gstab = data.frame(clusgap$Tab, k = 1:nrow(clusgap$Tab))
  p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p = p + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim))
  p = p + ggtitle(title)
  return(p)
}

gs = gap_statistic_ordination(ent.ca, "pam1", B = 50, verbose = FALSE)
print(gs, method = "Tibs2001SEmax")

plot_clusgap(gs)

plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("k = 2 is best ... but  k = 3  pretty close")



