# MEG_R_metagenomic_analysis
Code to combine microbiome (kraken2 or Qiime2) and resistome results (AmrPlusPlus) for statistical analysis with metagenomeSeq.

# How to use this code:
* Start with the "staging_script.R"
* Once this is filled out, it calls scripts from the "scripts" directory.


## Statistical analysis
---

* Statistical analysis is undoubtedly the most complex component required for characterizing metagenomic sequencing samples. We adopt a lot of the techniques developed for ecology to analyze the multivariate data representing the microbiome and resistome.

* We use the R programming software with various R software packages and you can see the main code we use at this repository (https://github.com/EnriqueDoster/MEG_R_metagenomic_analysis).
  * Diversity metrics and ordination : vegan
  * Count normalization and ZIG model: metagenomeSeq
  * data handling and manipulation: data.table
  * plotting : ggplot2
  * handling Qiime2 results and unifrac distances: phyloseq

### Count normalization

* The count of classified reads will first need to normalized to account for differences in sequencing performance between samples. There are many ways to normalize counts and we use “Cumulative sum scaling”(CSS) developed by Paulson et al 2013 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010126/).
* Then, the CSS-normalized counts are used for all the following methods of summarizing metagenomic samples (and check out this helpful website https://mb3is.megx.net/gustame/home):

### Alpha diversity

* Alpha diversity
  * Indices like “Richness”, “Shannon’s index” or "Simpson's index" to summarize the microbiome/resistome with a single value that represents the unique number of features and how evenly distributed the counts were, respectively.
* Plotting alpha diversity values with boxplots.


### Differential feature abundance

* Feature abundance in metagenomic sequencing projects
  * Comparing relative abundances is troublesome because of issues inherent to compositional data (ie. the proportion of a feature is directly affected by changes in proportion of other features) and because count tables contain many zeroes (sparse).
      * Therefore, we use a “Zero-inflated Gaussian model” (ZIG) that allows us to combine two different distributions and fit a model that better represents metagenomic count data. Log fold changes in abundance are estimated and the Benjamini-Hochberg method is employed to adjust p-values for multiple testing.
        * Take a look at the “stats” directory and look for the file named “Group_Microbiome_Phylum_GroupNorm-GroupProb_Model_Contrasts.csv”.
         * These results have to be interpreted carefully because low abundance features will often be significantly different between sample groups. In combination with your knowledge of which features are most abundant in your samples, pick a threshold value for “Avg. Expression” and report everything above that threshold. For example, you might say something like out of 27 compared, X phyla had average expressions > 1 and X were significantly different between sample groups (p-value < 0.05).

### Ordination

* Ordination
  * Non-metric multidimensional scaling


