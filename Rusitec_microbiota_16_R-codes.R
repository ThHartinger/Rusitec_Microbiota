
-
Title: "Rusitec_Microbiota_16S"
Author: "Thomas Hartinger"
Date: "2nd May 2019"
-

    
#Load packages
library(ggplot2)
library(ape)
library(vegan)
library(phyloseq)
library(microbiome)
library(mgcv)
library(class)
library(codetools)
library(picante)

#Biomfile, trefile, mapping file -> Create phyloseq-object
biomfile <- import_biom("Rusitecmicrobiota.biom")
biomfile_mapfile <- import_qiime_sample_data("Rusitecmicrobiota_mapping.csv.txt") 
biomfile_merged <- merge_phyloseq(biomfile,biomfile_mapfile)
biomfile_tree <- read.tree("Rusitecmicrobiota.tre")
biomfile_merged <-merge_phyloseq(biomfile_merged,biomfile_tree)

#Clean the phyloseq object
colnames(tax_table(biomfile_merged)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
tax_table(biomfile_merged)[tax_table(biomfile_merged)[, "Domain"] == "NA", "Domain"] <- "Unknown_Domain"
tax_table(biomfile_merged)[tax_table(biomfile_merged)[, "Phylum"] == "p__", "Phylum"] <- "Unknown_Phylum"
tax_table(biomfile_merged)[tax_table(biomfile_merged)[, "Class"] == "c__", "Class"] <- "Unknown_Class"
tax_table(biomfile_merged)[tax_table(biomfile_merged)[, "Order"] == "o__", "Order"] <- "Unknown_Order"
tax_table(biomfile_merged)[tax_table(biomfile_merged)[, "Family"] == "f__", "Family"] <- "Unknown_Family"
tax_table(biomfile_merged)[tax_table(biomfile_merged)[, "Genus"] == "g__", "Genus"] <- "Unknown_Genus"
tax_table(biomfile_merged)[, colnames(tax_table(biomfile_merged))] <- gsub(tax_table(biomfile_merged)[,colnames(tax_table(biomfile_merged))], pattern = "[a-z]__", replacement = "")

#Subset the samples and prune 
ps1_unpruned <- subset_samples(biomfile_merged, Exclude != "yes")
ps1 <- prune_taxa(taxa_sums(ps1_unpruned) > 0, ps1_unpruned)
ps1 <- subset_taxa(ps1, Class != "Chloroplast")
ps1 <- subset_taxa(ps1, Family != "Mitochondria")
ps1_relabund <- transform_sample_counts(ps1, function(x)x/sum(x))

#Only LAB (liquid-associated microbes) and SAB (solid-associated microbes) for the respective time points (B = early and C = late) as relative abundances
ps1_LAB <- subset_samples(ps1, Phase != "Solid")
ps1_LAB_relabund <- transform_sample_counts(ps1_LAB, function(x)x/sum(x))
ps1_B_LAB <- subset_samples(ps1_LAB, Timepoint == "Early")
ps1_B_LAB_relabund <- transform_sample_counts(ps1_B_LAB, function(x)x/sum(x))
ps1_C_LAB <- subset_samples(ps1_LAB, Timepoint == "Late")
ps1_C_LAB_relabund <- transform_sample_counts(ps1_C_LAB, function(x)x/sum(x))
ps1_SAB <- subset_samples(ps1, Phase != "Liquid")
ps1_SAB_relabund <- transform_sample_counts(ps1_SAB, function(x)x/sum(x))
ps1_B_SAB <- subset_samples(ps1_SAB, Timepoint == "Early")
ps1_B_SAB_relabund <- transform_sample_counts(ps1_B_SAB, function(x)x/sum(x))
ps1_C_SAB <- subset_samples(ps1_SAB, Timepoint == "Late")
ps1_C_SAB_relabund <- transform_sample_counts(ps1_C_SAB, function(x)x/sum(x))
ps1_B <- subset_samples(ps1, Timepoint == "Early adaption time point")
ps1_B_relabund <- transform_sample_counts(ps1_B, function(x)x/sum(x))
ps1_C <- subset_samples(ps1, Timepoint == "Late adaption time point")
ps1_C_relabund <- transform_sample_counts(ps1_C, function(x)x/sum(x))

#Beta-diversity
#unifrac weighted
plot_ordination(ps1_relabund, ordinate(ps1_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="Phase", shape = "Timepoint", title = "PCoA unifrac: Phase and Timepoint", label="Phase")
#LAB_B samples unifrac
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: LAB_B DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: LAB_B WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: LAB_B SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: LAB_B DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: LAB_B DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: LAB_B WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#LAB_C samples unifrac
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: LAB_C DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: LAB_C WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: LAB_C SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: LAB_C DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: LAB_C DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: LAB_C WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_B samples unifrac
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: SAB_B DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: SAB_B WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: SAB_B SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: SAB_B DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: SAB_B DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: SAB_B WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_C samples unifrac
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: SAB_C DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: SAB_C WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: SAB_C SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: SAB_C DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: SAB_C DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: SAB_C WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)

#unifrac unweighted
plot_ordination(ps1_relabund, ordinate(ps1_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="Phase", shape = "Timepoint", title = "PCoA unifrac: Phase and Timepoint", label="Phase")
#LAB_B samples uunifrac
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: LAB_B DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: LAB_B WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: LAB_B SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: LAB_B DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: LAB_B DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: LAB_B WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#LAB_C samples uunifrac
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: LAB_C DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: LAB_C WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: LAB_C SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: LAB_C DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: LAB_C DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: LAB_C WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_B samples unifrac
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: SAB_B DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: SAB_B WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: SAB_B SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: SAB_B DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: SAB_B DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: SAB_B WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_C samples uunifrac
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", title = "PCoA unifrac: SAB_C DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", title = "PCoA unifrac: SAB_C WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", title = "PCoA unifrac: SAB_C SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", title = "PCoA unifrac: SAB_C DMxSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", title = "PCoA unifrac: SAB_C DMxWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", title = "PCoA unifrac: SAB_C WIxSA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)

#Test beta-diversity for significance by adonis for weighted unifrac
#Insert the respective phyloseq object (e.g. ps1_LAB) and the variable of interest
rel_distance <- distance(ps1_B_LAB, "unifrac")
adonis(rel_distance ~ DM, as(sample_data(ps1_B_LAB), "data.frame"))

#Test beta-diversity for significance by adonis for unweighted unifrac
#Insert the respective phyloseq object (e.g. ps1_LAB) and the variable of interest
rel_distance <- distance(ps1_B_LAB, "uunifrac")
adonis(rel_distance ~ DM, as(sample_data(ps1_B_LAB), "data.frame"))

#alpha-diversity
#Alpha diversity plot for Timepoint B -> Insert the respective variable
ps1_otu_table <- as.data.frame(ps1_B@otu_table)
ps1_metadata_table <- as.data.frame(ps1_B@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(Phase, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = Phase)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
print(Alphadiversity_plot)  
#Alpha diversity plot for Timepoint C -> Insert the respective variable
ps1_otu_table <- as.data.frame(ps1_C@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(Phase, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = Phase)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
print(Alphadiversity_plot)  

#Alpha diversity plot for LAB_B -> insert the respective variable
ps1_otu_table <- as.data.frame(ps1_B_LAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_B_LAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
print(Alphadiversity_plot)
print(Alphadiversity_table$PD)
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_LAB_B.csv")
#Alpha diversity plot for LAB_C -> insert the respective variable
ps1_otu_table <- as.data.frame(ps1_C_LAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C_LAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
print(Alphadiversity_plot)
print(Alphadiversity_table$PD)
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_LAB_C.csv") 
#Alpha diversity plot for SAB_B -> insert the respective variable
ps1_otu_table <- as.data.frame(ps1_B_SAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_B_SAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
print(Alphadiversity_plot)  
print(Alphadiversity_table$PD)
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_SAB_B.csv")
#Alpha diversity plot for SAB_C -> insert the respective variable
ps1_otu_table <- as.data.frame(ps1_C_SAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C_SAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
print(Alphadiversity_plot)
print(Alphadiversity_table$PD)
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_SAB_C.csv")

#Test alpha-diversity for significance by Kruskal-Wallis-test
ps1.adiv <- estimate_richness(ps1_C_SAB, measures = c("Chao1", "Shannon", "Observed", "InvSimpson"))
ps1.metadata <- as(sample_data(ps1_C_SAB), "data.frame")
head(ps1.metadata)

#Insert the respective phyloseq-object subset to load the appropriate OTU and metadata table
ps1_otu_table <- as.data.frame(ps1_C_SAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C_SAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)

#Add the alpha diversity indices to the metadata table 
ps1.metadata$Observed <- ps1.adiv$Observed 
ps1.metadata$Shannon <- ps1.adiv$Shannon
ps1.metadata$InvSimpson <- ps1.adiv$InvSimpson
ps1.metadata$Phyogenetic_diversity <- Alphadiversity_table$PD
#check if the last three rows are alpha diversity indices, i.e. Shannon, InvSimpson, Phylogenetic_diversity
colnames(ps1.metadata)

#Compare difference between samples/groups using richness
kruskal.pd <- kruskal.test(ps1.metadata$Phyogenetic_diversity ~ ps1.metadata$SilageType)
print(kruskal.pd)
#You can also compare two groups with adjustmen for multiple testing -> insert the respective variable, e.g. DM, WI, SA, etc. 
pairwise.wilcox.test(ps1.metadata$Phyogenetic_diversity, ps1.metadata$DMSA, p.adj = "BH")

#Alpha-diversity facetted plots
ps1_metadata_table <- as.data.frame(ps1@sam_data)
ps1_otu_table <- as.data.frame(ps1@otu_table)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD

plot.pd <- ggplot(ps1_metadata_table, aes(x=Timepoint, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = Phase), alpha=0.3) +
  scale_shape_manual(values = c("gold", "lightpink1", "springgreen3", "skyblue1", "darkorange", "brown",
                                "gray60", "gray40", "black", "dodgerblue", "purple")) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_classic() +
  labs(x = "Adaption time point", y = "Phylogenetic diversity") + facet_wrap(~Phase,scales = "free_y") + ylim(13,16)

plot.pd

#Barplot of genera with relative abundance >1% for LAB_B
OTU_filtered = prune_taxa(taxa_sums(ps1_B_LAB_relabund)>=1, ps1_B_LAB_relabund)
OTU_filtered_relabund <- transform_sample_counts(OTU_filtered,function(x)x/sum(x)*100)
p <- plot_bar(OTU_filtered_relabund, fill = "Genus") 
p <- p + facet_wrap(~Treatment, scales = "free_x", nrow = 1) + ylab("Relative abundance (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(panel.background = element_rect(fill = "white",
                                                                              colour = "white", size = 0.5, linetype = "solid"),
                                              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                              colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                                                                colour = "white"))
plot(p)

#Barplot of genera with relative abundance >1% for LAB_C
OTU_filtered = prune_taxa(taxa_sums(ps1_C_LAB_relabund)>=1, ps1_C_LAB_relabund)
OTU_filtered_relabund <- transform_sample_counts(OTU_filtered,function(x)x/sum(x)*100)
p <- plot_bar(OTU_filtered_relabund, fill = "Genus") 
p <- p + facet_wrap(~Treatment, scales = "free_x", nrow = 1) + ylab("Relative abundance (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(panel.background = element_rect(fill = "white",
                                                                              colour = "white", size = 0.5, linetype = "solid"),
                                              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                              colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                                                                colour = "white"))
plot(p)

#Barplot of genera with relative abundance >1% for SAB_B
OTU_filtered = prune_taxa(taxa_sums(ps1_B_SAB_relabund)>=1, ps1_B_SAB_relabund)
OTU_filtered_relabund <- transform_sample_counts(OTU_filtered,function(x)x/sum(x)*100)
p <- plot_bar(OTU_filtered_relabund, fill = "Genus") 
p <- p + facet_wrap(~Treatment, scales = "free_x", nrow = 1) + ylab("Relative abundance (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(panel.background = element_rect(fill = "white",
                                                                              colour = "white", size = 0.5, linetype = "solid"),
                                              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                              colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                                                                colour = "white"))
plot(p)

#Barplot of genera with relative abundance >1% for SAB_C
OTU_filtered = prune_taxa(taxa_sums(ps1_C_SAB_relabund)>=1, ps1_C_SAB_relabund)
OTU_filtered_relabund <- transform_sample_counts(OTU_filtered,function(x)x/sum(x)*100)

p <- plot_bar(OTU_filtered_relabund, fill = "Genus") 
p <- p + facet_wrap(~Treatment, scales = "free_x", nrow = 1) + ylab("Relative abundance (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(panel.background = element_rect(fill = "white",
                                                                              colour = "white", size = 0.5, linetype = "solid"),
                                              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                              colour = "white"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                                                                colour = "white"))
plot(p)

#ß-diversity plots for interactions of DMxSA and WIxSA at late adaption time point
#LAB_C
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA")+ geom_point(aes(color=DMSA),  size = 3) + theme(legend.title = element_blank()) + xlab("PCoA 1 (26.8%)") + ylab("PCoA 2 (19.7%)")
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA")+ geom_point(aes(color=WISA),  size = 3) + theme(legend.title = element_blank()) + xlab("PCoA 1 (26.8%)") + ylab("PCoA 2 (19.7%)")

#SAB_C
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA")+ geom_point(aes(color=DMSA),  size = 3) + theme(legend.title = element_blank()) + xlab("PCoA 1 (29.5%)") + ylab("PCoA 2 (15.2%)")
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA")+ geom_point(aes(color=WISA),  size = 3) + theme(legend.title = element_blank()) + xlab("PCoA 1 (29.5%)") + ylab("PCoA 2 (15.2%)")

