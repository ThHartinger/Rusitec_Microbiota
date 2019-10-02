
-
Title: "Rusitec_Microbiota_16S"
Author: "Thomas Hartinger"
Date: "2nd October 2019"
-


#####Load packages####
library(ggplot2)
library(ape)
library(vegan)
library(phyloseq)
library(microbiome)
library(mgcv)
library(class)
library(codetools)
library(picante)

####Biomfile, trefile, mapping file -> Create phyloseq-object####
biomfile <- import_biom("Rusitecmicrobiota.biom1")
biomfile_mapfile <- import_qiime_sample_data("Rusitecmicrobiota_mapping.txt") 
biomfile_merged <- merge_phyloseq(biomfile,biomfile_mapfile)
biomfile_tree <- read.tree("Rusitecmicrobiota.tre")
biomfile_merged <-merge_phyloseq(biomfile_merged,biomfile_tree)

####Prepare files for analysis####
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
ps1_B <- subset_samples(ps1, Timepoint == "Early")
ps1_B_relabund <- transform_sample_counts(ps1_B, function(x)x/sum(x))
ps1_C <- subset_samples(ps1, Timepoint == "Late")
ps1_C_relabund <- transform_sample_counts(ps1_C, function(x)x/sum(x))


####Alpha-diversity####
#Alpha diversity for LAB_B -> insert the respective variable/treatment
ps1_otu_table <- as.data.frame(ps1_B_LAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_B_LAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_LAB_B.csv") 
#Alpha diversity for LAB_C -> insert the respective variable/treatment
ps1_otu_table <- as.data.frame(ps1_C_LAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C_LAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_LAB_C.csv") 
#Alpha diversity for SAB_B -> insert the respective variable/treatment
ps1_otu_table <- as.data.frame(ps1_B_SAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_B_SAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_SAB_B.csv")
#Alpha diversity for SAB_C -> insert the respective variable/treatment
ps1_otu_table <- as.data.frame(ps1_C_SAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C_SAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
datatable(Alphadiversity_table)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
Alphadiversity_plot <- ggplot(ps1_metadata_table, aes(SilageType, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = SilageType)) + geom_point(size = 2) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw()
write.table(Alphadiversity_table$PD,file="Phylogenetic_diversity_SAB_C.csv")

#Test alpha-diversity for significance by Kruskal-Wallis-test
#Insert the respective phyloseq object (e.g. ps1_LAB) and the variable/treatment of interest
ps1.adiv <- estimate_richness(ps1_C_SAB, measures = c("Chao1", "Shannon", "Observed", "InvSimpson"))
ps1.metadata <- as(sample_data(ps1_C_SAB), "data.frame")
head(ps1.metadata)
#Insert the respective phyloseq-object/subset to load the appropriate OTU and metadata table
ps1_otu_table <- as.data.frame(ps1_C_SAB@otu_table)
ps1_metadata_table <- as.data.frame(ps1_C_SAB@sam_data)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
#Add the alpha diversity indices to the metadata table 
ps1.metadata$Observed <- ps1.adiv$Observed 
ps1.metadata$Shannon <- ps1.adiv$Shannon
ps1.metadata$InvSimpson <- ps1.adiv$InvSimpson
ps1.metadata$Phyogenetic_diversity <- Alphadiversity_table$PD
#Compare difference between samples/groups using richness
kruskal.pd <- kruskal.test(ps1.metadata$Phyogenetic_diversity ~ ps1.metadata$SilageType)
print(kruskal.pd)

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
  labs(x = "Adaption time point", y = "Phylogenetic diversity") + theme(legend.position="none") + facet_wrap(~Phase,scales = "free_y") + ylim(13,16)
plot.pd

ps1_metadata_table <- as.data.frame(ps1@sam_data)
ps1_otu_table <- as.data.frame(ps1@otu_table)
Alphadiversity_table <- pd(t(ps1_otu_table), biomfile_tree,include.root=F)
ps1_metadata_table$Phylogenetic_diversity <- Alphadiversity_table$PD
plot.pd <- ggplot(ps1_metadata_table, aes(x=Treatment, Phylogenetic_diversity)) + 
  geom_boxplot(aes(fill = Phase), alpha=0.3) +
  scale_shape_manual(values = c("gold", "lightpink1", "springgreen3", "skyblue1", "darkorange", "brown",
                                "gray60", "gray40", "black", "dodgerblue", "purple")) + 
  theme(axis.text.x = element_text(size=14, angle = 90)) + theme_classic() +
  labs(x = "Alfalfa silage", y = "Phylogenetic diversity") + facet_wrap(~PhaseTimepoint,scales = "free_y") + theme(legend.position="none")
plot.pd


####Beta-diversity####
#All samples unifrac weighted
plot_ordination(ps1, ordinate (ps1, distance = "unifrac", method="PCoA"), "samples", color="Timepoint", shape="Phase_Sucrose") +
  scale_shape_manual(values=c(0,15,2,17)) + theme(legend.position="right") + theme_light() + theme(legend.title=element_blank()) + geom_point(size = 2.0)

#All samples unifrac_unweighted
plot_ordination(ps1, ordinate (ps1, distance = "uunifrac", method="PCoA"), "samples", color="Timepoint", shape="Phase_Sucrose") +
  scale_shape_manual(values=c(0,15,2,17)) + theme(legend.position="right") + theme_light() + theme(legend.title=element_blank()) + geom_point(size = 2.0)

#ß-diversity plots for interactions of DMxSA and WIxSA at late adaption time point
#LAB_C
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA")+ geom_point(aes(color=DMSA),  size = 3) + theme(legend.title = element_blank()) 
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA")+ geom_point(aes(color=WISA),  size = 3) + theme(legend.title = element_blank())
#SAB_C
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA")+ geom_point(aes(color=DMSA),  size = 3) + theme(legend.title = element_blank())
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA")+ geom_point(aes(color=WISA),  size = 3) + theme(legend.title = element_blank()) 

#LAB_B samples unifrac all treatments
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#LAB_C samples unifrac all treatments
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_B samples unifrac all treatments
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_C samples unifrac all treatments
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "unifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#LAB_B samples uunifrac all treatments
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_LAB_relabund, ordinate(ps1_B_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#LAB_C samples uunifrac all treaments
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_LAB_relabund, ordinate(ps1_C_LAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_B samples uunifrac all treatments 
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_B_SAB_relabund, ordinate(ps1_B_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)
#SAB_C samples uunifrac all treatments
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DM", label="DM")+ geom_point(aes(color=DM),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WI", label="WI")+ geom_point(aes(color=WI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="SA", label="SA")+ geom_point(aes(color=SA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMSA", label="DMSA")+ geom_point(aes(color=DMSA),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="DMWI", label="DMWI")+ geom_point(aes(color=DMWI),  size = 3)
plot_ordination(ps1_C_SAB_relabund, ordinate(ps1_C_SAB_relabund, distance = "uunifrac", method="PCoA"),
                "samples", color="WISA", label="WISA")+ geom_point(aes(color=WISA),  size = 3)

#Test beta-diversity for significance by adonis for weighted unifrac
#Insert the respective phyloseq object (e.g. ps1_LAB) and the variable/treatment of interest, no p-value adjustment due to multiple testing is needed
rel_distance <- distance(ps1_B_LAB, "unifrac")
adonis(rel_distance ~ DM, as(sample_data(ps1_B_LAB), "data.frame"))

#Test beta-diversity for significance by adonis for unweighted unifrac
#Insert the respective phyloseq object (e.g. ps1_LAB) and the variable/treatment of interest, no p-value adjustment due to multiple testing is needed
rel_distance <- distance(ps1_B_LAB, "uunifrac")
adonis(rel_distance ~ DM, as(sample_data(ps1_B_LAB), "data.frame"))

#Remarks:
#Figure 1 was created using Excel from ps1_relabund
#Figures 5 and 6 were created using Canoco 5.11