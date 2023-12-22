# 16S sequencing results for pig skin study - took swabs from the pig skin to find out what is there  

library(phyloseq)
library(ggplot2)
library(microbiome)
library(tidyverse)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(decontam)
library(Maaslin2)
library(vegan)
library(metagMisc)
library(dplyr)
library(cowplot)
library(reshape2)
library(simarioV2)

# ps.NAME = a phyloseq object --> this is importing the datasets and associated metadata 
ps.Pig16S <- qza_to_phyloseq(features = "/Users/karindadelacruz/LK16S006/16SAnalysis/table-PigSkin2-dada2.qza",
                          taxonomy = "/Users/karindadelacruz/LK16S006/16SAnalysis/taxonomy-PigSkinSilva.qza",
                          metadata = "/Users/karindadelacruz/LK16S006/16SAnalysis/16S_Manufest_PigSkin.txt",
                          tree = "/Users/karindadelacruz/LK16S006/16SAnalysis/PigSkin-rooted-tree.qza")

# create dataframe with all OTUs that are in negative sample
ps.negativeotus16S <- subset_samples(ps.Pig16S, Dorsal_Ventral == "Neg")
ps.negativeotus16S <- subset_samples(ps.negativeotus16S, SampleType != "PBS no gauze")
negativeotus16S <- data.frame(ps.negativeotus16S@otu_table)
negativeotus16S <- negativeotus16S[rowSums(negativeotus16S[])>0,]
OTUstoRemove16S <- as.character(rownames(negativeotus16S))

# Removing contaminants identified in the negative swabs
ps.Decontam <- ps.Pig16S
OTUDoNotRemove16S <- setdiff(taxa_names(ps.Decontam), OTUstoRemove16S)
ps.Decontam <- prune_taxa(OTUDoNotRemove16S, ps.Decontam)

ps.Decontam <- subset_samples(ps.Decontam, Dorsal_Ventral != "Neg")
ps.Decontam <- subset_samples(ps.Decontam,Dorsal_Ventral != "Pos") 
ps.Decontam <- subset_taxa(ps.Decontam, !is.na(Phylum))
ps.Decontam <- subset_taxa(ps.Decontam, !is.na(Genus))
ps.Decontam <- subset_taxa(ps.Decontam, Kingdom != "d__Eukaryota")
ps.Decontam <- subset_taxa(ps.Decontam, Kingdom != "d__Archaea")
ps.Decontam <- subset_taxa(ps.Decontam, Genus != "Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Genus != "Mitochondria")
ps.Decontam <- subset_samples(ps.Decontam, SampleType == "Swab")
ps.Decontam <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = NULL)

write.csv(ps.Decontam@tax_table, "/Users/karindadelacruz/LK16S006/R_16S/TaxTable.csv")

ps.Decontam # THE CLEANED FILE 
sample_names(ps.Decontam)

# For maaslin
Pig16S.meta <- data.frame(sample_data(ps.Decontam))
Pig16S.meta$SampleID <- row.names(Pig16S.meta)

# Making compiled files of the total and relative abundance 
Pig16Sout <- psmelt(ps.Decontam)
write.csv(Pig16Sout, "/Users/karindadelacruz/LK16S006/R_16S/Pig16SabundanceTotalsOut.csv")
head(Pig16Sout)

relAll <- transform(ps.Decontam, "compositional")
relAllout <- psmelt(relAll)
write.csv(relAllout, "/Users/karindadelacruz/LK16S006/R_16S/Pig16SRelativeOUT.csv")

# Grouping
ps.ExGroup16S <- merge_samples(ps.Decontam, "ExGroup", fun = mean)

Combined_Man <- as.data.frame(X16S_CombinedManufest_PigSkin)
rownames(Combined_Man) <- Combined_Man$ExGroup

sample_data(ps.ExGroup16S) <- sample_data(Combined_Man)
relCombined <- transform(ps.ExGroup16S, "compositional")

# Relative abundance for WHOLE dataset
relCombined@sam_data$Timepoint <- factor(relCombined@sam_data$Timepoint, levels = c("1", "2","3","4", "5"))
plot_bar(relCombined,  x = "Timepoint", fill = "Phylum") + facet_wrap(~Pig + Dorsal_Ventral, scales = "free_x", nrow = 2)

# Relative abundance plots for TP3 
relCombinedAbudance <-psmelt(relCombined)
write.csv(relCombinedAbudance, "/Users/karindadelacruz/LK16S006/R_16S/HighTx0.1_OUT.csv")

#relative abundance plot for phyla >0.5% reads
abundance_other<- relCombinedAbudance
abundance_other <- subset.data.frame(abundance_other, Abundance > 0.0025)
abundance_other<- abundance_other %>% mutate(Phylum = ifelse(Abundance < 0.005, "Other", Phylum))
write.csv(abundance_other, "/Users/karindadelacruz/LK16S006/R_16S/abundance_other.csv")

abundance_other$Phylum <- factor(abundance_other$Phylum, levels = c("Other", "Actinobacteriota", "Bacteroidota", "Campilobacterota", "Desulfobacterota", "Elusimicrobiota", "Firmicutes", "Fusobacteriota", "Patescibacteria","Proteobacteria", "Spirochaetota", "Verrucomicrobiota"))
phylumcolors <- c("grey", "#665191", "#a05195", "#d45087", "#D48BA9", "#D4AEBE", "#003F5C", "#f95d6a", "#F999A1", "#ff7c43", "#FFA601", "#FFE2AB")(16)
abundance_other$"Dorsal_Ventral" <- factor(abundance_other$"Dorsal_Ventral", levels = c("D", "V", "F"))
abundance_other$"Dorsal_Ventral" <- fct_recode(abundance_other$"Dorsal_Ventral", "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")
                                     
ggplot(abundance_other, aes(x=Pig, y = Abundance, fill = Phylum, color = Phylum))+
  labs(x="Swab Type", element_text(size=8))+
  geom_bar(aes(fill = Phylum, color = Phylum), stat="identity", position="fill") +
  scale_x_discrete(labels=c("Pig 1", "Pig 2", "Pig 1", "Pig 2", "Pig 1", "Pig 2"))+
  facet_wrap("Dorsal_Ventral", scales = "free_x", ncol = 3)+
  theme_bw(base_size = 8)+
  theme(legend.title=element_text(size=8), 
        legend.text=element_text(size=6))+
  scale_fill_manual(values = phylumcolors) +
  scale_color_manual(values = phylumcolors)+
  ylab("Relative Abundance")+
  theme_light()
  
#relative abundance plot for genus >0.5% reads
abundance_other<- relCombinedAbudance
abundance_other <- subset.data.frame(abundance_other, Abundance > 0.0025)
abundance_other<- abundance_other %>% mutate(Genus = ifelse(Abundance < 0.005, "Other", Genus))
write.csv(abundance_other, "/Users/karindadelacruz/LK16S006/R_16S/abundance_othergenus.csv")

genusabudancetable<-table(abundance_other$Genus, abundance_other$Phylum) # table of genera in each phyla
write.csv(genusabudancetable, "/Users/karindadelacruz/LK16S006/R_16S/GenusAbundanceTable.csv")


genus_colors <- c("gray", "black",
                  colorRampPalette(c("#665191", "#AA87F2", "#E1D9F2"))(10),
                  colorRampPalette(c("#a05195", "#E172D1", "#E1C0DC"))(14), 
                  colorRampPalette(c("#d45087", "#D46E98"))(2), 
                  colorRampPalette(c("#D48BA9"))(1),
                  colorRampPalette(c("#D4AEBE"))(1),
                  colorRampPalette(c("#003F5C", "#009DE5", "#BAD7E5", "#2f4b7c", "#5B92EF", "#BFD1EF"))(52),
                  colorRampPalette(c("#f95d6a"))(1),
                  colorRampPalette(c("#F999A1"))(1),
                  colorRampPalette(c("#ff7c43", "#FFA67F", "#FFD8C7"))(6),
                  colorRampPalette(c("#FFA601"))(1), 
                  colorRampPalette(c("#FFE2AB"))(1)
                  )

abundance_other$"Dorsal_Ventral" <- factor(abundance_other$"Dorsal_Ventral", levels = c("D", "V", "F"))
abundance_other$"Dorsal_Ventral" <- fct_recode(abundance_other$"Dorsal_Ventral", "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")

abundance_other$Genus<- factor(abundance_other$Genus, levels = c("Other", "uncultured",
                                                                           "Bifidobacterium", "Brachybacterium", "Collinsella", "Corynebacterium", "Dermabacter", "Enterorhabdus", "Kocuria", "Micrococcus", "Olsenella", "Rothia", #Actinobacteriota (10)
                                                                           "Alistipes", "Alloprevotella", "Bacteroides", "Chryseobacterium", "dgA-11_gut_group", "Empedobacter", "Muribaculaceae", "Parabacteroides", "Porphyromonas", "Prevotella", "Prevotellaceae_NK3B31_group", "Prevotellaceae_UCG-001", "Prevotellaceae_UCG-003", "Rikenellaceae_RC9_gut_group", #Bacteroidota (14)
                                                                           "Arcobacter", "Campylobacter", #Campilobacterota (2)
                                                                 "Desulfovibrio", #Desulfobacterota (1)
                                                                 "Elusimicrobium", #Elusimicrobiota (1)
                                                                           "[Eubacterium]_coprostanoligenes_group", "[Eubacterium]_ruminantium_group", "[Ruminococcus]_torques_group", "Acidaminococcus", "Aerococcus", "Anaerococcus", "Anaerovibrio", "Bacillus", "Candidatus_Soleaferrea", "Christensenellaceae_R-7_group", "Clostridia_UCG-014", "Clostridia_vadinBB60_group", "Clostridium_sensu_stricto_1", "Coprococcus", "Dialister", "Dorea", "Enterococcus",  "Erysipelatoclostridium", "Erysipelotrichaceae_UCG-009", "Facklamia", "Faecalicoccus", "Family_XIII_AD3011_group", "Fenollaria", "Frisingicoccus", "Incertae_Sedis", "Jeotgalibaca", "Jeotgalicoccus", "Kurthia", "Lachnospiraceae_NK4A136_group", "Lactobacillus", "Megasphaera", "Mitsuokella", "NK4A214_group", "Phascolarctobacterium",  "Romboutsia", "Roseburia", "Ruminococcus", "S5-A14a", "Solibacillus", "Solobacterium", "Staphylococcus", "Streptococcus",  "Terrisporobacter", "Trichococcus", "Turicibacter", "UCG-001", "UCG-002",  "UCG-005", "UCG-008", "UCG-009", "UCG-010", "W5053", #Firmicutes (52)
                                                                           "Fusobacterium", #Fusobacteriota (1)
                                                                 "TM7a", #Patescibacteria (1)
                                                                           "Acinetobacter", "Escherichia-Shigella", "Lysobacter", "Moraxella", "Psychrobacter", "Succinivibrio", #Proteobacteria (6)
                                                                           "Treponema", #Spirochaetota (1)
                                                                 "WCHB1-41"))#Verrucomicrobiota(1)
                                                 
ggplot(abundance_other, aes(x = Pig , y = Abundance, fill = Genus, color = Genus))+
  geom_bar(aes(fill = Genus, color = Genus), stat= "identity", position = "fill")+
  labs(x="Swab Type", element_text(size=8))+
  theme(legend.title=element_text(size=8), 
        legend.text=element_text(size=6))+
  scale_x_discrete(labels=c("Pig 1", "Pig 2", "Pig 1", "Pig 2", "Pig 1", "Pig 2"))+
  facet_wrap("Dorsal_Ventral", scales = "free_x", ncol = 3)+
  theme_bw(base_size = 8)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = genus_colors) +
  scale_color_manual(values = genus_colors)

# FOR ALPHA and BETA TABLES
# Create subset for swab samples
ps.Swabs <- subset_samples(relAll, SampleType != "Gauze in PBS")

# ALPHA ABUNDANCE - Swabs
ps.Swabs@sam_data$Timepoint <- factor(ps.Swabs@sam_data$Timepoint, levels = c("1", "2","3","4", "5"))
plot_richness(ps.Swabs, x = "Timepoint", measures=c("Shannon"), color = "Pig") + 
  facet_wrap("Dorsal_Ventral", scales = "free_x", ncol = 2) + 
  theme_bw()

tabSwab <-microbiome::alpha(ps.Swabs, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                  "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                  "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
write.csv(tabSwab, "/Users/karindadelacruz/LK16S006/R_16S/Pig16SSwabalphaTable.csv")

# MAKE A METADATA TABLE with the alpha diversity metrics and import it - Swab
AlphaSwab <- Pig16SSwabalphaTable # the alpha diversity table from above with incorporated metadata
AlphaSwab$Timepoint <- factor(AlphaSwab$Timepoint, levels = c("1", "2","3","4", "5"))

#use this alpha diversity
AlphaSwab$Dorsal_Ventral <- factor(AlphaSwab$Dorsal_Ventral , levels = c("D", "V", "F"))
AlphaSwab$Dorsal_Ventral <- fct_recode(AlphaSwab$Dorsal_Ventral, "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")
ggplot(AlphaSwab, aes (x = Dorsal_Ventral, y = diversity_shannon, color = Pig, fill = Pig))+
  geom_point(size = 3, position = position_dodge(width = 0.4))+
  stat_summary(geom = "point", shape = 3, fun = "mean", size = 8, position = position_dodge(width = 0.4))+
  scale_color_manual(values = c("#7A5195", "#EE5675"))+
  theme_bw(base_size = 8)+
  ylab("Total Community\nShannon Alpha Diversity") +
  xlab("Swab Type") +
  theme(legend.position = "bottom")

# Beta diversity 
ps.SwabsBeta <- subset_samples(ps.Decontam, SampleType != "Gauze in PBS")

min(sample_sums(ps.SwabsBeta))# minimum sample read is 8675
median(sample_sums(ps.SwabsBeta)) #18637

table(sample_sums(ps.SwabsBeta))

tabSwab16S <- t(otu_table(ps.SwabsBeta))

ps.rarefied = rarefy_even_depth(ps.SwabsBeta, rngseed=1, sample.size=15000, replace=F) # THIS IS THE ONE TO GO WITH

# Beta diversity - Bray Curtis - use!
GP = ps.rarefied
GP.ord <- ordinate(GP, "PCoA",  "unifrac", weighted=T) 
allGroupsColors <- c("#003F5C", "#bc5090", "#ffa600")

sample_data(GP)$Dorsal_Ventral <- factor(sample_data(GP)$Dorsal_Ventral , levels = c("D", "V", "F"))
sample_data(GP)$Dorsal_Ventral <- fct_recode(sample_data(GP)$Dorsal_Ventral, "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")

plot_ordination(GP, GP.ord, type="samples", color="Dorsal_Ventral") + 
  stat_ellipse(type = "t", level = 0.95)+
  geom_point(aes(shape =Pig),size=3) + #ggtitle("Bray Curtis Beta Diversity")+
  #scale_shape_manual(values = c(8,21))+
  theme_bw(base_size = 8)+
  scale_color_manual(values = allGroupsColors,
                     name = "Swab Type") +
  scale_shape_discrete(name="Animal",
                       labels=c("Pig 1", "Pig 2")) +
  theme(legend.position = "bottom", legend.box = "vertical")

# permanova
sampledf <- data.frame(sample_data(ps.rarefied))
pig16sunifrac <- phyloseq::distance(ps.rarefied, method = "unifrac")
unifrac_clust <-hclust(pig16sunifrac)
plot(unifrac_clust)
adonis2(pig16sunifrac ~ Dorsal_Ventral, by= "margin", data = sampledf, permutations = 9999) #
Perm1 <- adonis2(pig16sunifrac ~ Dorsal_Ventral, by= "margin", data = sampledf, permutations = 9999) #
perm1df <- as.data.frame(Perm1$aov.tab)
write.csv(perm1df, "/Users/karindadelacruz/LK16S006/R_16S/PERMANOVA/Dorsal_Ventral_PERMANOVA.csv")

sampledf2 <- data.frame(sample_data(ps.rarefied))
sampledf2 <-subset.data.frame(sampledf2, SampleType == "Swab")
sampledf2$Dorsal_Ventral <- factor(sampledf2$Dorsal_Ventral, levels = c("D", "V", "F"))
sampledf2 <- sampledf2 %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))
pig16sunifrac2 <- phyloseq::distance(ps.rarefied, method = "unifrac")
unifrac_clust2 <-hclust(pig16sunifrac2)
plot(unifrac_clust2)
adonis2(pig16sunifrac2 ~ SkinFecal, by= "margin", data = sampledf2, permutations = 9999) #
Perm2 <- adonis2(pig16sunifrac2 ~ SkinFecal, by= "margin", data = sampledf2, permutations = 9999) #
perm2df <- as.data.frame(Perm2$aov.tab)
write.csv(perm1df, "/Users/karindadelacruz/LK16S006/R_16S/PERMANOVA/Skin_Fecal_PERMANOVA.csv")

# MAASLIN
# Running Maaslin2 for determining differential abundance - similar to above but theretically more accurate
# https://huttenhower.sph.harvard.edu/maaslin/
# input = two .tsv (or .csv) files. the first is the relative abundance table (####taxaRElAbundance.tsv).
# the seccond file is the .tsv of the metadata file (16S_Manufest_PigSkin.tsv)

input_OTU_data <- as.data.frame(Pig16SOTUs)
rownames(input_OTU_data) <- input_OTU_data[,1]

input_genus_data <- as.data.frame(Pig16Sgenus)
rownames(input_genus_data) <-input_genus_data[,1]

input_phyla_data <- as.data.frame(Pig16SPhyla)
rownames(input_phyla_data) <-input_phyla_data[,1]

input_meta_file <- Pig16S.meta
input_meta_file <-subset.data.frame(input_meta_file, SampleType == "Swab")
input_meta_file$Dorsal_Ventral <- factor(input_meta_file$Dorsal_Ventral, levels = c("D", "V", "F"))
input_meta_file <- input_meta_file %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))

SkinOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Skin")
FecalOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Fecal")

# evaluating fecal vs. skin 
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = input_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/SkinFecal_OTU_maaslin",
  fixed_effects = c("SkinFecal"),
  reference = c("Fecal", "Skin"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = input_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/SkinFecal_genus_maaslin",
  fixed_effects = c("SkinFecal"),
  reference = c("Fecal", "Skin"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = input_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/SkinFecal_phyla_maaslin",
  fixed_effects = c("SkinFecal"),
  reference = c("Fecal", "Skin"))

# evaluating dorsal vs. ventral 
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/DorsalVentral_OTU_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/DorsalVentral_genus_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/DorsalVentral_phyla_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

# evaluating pig 1 skin vs. pig 2 skin
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/PigSkinPigSkin_OTU_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/PigSkinPigSkin_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/PigSkinPigSkin_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

# evaluating pig 1 fecal vs. pig 2 fecal
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = FecalOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/PigFecalPigFecal_OTU_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = FecalOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/PigFecalPigFecal_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = FecalOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_16S/maaslin/PigFecalPigFecal_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

## SPIEC-EASI - skin only
# load packages necessary
library(devtools)
library(SpiecEasi)
library(igraph)
library(seqtime)
library(taxonomizr)
library(phylosmith)

# subset the samples to only include skin only and remove low abundance phylum (as specified in the abudance plot)
ps.skinonly <- subset_samples(relCombined, Dorsal_Ventral %in% c("D","V"))

# extract the otu and taxonomy tables from the final phyloseq object
ps.skinonly #finalized phyloseq object
otus16S=otu_table(ps.skinonly)
taxa16S=tax_table(ps.skinonly)
otus16S <- t(otus16S)

# filter the OTU matrix to remove zeros, but keep the sum to not change the overall sample counts
filterobj16S=filterTaxonMatrix(otus16S, minocc=1, keepSum = TRUE, return.filtered.indices = TRUE) # otus16S must have taxa as rows and samples as columns
otus16S.f=filterobj16S$mat

# introduce dummy taxonomic information for the last row of the Taxa table
taxa16S.f=taxa16S[setdiff(1:nrow(taxa16S),filterobj16S$filtered.indices),]
#dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
#taxa16S.f=rbind(taxa16S.f,dummyTaxonomy)
rownames(taxa16S.f)[nrow(taxa16S.f)]="0"
rownames(otus16S.f)[nrow(otus16S.f)]="0"

# assemble a new phyloseq object with the filtered OTU and taxonomy table
updatedotus16S=otu_table(otus16S.f, taxa_are_rows = TRUE)
view(updatedotus16S)
updatedtaxa16S=tax_table(taxa16S.f)
view(updatedtaxa16S)
phyloseqobj16S.f=phyloseq(updatedotus16S, updatedtaxa16S)

# make a phyloseq object that contains enough data but not too much
ps.spieceasi16S <- phyloseq_filter_prevalence(phyloseqobj16S.f, prev.trh = 0.1, abund.trh = 10, threshold_condition = "OR")

# apply SPIEC-EASI - 0.01
spiec16S.out2=spiec.easi(ps.spieceasi16S, method="mb", icov.select.params=list(rep.num=20, thresh=0.01))
# convert the output of SPIEC-EASI into a network and plot it with colors on a class level
spiec16S.graph2=adj2igraph(getRefit(spiec16S.out2), vertex.attr=list(name=taxa_names(ps.spieceasi16S)))
plot_network(spiec16S.graph, ps.spieceasi16S, type="taxa", color="Phylum", label=NULL, point_size = 4, line_weight=0.5)

# count positive and negative edges to extract the regression coefficients from the output
betaMat16S=as.matrix(symBeta(getOptBeta(spiec16S.out2)))

# obtain the number of edges in the network from the matrix of regression coefficients
positive16S=length(betaMat16S[betaMat16S>0])/2
negative16S=length(betaMat16S[betaMat16S<0])/2
total16S=length(betaMat16S[betaMat16S!=0])/2

# cluster SPIEC-EASI network and list the taxa present in each cluster
clusters16S=cluster_fast_greedy(spiec16S.graph2)
clusterOneIndices16S=which(clusters16S$membership==1)
clusterOneOtus16S=clusters16S$names[clusterOneIndices16S]
clusterTwoIndices16S=which(clusters16S$membership==2)
clusterTwoOtus16S=clusters16S$names[clusterTwoIndices16S]

# summarize and sort the representatives of the different classes
sort(table(getTaxonomy(clusterOneOtus16S,taxa16S.f,useRownames = TRUE)),decreasing = TRUE)
sort(table(getTaxonomy(clusterTwoOtus16S,taxa16S.f,useRownames = TRUE)),decreasing = TRUE)

# export while preserving the OTU identifiers
write.graph(spiec16S.graph2,file=file.path("/Users/karindadelacruz/LK16S006/spieceasinetworks/spieceasi16Snetwork.txt"),format="ncol")

# export the taxonomic information
write.table(taxa16S.f,file=file.path("/Users/karindadelacruz/LK16S006/spieceasinetworks/spieceasi16S.txt"),sep="\t", quote=FALSE)

# open in cytoscape
