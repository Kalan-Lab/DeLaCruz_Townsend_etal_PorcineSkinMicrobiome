# De La Cruz and Townsend. Porcine Skin Bacterial Fungal Interactions R code 

# Packages 
library(ggplot2)
library(pheatmap)
library(phyloseq)
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

# working directory
setwd(dir = "/Users/liztown/Documents/KalanLab/Papers/PigBactFun.Karinda/Final.Rcode")


# Figure 1 A-D porcine mycobiome (fungi)
ITSSingle_CombinedManufest_PigSkin<-as.data.frame(read.csv("./ITS/ITSSingle_CombinedManufest_PigSkin.csv"))
ps.PigITS <- qza_to_phyloseq(features = "./ITS/table-PigSkinITSSingle-dada2.qza",
                             taxonomy = "./ITS/taxonomy-PigSkinITSSingle.qza",
                             metadata = "./ITS/ITSSingle_Manufest_PigSkin.txt",
                             tree = "./ITS/PigSkinITSSingle-rooted-tree.qza")

ps.negativeotusITS <- subset_samples(ps.PigITS, Dorsal_Ventral == "Neg")
ps.negativeotusITS <- subset_samples(ps.negativeotusITS, SampleType != "PBS no gauze")
negativeotusITS <- data.frame(ps.negativeotusITS@otu_table)
negativeotusITS <- negativeotusITS[rowSums(negativeotusITS[])>0,]
ITSOTUstoRemove <- as.character(rownames(negativeotusITS))

ps.DecontamITS <- ps.PigITS
ITSOTUDoNotRemove <- setdiff(taxa_names(ps.DecontamITS), ITSOTUstoRemove)
ps.DecontamITS <- prune_taxa(ITSOTUDoNotRemove, ps.DecontamITS)
ps.DecontamITS <- subset_samples(ps.DecontamITS, Dorsal_Ventral != "Neg")
ps.DecontamITS <- subset_samples(ps.DecontamITS,Dorsal_Ventral != "Pos") 
ps.DecontamITS# THE CLEANED FILE 

ps.DecontamITS <- subset_taxa(ps.DecontamITS, !is.na(Phylum))
ps.DecontamITS <- subset_taxa(ps.DecontamITS, !is.na(Genus))

PigITSout <- psmelt(ps.DecontamITS)
#write.csv(PigITSout, "./ITS/PigITSabundanceTotalsOut.csv")

relAllITS <- microbiome::transform(ps.DecontamITS, "compositional")
relAlloutITS <- psmelt(relAllITS)
#write.csv(relAlloutITS, "./ITS/PigITSRelativeOUT.csv")

ps.ExGroupITS <- merge_samples2(ps.DecontamITS, "ExGroup", fun_otu = mean)

Combined_ManITS <- as.data.frame(ITSSingle_CombinedManufest_PigSkin)
rownames(Combined_ManITS) <- Combined_ManITS$ExGroup

sample_data(ps.ExGroupITS) <- sample_data(Combined_ManITS)
relCombinedITS <- microbiome::transform(ps.ExGroupITS, "compositional")

#Figure 2A. fungi relative abundance plot for phyla >0.5% reads
relCombinedITSAbundance <-psmelt(relCombinedITS)
relCombinedITS_other<- relCombinedITSAbundance %>% mutate(Phylum = ifelse(Abundance < 0.005, "Other", Phylum))
relCombinedITS_other$Phylum <- factor(relCombinedITS_other$Phylum, levels = c("Other", "unidentified", "Ascomycota", "Basidiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota"))
phylumITScolors <- c("grey", "black", "#003f5c", "#665191", "#d45087","#ff7c43", "#ffa600")
relCombinedITS_other$"Dorsal_Ventral" <- factor(relCombinedITS_other$"Dorsal_Ventral", levels = c("D", "V", "F"))
relCombinedITS_other$"Dorsal_Ventral" <- fct_recode(relCombinedITS_other$"Dorsal_Ventral", "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")

ggplot(relCombinedITS_other, aes(x=Pig, y = Abundance, fill = Phylum, color = Phylum))+
  geom_bar(aes(fill = Phylum, color = Phylum),stat="identity", position="fill") +
  labs(x="Swab Type", element_text(size=8))+
  facet_wrap("Dorsal_Ventral", scales = "free_x", ncol = 3)+
  scale_x_discrete(labels=c("Pig 1", "Pig 2", "Pig 1", "Pig 2", "Pig 1", "Pig 2"))+
  theme_bw(base_size = 8)+
  theme(legend.title=element_text(size=8), 
        legend.text=element_text(size=6))+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = phylumITScolors) +
  scale_color_manual(values = phylumITScolors)

# figure 2B fungi relative abundance plot for genus >0.5% reads
relCombinedITSGenus <- relCombinedITSAbundance
relCombinedITSGenus<- relCombinedITSGenus %>% mutate(Genus = ifelse(Abundance < 0.005, "Other", Genus))
table(relCombinedITSGenus$Genus, relCombinedITSGenus$Phylum) # table of genera in each phyla

genusITScolors <- c("gray","black",
                    colorRampPalette(c("#003F5C", "#009DE5", "#BAD7E5", "#2f4b7c", "#5B92EF", "#BFD1EF"))(33), 
                    colorRampPalette(c("#665191", "#AA87F2", "#E1D9F2", "#a05195", "#E172D1", "#E1C0DC"))(16), 
                    colorRampPalette(c("#d45087"))(1))

relCombinedITSGenus$"Dorsal_Ventral" <- factor(relCombinedITSGenus$"Dorsal_Ventral", levels = c("D", "V", "F"))
relCombinedITSGenus$Genus<- factor(relCombinedITSGenus$Genus, levels = c("Other", "unidentified",
                                                                         "Alternaria", "Aspergillus", "Aureobasidium", "Bipolaris", "Candida", "Cephaliophora", "Cladosporium", "Cleistothelebolus", "Coniothyrium", "Debaryomyces", "Dipodascus", "Epicoccum", "Fusarium", "Hypomyces", "Kazachstania", "Lectera", "Meyerozyma", "Neodevriesia", "Neosetophoma", "Nigrospora", "Ophiosphaerella", "Paraphaeosphaeria", "Penicillium", "Saccharomyces", "Sarocladium", "Sclerotinia", "Talaromyces", "Torulaspora", "Trichoderma", "Verrucocladosporium", "Wickerhamiella", "Wickerhamomyces", "Xeromyces", #Ascomycota (33)
                                                                         "Apiotrichum", "Cerrena", "Cryptococcus", "Cutaneotrichosporon", "Hannaella", "Holtermanniella", "Lycoperdon", "Malassezia", "Naganishia", "Rigidoporus", "Saitozyma", "Sistotrema", "Sterigmatomyces", "Trichosporon", "Vishniacozyma", "Wallemia", #Basidiomycota (16)
                                                                         "Mortierella" #Mortierellomycota (1)
)) 
relCombinedITSGenus$"Dorsal_Ventral" <- fct_recode(relCombinedITSGenus$"Dorsal_Ventral", "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")

ggplot(relCombinedITSGenus, aes(x = Pig, y = Abundance, fill = Genus, color = Genus))+
  geom_bar(aes(fill = Genus, color = Genus), stat= "identity", position = "fill")+
  labs(x="Swab Type", element_text(size=8))+
  theme(legend.title=element_text(size=8), 
        legend.text=element_text(size=6))+
  scale_x_discrete(labels=c("Pig 1", "Pig 2", "Pig 1", "Pig 2", "Pig 1", "Pig 2"))+
  facet_wrap("Dorsal_Ventral", scales = "free_x", ncol = 3)+
  theme_bw(base_size = 8)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = genusITScolors)+
  scale_color_manual(values = genusITScolors)

# Figure 2C beta diversity 
ps.SwabsBetaITS <- subset_samples(ps.DecontamITS, SampleType != "Gauze in PBS")
min(sample_sums(ps.SwabsBetaITS))# minimum sample read is 29
median(sample_sums(ps.SwabsBetaITS)) #2185
table(sample_sums(ps.SwabsBetaITS))
tabSwabITS <- t(otu_table(ps.SwabsBetaITS))
ps.rarefiedITS = rarefy_even_depth(ps.SwabsBetaITS, rngseed=1, sample.size=600, replace=F)

GPITS = ps.rarefiedITS
GP.ordITS <- ordinate(GPITS, "NMDS", "bray") # For Bray Curtis of the rareified dataset
allGroupsColors <- c("#003F5C", "#bc5090", "#ffa600")

sample_data(GPITS)$Dorsal_Ventral <- factor(sample_data(GPITS)$Dorsal_Ventral , levels = c("D", "V", "F"))
sample_data(GPITS)$Dorsal_Ventral <- fct_recode(sample_data(GPITS)$Dorsal_Ventral, "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")

GP.ordITS <- ordinate(GPITS, "PCoA",  "unifrac", weighted = T) 
plot_ordination(GPITS, GP.ordITS, type="samples", color="Dorsal_Ventral") + 
  stat_ellipse(type = "t", level = 0.8)+
  geom_point(aes(shape =Pig),size=3) + #ggtitle("Bray Curtis Beta Diversity")+
  #scale_shape_manual(values = c(8,21))+
  theme_bw(base_size = 8)+
  scale_color_manual(values = allGroupsColors,
                     name = "Swab Type") +
  scale_shape_discrete(name="Animal",
                       labels=c("Pig 1", "Pig 2")) +
  theme(legend.position = "bottom", legend.box = "vertical")

sampledfITS <- data.frame(sample_data(ps.rarefiedITS))
pigITSunifrac <- phyloseq::distance(ps.rarefiedITS, method = "unifrac")
adonis2(pigITSunifrac ~ Dorsal_Ventral, by= "margin", data = sampledfITS, permutations = 9999) #

sampledf2ITS <- data.frame(sample_data(ps.rarefiedITS))
sampledf2ITS <-subset.data.frame(sampledf2ITS, SampleType == "Swab")
sampledf2ITS$Dorsal_Ventral <- factor(sampledf2ITS$Dorsal_Ventral, levels = c("D", "V", "F"))
sampledf2ITS <- sampledf2ITS %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))
pigITSunifrac2 <- phyloseq::distance(ps.rarefiedITS, method = "unifrac")
adonis2(pigITSunifrac2 ~ SkinFecal, by= "margin", data = sampledf2ITS, permutations = 9999) #

## Fig 2D Alpha abundance for fungal communities
ps.Swabs.ITS <- subset_samples(relAllITS, SampleType != "Gauze in PBS")
ps.ITSs.meta<-as.data.frame(sample_data(ps.Swabs.ITS))
tabSwabITS <- microbiome::alpha(ps.Swabs.ITS, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                     "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                     "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))

AlphaSwabITS <- merge(tabSwabITS,ps.ITSs.meta, by = "row.names")

AlphaSwabITS$Dorsal_Ventral <- factor(AlphaSwabITS$Dorsal_Ventral , levels = c("D", "V", "F"))
AlphaSwabITS$Dorsal_Ventral <- fct_recode(AlphaSwabITS$Dorsal_Ventral, "Dorsal" = "D", "Ventral" = "V", "Fecal" = "F")

ggplot(AlphaSwabITS, aes (x = Dorsal_Ventral, y = diversity_shannon, color = Pig, fill = Pig))+
  geom_point(size = 4, position = position_dodge(width = 0.4))+
  stat_summary(geom = "point", shape = 3, fun = "mean", size = 8, position = position_dodge(width = 0.4))+
  scale_color_manual(values = c("#7A5195", "#EE5675"))+
  theme_bw(base_size = 8)+
  ylab("Total Community\nShannon Alpha Diversity") +
  xlab("Swab Type")+
  theme(legend.position = "bottom")


# Supplemental Figure  2A (and additional maaslin analysis on the fungal data set )
PigITS.meta <- data.frame(sample_data(ps.DecontamITS))
PigITS.meta$SampleID <- row.names(PigITS.meta)

PigITSGenus <-read.csv("./PigITSGenus.csv")
input_genus_data <- as.data.frame(PigITSGenus)
rownames(input_genus_data) <-input_genus_data[,1]

PigITSPhyla <-read.csv("./PigITSPhyla.csv")
input_phyla_data <- as.data.frame(PigITSPhyla)
rownames(input_phyla_data) <-input_phyla_data[,1]

input_meta_file <- PigITS.meta
input_meta_file <-subset.data.frame(input_meta_file, SampleType == "Swab")
input_meta_file$Dorsal_Ventral <- factor(input_meta_file$Dorsal_Ventral, levels = c("D", "V", "F"))
input_meta_file <- input_meta_file %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))

SkinOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Skin")
FecalOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Fecal")

# evaluating dorsal vs. ventral 
fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "./ITS/maaslin/DorsalVentral_genus_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "./ITS/maaslin/DorsalVentral_phyla_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

# evaluating pig 1 skin vs. pig 2 skin
fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "./ITS/maaslin/PigSkinPigSkin_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "./ITS/maaslin/PigSkinPigSkin_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

# evaluating pig 1 fecal vs. pig 2 fecal
fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = FecalOnly_meta_file,
  output = "./ITS/maaslin/PigFecalPigFecal_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = FecalOnly_meta_file,
  output = "./ITS/maaslin/PigFecalPigFecal_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))


# Figure 2 E-H porcine microbiome (bacteria)
Combined_Man <-read.csv("./16S/16S_CombinedManufest_PigSkin.csv")
Combined_Man <- as.data.frame(Combined_Man)
rownames(Combined_Man) <- Combined_Man$ExGroup

ps.Pig16S <- qza_to_phyloseq(features = "./16S/table-PigSkin2-dada2.qza",
                             taxonomy = "./16S/taxonomy-PigSkinSilva.qza",
                             metadata = "./16S/16S_Manufest_PigSkin.txt",
                             tree = "./16S/PigSkin-rooted-tree.qza")

ps.negativeotus16S <- subset_samples(ps.Pig16S, Dorsal_Ventral == "Neg")
ps.negativeotus16S <- subset_samples(ps.negativeotus16S, SampleType != "PBS no gauze")
negativeotus16S <- data.frame(ps.negativeotus16S@otu_table)
negativeotus16S <- negativeotus16S[rowSums(negativeotus16S[])>0,]
OTUstoRemove16S <- as.character(rownames(negativeotus16S))

ps.Decontam <- ps.Pig16S # Removing contaminants identified in the negative swabs
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
ps.Decontam@sam_data$PigDVF <-paste(ps.Decontam@sam_data$Pig, ps.Decontam@sam_data$Dorsal_Ventral)
ps.Decontam <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = NULL)

ps.Decontam # THE CLEANED FILE 

Pig16Sout <- psmelt(ps.Decontam)# Making compiled files of the total and relative abundance 
#write.csv(Pig16Sout, "./16S/Pig16SabundanceTotalsOut.csv")
relAll <- microbiome::transform(ps.Decontam, "compositional")
relAllout <- psmelt(relAll)
#write.csv(relAllout, ".16S/Pig16SRelativeOUT.csv")

ps.ExGroup16S <- merge_samples2(ps.Decontam,group = "ExGroup", fun_otu = mean)# Grouping
sample_data(ps.ExGroup16S) <- sample_data(Combined_Man)
relCombined <- microbiome::transform(ps.ExGroup16S, "compositional")

relCombinedAbudance <-psmelt(relCombined)
#write.csv(relCombinedAbudance, "/Users/karindadelacruz/LK16S006/R_16S/HighTx0.1_OUT.csv")

# Figure 2E. relative abundance plot for phyla >0.5% reads
abundance_other<- relCombinedAbudance
abundance_other <- subset.data.frame(abundance_other, Abundance > 0.0025)
abundance_other<- abundance_other %>% mutate(Phylum = ifelse(Abundance < 0.005, "Other", Phylum))
abundance_other$Phylum <- factor(abundance_other$Phylum, levels = c("Other", "Actinobacteriota", "Bacteroidota", "Campilobacterota", "Desulfobacterota", "Elusimicrobiota", "Firmicutes", "Fusobacteriota", "Patescibacteria","Proteobacteria", "Spirochaetota", "Verrucomicrobiota"))
phylumcolors <- c("gray", "#665191", "#a05195", "#d45087", "#D48BA9", "#D4AEBE", "#003F5C", "#f95d6a", "#F999A1", "#ff7c43", "#FFA601", "#FFE2AB")
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

# Figure 2F. relative abundance plot for genus >0.5% reads
abundance_other<- relCombinedAbudance
abundance_other <- subset.data.frame(abundance_other, Abundance > 0.0025)
abundance_other<- abundance_other %>% mutate(Genus = ifelse(Abundance < 0.005, "Other", Genus))
#write.csv(abundance_other, "./16S/abundance_othergenus.csv")

genusabudancetable<-table(abundance_other$Genus, abundance_other$Phylum) # table of genera in each phyla
#write.csv(genusabudancetable, "./16S/GenusAbundanceTable.csv")

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

# Figure 2G - bacteria beta diversity 
ps.SwabsBeta <- subset_samples(ps.Decontam, SampleType != "Gauze in PBS")
min(sample_sums(ps.SwabsBeta))# minimum sample read is 8675
median(sample_sums(ps.SwabsBeta)) #18637
table(sample_sums(ps.SwabsBeta))
tabSwab16S <- t(otu_table(ps.SwabsBeta))
ps.rarefied = rarefy_even_depth(ps.SwabsBeta, rngseed=1, sample.size=15000, replace=F) # THIS IS THE ONE TO GO WITH

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

sampledf <- data.frame(sample_data(ps.rarefied)) # Permanovas
pig16sunifrac <- phyloseq::distance(ps.rarefied, method = "unifrac")
adonis2(pig16sunifrac ~ Dorsal_Ventral, by= "margin", data = sampledf, permutations = 9999) #

sampledf2 <- data.frame(sample_data(ps.rarefied))
sampledf2 <-subset.data.frame(sampledf2, SampleType == "Swab")
sampledf2$Dorsal_Ventral <- factor(sampledf2$Dorsal_Ventral, levels = c("D", "V", "F"))
sampledf2 <- sampledf2 %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))
pig16sunifrac2 <- phyloseq::distance(ps.rarefied, method = "unifrac")
adonis2(pig16sunifrac2 ~ SkinFecal, by= "margin", data = sampledf2, permutations = 9999) #

# Figure 2H. Bacteria alpha diversity 
ps.Swabs <- subset_samples(relAll, SampleType != "Gauze in PBS")
ps.swabs.meta<-as.data.frame(sample_data(ps.Swabs))
ps.Swabs@sam_data$Timepoint <- factor(ps.Swabs@sam_data$Timepoint, levels = c("1", "2","3","4", "5"))

tabSwab <-microbiome::alpha(ps.Swabs, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
AlphaSwab<-merge(tabSwab,ps.swabs.meta, by = "row.names")
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


# Supplemental Figure 2B-I MAASLIN for bacteria 
Pig16Sgenus<-read.csv("./16S/Pig16Sgenus.csv")
input_genus_data <- as.data.frame(Pig16Sgenus)
rownames(input_genus_data) <-input_genus_data[,1]

Pig16SPhyla<-read.csv("./16S/Pig16SPhyla.csv")
input_phyla_data <- as.data.frame(Pig16SPhyla)
rownames(input_phyla_data) <-input_phyla_data[,1]

Pig16S.meta <- data.frame(sample_data(ps.Decontam))
Pig16S.meta$SampleID <- row.names(Pig16S.meta)

input_meta_file <- Pig16S.meta
input_meta_file <-subset.data.frame(input_meta_file, SampleType == "Swab")
input_meta_file$Dorsal_Ventral <- factor(input_meta_file$Dorsal_Ventral, levels = c("D", "V", "F"))
input_meta_file <- input_meta_file %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))

SkinOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Skin")
FecalOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Fecal")

# evaluating dorsal vs. ventral 
fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "./16S/maaslin/DorsalVentral_genus_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))
fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "./16S/maaslin/DorsalVentral_phyla_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

# evaluating pig 1 skin vs. pig 2 skin
fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "./16S/maaslin/PigSkinPigSkin_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))
fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "./16S/maaslin/PigSkinPigSkin_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

# evaluating pig 1 fecal vs. pig 2 fecal
fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = FecalOnly_meta_file,
  output = "./16S/maaslin/PigFecalPigFecal_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))
fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = FecalOnly_meta_file,
  output = "./16S/maaslin/PigFecalPigFecal_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))



# Figure 2.A and B
iso<-read.csv("./IsolatesPerBodySite.csv")
iso<-as.data.frame(iso)

iso$Genus<- factor(iso$Genus, levels = rev(c("Brachybacterium", "Corynebacterium", "Kocuria",  "Luteococcus", "Micrococcus", "Rothia", "Paenibacillus",  "Staphylococcus", "Viridibacillus", "Acinetobacter", "Pseudomonas")))

ggplot(iso, aes(x= "", y = NumberOfIsolates, fill = Genus))+
  geom_bar(stat="identity", position = "fill", width = 1)+
  facet_wrap(~BodySite)+
  scale_fill_manual(values =rev(c( "#56396A",  "#825C9C",  "#BEA8CA", "#6D3765","#A05195", "#CBA2C7",  "#374C80", "#5D81D9",  "#A9B7DA",  "#FF764A",  "#FFC2AD")))+
  coord_polar("y", start=0) +
  theme_void()

ggplot(iso, aes(x= BodySite, y = NumberOfIsolates))+
  geom_bar(aes(fill = Genus), position = "stack",  stat="identity")+
  scale_y_continuous(limits = c(0,30), breaks = c(0,1,5,10,15,20,25,30))+ 
  scale_fill_manual(values =rev(c( "#56396A",  "#825C9C",  "#BEA8CA", "#6D3765","#A05195", "#CBA2C7",  "#374C80", "#5D81D9",  "#A9B7DA",  "#FF764A",  "#FFC2AD")))+
  theme_light()

# Figure 2.C and D
biotable = read.csv("./bioassaytable.csv")
bioassay = as.data.frame(biotable[1:92,])
bioassay[is.na(bioassay)] = -1
df_num = as.matrix(bioassay[1:92,2:23])
rownames(df_num) = bioassay$ID

annotation_df = data.frame("Genera" = bioassay$Genera, "Phyla" = bioassay$Phyla)
rownames(annotation_df) = rownames(df_num) # name matching

column_df <- data.frame(read.csv("./PathogenGroup.csv"))
row.names(column_df) <- column_df$Pathogen

annotation_df$Genera <- factor(annotation_df$Genera, levels = c("Brachybacterium", "Corynebacterium", "Kocuria", "Luteococcus", "Micrococcus", "Rothia", "Paenibacillus", "Pseudomonas", "Staphylococcus", "Viridibacillus", "Acinetobacter"))
annotationcolors = list(
  Group = c(Bacteria_GramPositive = "#7A5195", Bacteria_GramNegative = "#d45087", Fungi = "#ffa600"),
  Phyla = c(Actinobacteriota = "#7A5195", Firmicutes = "#374C80", Proteobacteria = "#FF764A"),
  Genera = c(Brachybacterium = "#56396A", Corynebacterium = "#825C9C", Kocuria = "#D0C0D8", Luteococcus = "#6D3765", Micrococcus = "#A05195", Rothia = "#CBA2C7", Paenibacillus = "#374C80", Staphylococcus = "#5D81D9", Viridibacillus = "#A9B7DA", Acinetobacter = "#FF764A", Pseudomonas = "#FFC2AD"))

pheatmap(df_num, fontsize_row = 3, annotation_row = annotation_df, annotation_col = column_df, annotation_colors = annotationcolors, color=c("#Cccccc", "#EBF0FF", "#8DBDD9", "#357EB9", "#0C3281"))


bioasssayInScore <- biotable[1:92,c(1,24:28)]
bioasssayInS.Fun<-bioasssayInScore 
bioasssayInS.Fun$FractionalInhibition <- bioasssayInS.Fun$Fungi
bioasssayInS.Fun$InhibitedGroup <- "Fungi"

bioasssayInS.GP<-bioasssayInScore 
bioasssayInS.GP$FractionalInhibition <- bioasssayInS.GP$Gram.Positive
bioasssayInS.GP$InhibitedGroup <- "Gram Positive Bacteria"

bioasssayInS.GN<-bioasssayInScore 
bioasssayInS.GN$FractionalInhibition <- bioasssayInS.GN$Gram.Negative
bioasssayInS.GN$InhibitedGroup <- "Gram Negative Bacteria"

bioassay.InS.Final <- rbind(bioasssayInS.Fun,bioasssayInS.GP,bioasssayInS.GN)
bioassay.InS.Final$InhibitedGroup<-factor(bioassay.InS.Final$InhibitedGroup, levels = c("Fungi","Gram Positive Bacteria", "Gram Negative Bacteria"))

ggplot(bioassay.InS.Final,aes(x = InhibitedGroup, y = FractionalInhibition, color = InhibitedGroup, fill = InhibitedGroup))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~Genera, nrow = 3, ncol = 4)+
  scale_fill_manual(values = c("#ffa600", "#7A5195", "#d45087"))+
  scale_color_manual(values = c("#B87700", "#56396A", "#B52C65"))+
  theme_light()

ggplot(bioassay.InS.Final,aes(x = InhibitedGroup, y = FractionalInhibition, color = InhibitedGroup, fill = InhibitedGroup))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values = c("#ffa600", "#7A5195", "#d45087"))+
  scale_color_manual(values = c("#B87700", "#56396A", "#B52C65"))+
  theme_light()


# Figure 4.A
fig4 <- read.csv('./isolate distribution in the body sites.txt')
fig4$Body.Site <- factor(fig4$Body.Site, levels = c("Dorsal", "Ventral", "Axilla", "Inguinal", "Nares", "Oral"))
fig4$Isolate<-factor(fig4$Isolate, levels = c("LK2510", "LK2521", "LK2536", "LK2537", "LK2567", "LK2582","LK2602", "LK2640", "LK2538","LK2569",
                                               "LK2598", "LK2625", "LK2660","LK2664","LK2561",
                                               "LK2588", "LK2492", "LK2545", "LK2570"))

fig1colors <- c(colorRampPalette(c("#003F5C", "#009DE5", "#BAD7E5", "#2F4B7B", "#5B92EF", "#BFD1EF"))(10),
               colorRampPalette(c("#665191", "#DED8E9"))(3),
               colorRampPalette(c("#A05195"))(1), 
               colorRampPalette(c("#D45087"))(1),
               colorRampPalette(c("#FA5D6B", "#FACED2"))(2),
               colorRampPalette(c("#FE7C43"))(1), 
               colorRampPalette(c("#FFA601"))(1))
ggplot(fig4, aes(x=`Body.Site`)) + 
  geom_bar(aes(fill=`Isolate`)) + 
  scale_fill_manual(values = fig1colors)

# Figure 4.B
coryne <- read.table('./AntiSMASH_Results_Corynebacterium.txt', header = F, sep = '\t')
names(coryne)[names(coryne) == "V4"] <- "BGC Type"
#coryne$V1 <- factor(coryne$V1, levels = c("LK2567", "LK2569", "LK2510", "LK2582", "LK2602", "LK2538", "LK2640", "LK2536", "LK2537", "LK2521"))
coryne$V1 <- factor(coryne$V1, levels = c("LK2569", "LK2510", "LK2538", "LK2640", "LK2521", "LK2536", "LK2537", "LK2567", "LK2582", "LK2602"))
names(coryne)[names(coryne) == "V8"] <- "Body Site"
coryne$Genus <- "Corynebacterium"

kocuria <- read.table('./AntiSMASH_Results_Kocuria.txt', header = F, sep = '\t')
names(kocuria)[names(kocuria) == "V4"] <- "BGC Type"
names(kocuria)[names(kocuria) == "V8"] <- "Body Site"
kocuria$Genus <- "Kocuria"

micro <- read.table('./AntiSMASH_Results_Micrococcus.txt', header = F, sep = '\t')
names(micro)[names(micro) == "V4"] <- "BGC Type"
names(micro)[names(micro) == "V8"] <- "Body Site"
micro$Genus <- "Micrococcus"

rothia <- read.table('./AntiSMASH_Results_Rothia.txt', header = F, sep = '\t')
names(rothia)[names(rothia) == "V4"] <- "BGC Type"
names(rothia)[names(rothia) == "V8"] <- "Body Site"
rothia$Genus <- "Rothia"

paeni <- read.table('./AntiSMASH_Results_for_Paenibacillus.txt', header = F, sep = '\t')
names(paeni)[names(paeni) == "V4"] <- "BGC Type"
names(paeni)[names(paeni) == "V8"] <- "Body Site"
paeni$Genus <- "Paenibacillus"

staph <- read.table('./AntiSMASH_Results_staphylococcus.txt', header = F, sep = '\t')
names(staph)[names(staph) == "V4"] <- "BGC Type"
names(staph)[names(staph) == "V8"] <- "Body Site"
staph$Genus <-"Staphylococcus"

viridibacillus <- read.table('./AntiSMASH_Results_Viridibacillus.txt', header = F, sep = '\t')
names(viridibacillus)[names(viridibacillus) == "V4"] <- "BGC Type"
names(viridibacillus)[names(viridibacillus) == "V8"] <- "Body Site"
viridibacillus$Genus <-"Viridibacillus"

isolateBGC <-rbind(coryne, rothia, kocuria, paeni, staph, micro,viridibacillus )
isolateBGC$`BGC Type`<-factor(isolateBGC$`BGC Type`, levels = c("betalactone", "cyclic-lactone-autoinducer", "ectoine", "lanthipeptide-class-iii", "lanthipeptide-class-v", "linaridin","NI-siderophore","opine-like-metallophore", "multi-type", "multi-type (with NRPS & PKS)", "NAPAA", "NRPS", "NRPS-like", "multi-type (with NRPS)", "T1PKS", "T3PKS", "transAT-PKS", "multi-type (with PKS)", "proteusin", "RiPP-like", "RRE-containing", "terpene"))

BGCcolors <- c(colorRampPalette(c("#003F5C", "#009DE5", "#BAD7E5", "#2F4B7B", "#5B92EF", "#BFD1EF"))(6),
               colorRampPalette(c("#665191", "#DED8E9"))(2),
               "gray","black",
               colorRampPalette(c("#A05195", "#EA77DA", "#93C1D7"))(4), 
               colorRampPalette(c("#D45087","#F45C9C", "#F4D6E3"))(4),
               colorRampPalette(c("#FA5D6B"))(1),
               colorRampPalette(c("#FE7C43", "#FEE0D2"))(2), 
               colorRampPalette(c("#FFA601"))(1))

ggplot(isolateBGC, aes(x=`Body Site`)) +
  geom_bar(aes(fill=`BGC Type`), position = "fill") + 
  facet_grid(. ~ V9, scales = "free_x", space = "free") + 
  scale_fill_manual(values = BGCcolors)+
  theme_light()

# Supplemental Figure 3
ggplot(isolateBGC, aes(x=V1)) +
  geom_bar(aes(fill=`BGC Type`)) + 
  facet_grid(. ~ Genus+ V9, scales = "free_x", space = "free") + 
  scale_fill_manual(values = BGCcolors)+
  theme_light()
ggplot(isolateBGC, aes(x=V1)) +
  geom_bar(aes(fill=`BGC Type`)) + 
  facet_grid(. ~ Genus, scales = "free_x", space = "free") + 
  scale_fill_manual(values = BGCcolors)+
  theme_light()



