# ITS sequencing results for pig skin study - took swabs from the pig skin to find out what is there  

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


# ps.NAME = a phyloseq object --> this is importing the datasets and associated metadata 
ps.PigITS <- qza_to_phyloseq(features = "/Users/karindadelacruz/LK16S006/ITSAnalysisSingle/table-PigSkinITSSingle-dada2.qza",
                          taxonomy = "/Users/karindadelacruz/LK16S006/ITSAnalysisSingle/taxonomy-PigSkinITSSingle.qza",
                          metadata = "/Users/karindadelacruz/LK16S006/ITSAnalysisSingle/ITSSingle_Manufest_PigSkin.txt",
                          tree = "/Users/karindadelacruz/LK16S006/ITSAnalysisSingle/PigSkinITSSingle-rooted-tree.qza")

# create dataframe with all OTUs that are in negative sample
ps.negativeotusITS <- subset_samples(ps.PigITS, Dorsal_Ventral == "Neg")
ps.negativeotusITS <- subset_samples(ps.negativeotusITS, SampleType != "PBS no gauze")
negativeotusITS <- data.frame(ps.negativeotusITS@otu_table)
negativeotusITS <- negativeotusITS[rowSums(negativeotusITS[])>0,]
ITSOTUstoRemove <- as.character(rownames(negativeotusITS))

# Removing contaminants identified via QIIME
ps.DecontamITS <- ps.PigITS

ITSOTUDoNotRemove <- setdiff(taxa_names(ps.DecontamITS), ITSOTUstoRemove)
ps.DecontamITS <- prune_taxa(ITSOTUDoNotRemove, ps.DecontamITS)

ps.DecontamITS <- subset_samples(ps.DecontamITS, Dorsal_Ventral != "Neg")
ps.DecontamITS <- subset_samples(ps.DecontamITS,Dorsal_Ventral != "Pos") 
ps.DecontamITS# THE CLEANED FILE 
sample_names(ps.DecontamITS)
ps.DecontamITS <- subset_taxa(ps.DecontamITS, !is.na(Phylum))
ps.DecontamITS <- subset_taxa(ps.DecontamITS, !is.na(Genus))

# For maaslin
PigITS.meta <- data.frame(sample_data(ps.DecontamITS))
PigITS.meta$SampleID <- row.names(PigITS.meta)

# Making compiled files of the total and relative abundance 
PigITSout <- psmelt(ps.DecontamITS)
write.csv(PigITSout, "/Users/karindadelacruz/LK16S006/R_ITS/PigITSabundanceTotalsOut.csv")
head(PigITSout)

relAllITS <- transform(ps.DecontamITS, "compositional")
relAlloutITS <- psmelt(relAllITS)
write.csv(relAlloutITS, "/Users/karindadelacruz/LK16S006/R_ITS/PigITSRelativeOUT.csv")

# Grouping
ps.ExGroupITS <- merge_samples(ps.DecontamITS, "ExGroup", fun = mean)

Combined_ManITS <- as.data.frame(ITSSingle_CombinedManufest_PigSkin)
rownames(Combined_ManITS) <- Combined_ManITS$ExGroup

sample_data(ps.ExGroupITS) <- sample_data(Combined_ManITS)
relCombinedITS <- transform(ps.ExGroupITS, "compositional")

## RELATIVE ABUDANCE PLOTS
# Relative abundance for whole dataset
relCombinedITS@sam_data$Timepoint <- factor(relCombinedITS@sam_data$Timepoint, levels = c("1", "2","3","4", "5"))
plot_bar(relCombinedITS,  x = "Timepoint", fill = "Phylum") + facet_wrap(~Pig + Dorsal_Ventral, scales = "free_x", nrow = 2)

# Relative abundance plots for an average by pig
relCombinedITSAbundance <-psmelt(relCombinedITS)
write.csv(relCombinedITSAbundance, "/Users/karindadelacruz/LK16S006/R_ITS/RelativeAbun_OUT.csv")

#relative abundance plot for phyla >0.5% reads
relCombinedITS_other<- relCombinedITSAbundance %>% mutate(Phylum = ifelse(Abundance < 0.005, "Other", Phylum))
write.csv(relCombinedITS_other, "/Users/karindadelacruz/LK16S006/R_ITS/relCombinedITS_other.csv")

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

#relative abundance plot for genus >0.5% reads
relCombinedITSGenus <- relCombinedITSAbundance
relCombinedITSGenus<- relCombinedITSGenus %>% mutate(Genus = ifelse(Abundance < 0.005, "Other", Genus))
write.csv(relCombinedITSGenus, "/Users/karindadelacruz/LK16S006/R_ITS/PrelCombinedITSGenus.csv")

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

## ALPHA AND BETA ABUDANCE
# create alpha abundance table
relAllITS@sam_data$Timepoint <- factor(relAllITS@sam_data$Timepoint, levels = c("1", "2", "3", "4", "5"))
plot_richness(relAllITS, x = "Timepoint", measures=c("Shannon"), color = "Pig") + 
  facet_wrap("Dorsal_Ventral", scales = "free_x", ncol = 2) + 
  theme_bw()

tabSwabITS <- microbiome::alpha(relAllITS, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                  "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                  "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
write.csv(tabSwabITS, "/Users/karindadelacruz/LK16S006/R_ITS/PigITSSwabalphaTable.csv")

#Make a metadata table with the alpha diversity metrics and import it
AlphaSwabITS <- PigITSSwabalphaTable # the alpha diversity table from above with incorporated metadata and imported into R
AlphaSwabITS$Timepoint <- factor(AlphaSwabITS$Timepoint, levels = c("1", "2","3","4", "5"))

#Alpha diversity
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

# Beta diversity 
ps.SwabsBetaITS <- subset_samples(ps.DecontamITS, SampleType != "Gauze in PBS")

min(sample_sums(ps.SwabsBetaITS))# minimum sample read is 29
median(sample_sums(ps.SwabsBetaITS)) #2185

table(sample_sums(ps.SwabsBetaITS))

tabSwabITS <- t(otu_table(ps.SwabsBetaITS))

ps.rarefiedITS = rarefy_even_depth(ps.SwabsBetaITS, rngseed=1, sample.size=600, replace=F)

#Beta diversity - weighted unifrac - use this! (axis 1 and axis 3)
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

#Beta diversity - weighted unifrac - supplemental! (axis 1 and axis 2)
GP.ordITS <- ordinate(GPITS, "PCoA",  "unifrac", weighted = T) 
plot_ordination(GPITS, GP.ordITS, type="samples", axes = c(1,2), color="Dorsal_Ventral") + 
  stat_ellipse(type = "t", level = 0.8)+
  geom_point(aes(shape =Pig),size=3) + #ggtitle("Bray Curtis Beta Diversity")+
  #scale_shape_manual(values = c(8,21))+
  theme_bw(base_size = 8)+
  scale_color_manual(values = allGroupsColors,
                     name = "Swab Type") +
  scale_shape_discrete(name="Animal",
                       labels=c("Pig 1", "Pig 2")) +
  theme(legend.position = "bottom", legend.box = "vertical")

#Beta diversity - weighted unifrac - supplemental! (axis 2 and axis 3)
GP.ordITS <- ordinate(GPITS, "PCoA",  "unifrac", weighted = T) 
plot_ordination(GPITS, GP.ordITS, type="samples", axes = c(2,3), color="Dorsal_Ventral") + 
  stat_ellipse(type = "t", level = 0.8)+
  geom_point(aes(shape =Pig),size=3) + #ggtitle("Bray Curtis Beta Diversity")+
  #scale_shape_manual(values = c(8,21))+
  theme_bw(base_size = 8)+
  scale_color_manual(values = allGroupsColors,
                     name = "Swab Type") +
  scale_shape_discrete(name="Animal",
                       labels=c("Pig 1", "Pig 2")) +
  theme(legend.position = "bottom", legend.box = "vertical")

# permanova
sampledfITS <- data.frame(sample_data(ps.rarefiedITS))
pigITSunifrac <- phyloseq::distance(ps.rarefiedITS, method = "unifrac")
unifrac_clustITS <-hclust(pigITSunifrac)
plot(unifrac_clustITS)
adonis2(pigITSunifrac ~ Dorsal_Ventral, by= "margin", data = sampledfITS, permutations = 9999) #
Perm1ITS <- adonis2(pigITSunifrac ~ Dorsal_Ventral, by= "margin", data = sampledfITS, permutations = 9999) #
perm1dfITS <- as.data.frame(Perm1ITS$aov.tab)

sampledf2ITS <- data.frame(sample_data(ps.rarefiedITS))
sampledf2ITS <-subset.data.frame(sampledf2ITS, SampleType == "Swab")
sampledf2ITS$Dorsal_Ventral <- factor(sampledf2ITS$Dorsal_Ventral, levels = c("D", "V", "F"))
sampledf2ITS <- sampledf2ITS %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))
pigITSunifrac2 <- phyloseq::distance(ps.rarefiedITS, method = "unifrac")
unifrac_clust2 <-hclust(pigITSunifrac2)
plot(unifrac_clust2)
adonis2(pigITSunifrac2 ~ SkinFecal, by= "margin", data = sampledf2ITS, permutations = 9999) #
Perm2ITS <- adonis2(pigITSunifrac2 ~ SkinFecal, by= "margin", data = sampledf2ITS, permutations = 9999) #
perm2dfITS <- as.data.frame(Perm2ITS$aov.tab)

## MAASLIN
# Running Maaslin2 for determining differential abundance - similar to above but theoretically more accurate
# https://huttenhower.sph.harvard.edu/maaslin/
# input = two .tsv (or .csv) files. the first is the relative abundance table (####taxaRElAbundance.tsv).
# the second file is the .tsv of the metadata file (16S_Manufest_PigSkin.tsv)

input_OTU_data <- as.data.frame(PigITSOTUs)
rownames(input_OTU_data) <- input_OTU_data[,1]

input_genus_data <- as.data.frame(PigITSGenus)
rownames(input_genus_data) <-input_genus_data[,1]

input_phyla_data <- as.data.frame(PigITSPhyla)
rownames(input_phyla_data) <-input_phyla_data[,1]

input_meta_file <- PigITS.meta
input_meta_file <-subset.data.frame(input_meta_file, SampleType == "Swab")
input_meta_file$Dorsal_Ventral <- factor(input_meta_file$Dorsal_Ventral, levels = c("D", "V", "F"))
input_meta_file <- input_meta_file %>% mutate(SkinFecal = ifelse(Dorsal_Ventral == "F" , "Fecal", "Skin"))

SkinOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Skin")
FecalOnly_meta_file <-subset.data.frame(input_meta_file, SkinFecal == "Fecal")

# evaluating fecal vs. skin 
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = input_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/SkinFecal_OTU_maaslin",
  fixed_effects = c("SkinFecal"),
  reference = c("Fecal", "Skin"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = input_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/SkinFecal_genus_maaslin",
  fixed_effects = c("SkinFecal"),
  reference = c("Fecal", "Skin"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = input_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/SkinFecal_phyla_maaslin",
  fixed_effects = c("SkinFecal"),
  reference = c("Fecal", "Skin"))

# evaluating dorsal vs. ventral 
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/DorsalVentral_OTU_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/DorsalVentral_genus_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/DorsalVentral_phyla_maaslin",
  fixed_effects = c("Dorsal_Ventral"),
  reference = c("Dorsal", "Ventral"))

# evaluating pig 1 skin vs. pig 2 skin
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/PigSkinPigSkin_OTU_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/PigSkinPigSkin_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = SkinOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/PigSkinPigSkin_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

# evaluating pig 1 fecal vs. pig 2 fecal
fit_data = Maaslin2(
  input_data = input_OTU_data,
  input_metadata = FecalOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/PigFecalPigFecal_OTU_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_genus_data,
  input_metadata = FecalOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/PigFecalPigFecal_genus_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

fit_data = Maaslin2(
  input_data = input_phyla_data,
  input_metadata = FecalOnly_meta_file,
  output = "/Users/karindadelacruz/LK16S006/R_ITS/maaslin/PigFecalPigFecal_phyla_maaslin",
  fixed_effects = c("Pig"),
  reference = c("WM1", "WM2"))

## SPIEC-EASI - skin only
# load packages necessary
install.packages("devtools") 
library(devtools)
install_github("zdk123/SpiecEasi", force=TRUE)
library(SpiecEasi)
library(igraph)

install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

install.packages("taxonomizr")
library(taxonomizr)

# subset the samples to only include skin data only
ps.skinonlyITS <- subset_samples(relCombinedITS, Dorsal_Ventral %in% c("D","V"))

# extract the otu and taxonomy tables from the final phyloseq object
ps.skinonlyITS #finalized phyloseq object
otusITS=otu_table(ps.skinonlyITS)
taxaITS=tax_table(ps.skinonlyITS)
otusITS <- t(otusITS)

# filter the OTU matrix to remove zeros, but keep the sum to not change the overall sample counts
filterobjITS=filterTaxonMatrix(otusITS, minocc=1, keepSum = TRUE, return.filtered.indices = TRUE) # otusITS must have taxa as rows and samples as columns
otusITS.f=filterobjITS$mat

# introduce dummy taxonomic information for the last row of the Taxa table
taxaITS.f=taxaITS[setdiff(1:nrow(taxaITS),filterobjITS$filtered.indices),]
#dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
#taxaITS.f=rbind(taxaITS.f,dummyTaxonomy)
rownames(taxaITS.f)[nrow(taxaITS.f)]="0"
rownames(otusITS.f)[nrow(otusITS.f)]="0"

# assemble a new phyloseq object with the filtered OTU and taxonomy table
updatedotusITS=otu_table(otusITS.f, taxa_are_rows = TRUE)
view(updatedotusITS)
updatedtaxaITS=tax_table(taxaITS.f)
view(updatedtaxaITS)
phyloseqobjITS.f=phyloseq(updatedotusITS, updatedtaxaITS)

# make a phyloseq object that contains enough data but not too much
ps.spieceasiITS <- phyloseq_filter_prevalence(phyloseqobjITS.f, prev.trh = 0.1, abund.trh = 10, threshold_condition = "OR")

# apply SPIEC-EASI - p level 0.01
spiecITS.out=spiec.easi(ps.spieceasiITS, method="mb", icov.select.params=list(rep.num=50, thresh=0.05))
# convert the output of SPIEC-EASI into a network and plot it with colors on a class level
spiecITS.graph=adj2igraph(getRefit(spiecITS.out), vertex.attr=list(name=taxa_names(ps.spieceasiITS)))
plot_network(spiecITS.graph, ps.spieceasiITS, type="taxa", color="Phylum", label=NULL, point_size = 4, line_weight=0.5)

# count positive and negative edges to extract the regression coefficients from the output
betaMatITS=as.matrix(symBeta(getOptBeta(spiecITS.out)))

# obtain the number of edges in the network from the matrix of regression coefficients
positiveITS=length(betaMatITS[betaMatITS>0])/2
negativeITS=length(betaMatITS[betaMatITS<0])/2
totalITS=length(betaMatITS[betaMatITS!=0])/2

# cluster SPIEC-EASI network and list the taxa present in each cluster
clustersITS=cluster_fast_greedy(spiecITS.graph)
clusterOneIndicesITS=which(clustersITS$membership==1)
clusterOneOtusITS=clustersITS$names[clusterOneIndicesITS]
clusterTwoIndicesITS=which(clustersITS$membership==2)
clusterTwoOtusITS=clustersITS$names[clusterTwoIndicesITS]

# summarize and sort the representatives of the different classes
sort(table(getTaxonomy(clusterOneOtusITS,taxaITS.f,useRownames = TRUE)),decreasing = TRUE)
sort(table(getTaxonomy(clusterTwoOtusITS,taxaITS.f,useRownames = TRUE)),decreasing = TRUE)

# export while preserving the OTU identifiers
write.graph(spiecITS.graph,file=file.path("/Users/karindadelacruz/LK16S006/spieceasinetworks/spieceasiITSnetwork.txt"),format="ncol")

# export the taxonomic information
write.table(taxaITS.f,file=file.path("/Users/karindadelacruz/LK16S006/spieceasinetworks/spieceasiITS.txt"),sep="\t", quote=FALSE)

# open in cytoscape
