#' ---
#' title: "Data analysis on mycobiota data of USEFUL project"
#' author: Francesco Vitali
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' ```{r global_options, echo = FALSE, include = FALSE}
#' options(width = 9999)
#' knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
#'                      cache = FALSE, tidy = FALSE, size = "small",
#'                      fig.height = 10, fig.width = 15,fig.align = "center")
#' ```

#' # INTRODUCTION
#' 
#' In this script we will perform classical alpha and beta diversity analysis on the results of mycobiota sequencing 
#' of Taleggio cheese. We will perform analysis on the R1 and R2 separately (i.e. using picking results from non merged reads)
#' then comparing results (i.e. Mantel test or procsutes rotation). Finally, we will explore sthe effect of some experimental metadata
#' 
#' 

#+ echo=FALSE, message = FALSE

# Load libraries
library(phyloseq)
library(ggplot2)
library(reshape2)
library(ape)
library(data.table)
library(ggsci)
library(ggpubr)
library(picante)
library(microbiome)
library(cowplot)
library(metagenomeSeq)
library(pals)
library(dplyr)
library(tidyr)
library(knitr)
library(pairwiseAdonis)
library(stringr)
library(patchwork)
library(ggrepel)

####################
## Loading data ####
####################

# Loading R1 data
metadata_Taleggio <- read.table(file = "./picking_results/taleggio_metadata.csv", header = T, sep = "\t", row.names = 1)
otutable_Taleggio <- read.table(file = "./picking_results/R1/otutable_R1.txt", header = T, sep = "", row.names = 1)
taxonomy_Taleggio <- read.table(file = "./picking_results/R1/taxa_R1_Rform.csv", sep = "\t", row.names = 1)
colnames(taxonomy_Taleggio) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# change rownames metadata to match colnames otutable
row.names(metadata_Taleggio) <- paste("X", row.names(metadata_Taleggio), sep = "")
summary(colnames(otutable_Taleggio) %in% row.names(metadata_Taleggio))
# mounting phyloseq object
tax_table_R1 <- tax_table(as.matrix(taxonomy_Taleggio))
otu_table_R1 <- otu_table(as.matrix(otutable_Taleggio), taxa_are_rows = T)
meta_table_R1 <- sample_data(metadata_Taleggio)
phyloseq_obj_initial_R1 <- merge_phyloseq(tax_table_R1,otu_table_R1,meta_table_R1)
phyloseq_obj_initial_R1


# Loading R2 data
metadata_Taleggio <- read.table(file = "./picking_results/taleggio_metadata.csv", header = T, sep = "\t", row.names = 1)
otutable_Taleggio <- read.table(file = "./picking_results/R2/otutable_R2.txt", header = T, sep = "", row.names = 1)
taxonomy_Taleggio <- read.table(file = "./picking_results/R2/taxa_R2_Rform.csv", sep = "\t", row.names = 1)
colnames(taxonomy_Taleggio) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# change rownames metadata to match colnames otutable
row.names(metadata_Taleggio) <- paste("X", row.names(metadata_Taleggio), sep = "")
row.names(metadata_Taleggio) <- gsub("R1", "R2", row.names(metadata_Taleggio)) # R converts the "-" in the otu table as "."
summary(colnames(otutable_Taleggio) %in% row.names(metadata_Taleggio))
# mounting phyloseq object
tax_table_R2 <- tax_table(as.matrix(taxonomy_Taleggio))
otu_table_R2 <- otu_table(as.matrix(otutable_Taleggio), taxa_are_rows = T)
meta_table_R2 <- sample_data(metadata_Taleggio)
phyloseq_obj_initial_R2 <- merge_phyloseq(tax_table_R2,otu_table_R2,meta_table_R2)
phyloseq_obj_initial_R2

######################
## Data transform ####
######################

# Creating rarefied object to verify if rarefaction affects alpha diversity
phyloseq_obj_raref_R1 <- rarefy_even_depth(phyloseq_obj_initial_R1, 
                                             rngseed=1234, 
                                             sample.size=round(0.99*min(sample_sums(phyloseq_obj_initial_R1))), 
                                             replace=F)

phyloseq_obj_raref_R2 <- rarefy_even_depth(phyloseq_obj_initial_R2, 
                                           rngseed=1234, 
                                           sample.size=round(0.99*min(sample_sums(phyloseq_obj_initial_R2))), 
                                           replace=F)

# Creating CSS object for beta diversity
doubleton <- genefilter_sample(phyloseq_obj_initial_R1, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj_initial_R1) 
# transforming 
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
phylo_obj_css_R1 <- phyloseq_obj_initial_R1
otu_table(phylo_obj_css_R1) <- otu_table(data.CSS, taxa_are_rows = T)

doubleton <- genefilter_sample(phyloseq_obj_initial_R2, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj_initial_R2) 
# transforming 
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
phylo_obj_css_R2 <- phyloseq_obj_initial_R2
otu_table(phylo_obj_css_R2) <- otu_table(data.CSS, taxa_are_rows = T)


######################
## Beta Diversity ####
######################

# We need to better undestand the samples set-up. It should be a repeated sampling, so pairwise sampling.

## NMDS ordination
library(broom)
set.seed(12387)
df <- as(sample_data(phylo_obj_css_R1), "data.frame")
distance.R1 <-  phyloseq::distance(phylo_obj_css_R1, "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df, Sample)
adonis_R1 <- adonis2(distance.R1 ~ Treatment + Locality + as.factor(Days_maturation) + Maturation_type, data = df, permutations = perm) # should also consider individual 

ordination_R1<- ordinate(phylo_obj_css_R1, "NMDS", distance = "bray")

# PLOT1
p_ordination_R1 <- plot_ordination(physeq = phylo_obj_css_R1, ordination_R1, axes = c(1,2))
p <- p_ordination_R1 + geom_point(size = 5, stroke = 0.8, aes(fill = Maturation_type, shape = Locality_short))
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) #+ facet_grid(. ~ Time)
p_1 <- p  + scale_shape_manual(values = c(21,22,23,24)) + 
  xlab("NMDS1") + ylab("NMDS2") + labs( fill = "Maturation Type", shape = "Locality") + theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23)))

annotations <- data.frame(
  xpos = c(-Inf,Inf,Inf),
  ypos =  c(Inf,Inf, Inf),
  annotateText = c(sprintf("NMDS Stress= %s", round(ordination_R1$stress, digits = 4))),
  hjustvar = c(-0.2) ,
  vjustvar = c(2)) #<- adjust

p_R1 <- p_1 + 
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) + ggtitle("Fungi NMDS Ordination")
p_R1

ggsave("NMDS1_colMaturation_shapeLocality.pdf",device = "pdf")

# PLOT2
p_ordination_R1 <- plot_ordination(physeq = phylo_obj_css_R1, ordination_R1, axes = c(1,2))
p <- p_ordination_R1 + geom_point(size = 5, stroke = 0.8, aes(fill = Locality_short, shape = Treatment))
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) #+ facet_grid(. ~ Time)
p_1 <- p  + scale_shape_manual(values = c(21,22,23,24)) + 
  xlab("NMDS1") + ylab("NMDS2") + labs( fill = "Locality", shape = "Treatment") + theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23)))

annotations <- data.frame(
  xpos = c(-Inf,Inf,Inf),
  ypos =  c(Inf,Inf, Inf),
  annotateText = c(sprintf("NMDS Stress= %s", round(ordination_R1$stress, digits = 4))),
  hjustvar = c(-0.2) ,
  vjustvar = c(2)) #<- adjust

p_R1 <- p_1 + 
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) + ggtitle("Fungi NMDS Ordination")
p_R1

ggsave("NMDS1_colLocality_shapeTreatment.pdf",device = "pdf")

## First variable seems to e VTA+CAL and BAL+PAN. Create a variable "group" and then subset to see effect of other variables

sample_data(phylo_obj_css_R1)$Group <- sample_data(phylo_obj_css_R1)$Locality_short
sample_data(phylo_obj_css_R1)$Group[sample_data(phylo_obj_css_R1)$Group == "BAL"] <- "Group2"
sample_data(phylo_obj_css_R1)$Group[sample_data(phylo_obj_css_R1)$Group == "PAN"] <- "Group2"
sample_data(phylo_obj_css_R1)$Group[sample_data(phylo_obj_css_R1)$Group == "CAL"] <- "Group1"
sample_data(phylo_obj_css_R1)$Group[sample_data(phylo_obj_css_R1)$Group == "VTA"] <- "Group1"

R1_group1 <- subset_samples(phylo_obj_css_R1, Group == "Group1")  %>% prune_taxa(taxa_sums(.) > 0.00001, .)
R1_group2 <- subset_samples(phylo_obj_css_R1, Group == "Group2")  %>% prune_taxa(taxa_sums(.) > 0.00001, .)

# Group 1
library(broom)
set.seed(12387)
df <- as(sample_data(R1_group1), "data.frame")
distance.R1 <-  phyloseq::distance(R1_group1, "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df, Sample)
df$Days_maturation_fct <- factor(x = df$Days_maturation, levels = c("35","60"))
adonis_R1 <- adonis2(distance.R1 ~ Treatment + Locality + Days_maturation + Maturation_type, data = df, permutations = perm) # should also consider individual 

# none is significant

ordination_R1_G1 <- ordinate(R1_group1, "NMDS", distance = "bray")

# PLOT1
p_ordination_R1_G1 <- plot_ordination(physeq = R1_group1, ordination_R1_G1, axes = c(1,2))
p <- p_ordination_R1_G1 + geom_point(size = 5, stroke = 0.8, aes(fill = Maturation_type, shape = Locality))
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) #+ facet_grid(. ~ Time)
p_1 <- p  + scale_shape_manual(values = c(21,22,23,24)) + 
  xlab("NMDS1") + ylab("NMDS2") + labs( fill = "Maturation Type", shape = "Locality") + theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23)))

annotations <- data.frame(
  xpos = c(-Inf,Inf,Inf),
  ypos =  c(Inf,Inf, Inf),
  annotateText = c(sprintf("NMDS Stress= %s", round(ordination_R1$stress, digits = 4))),
  hjustvar = c(-0.2) ,
  vjustvar = c(2)) #<- adjust

p_R1 <- p_1 + 
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) + ggtitle("Fungi NMDS Ordination - CAL + VTA sample")
p_R1 + facet_wrap(.~Treatment)

ggsave("NMDS1_CALeVTA.pdf",device = "pdf")


# Group 2
library(broom)
set.seed(12387)
df <- as(sample_data(R1_group2), "data.frame")
distance.R1 <-  phyloseq::distance(R1_group2, "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df, Sample)
df$Days_maturation_fct <- factor(x = df$Days_maturation, levels = c("35","60"))
adonis_R1 <- adonis2(distance.R1 ~ Treatment + Locality + Days_maturation + Maturation_type, data = df, permutations = perm) # should also consider individual 

ordination_R1_G2 <- ordinate(R1_group2, "NMDS", distance = "bray")

# PLOT1
p_ordination_R1_G2 <- plot_ordination(physeq = R1_group2, ordination_R1_G2, axes = c(1,2))
p <- p_ordination_R1_G2 + geom_point(size = 5, stroke = 0.8, aes(fill = Maturation_type, shape = Locality))
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) #+ facet_grid(. ~ Time)
p_1 <- p  + scale_shape_manual(values = c(21,22,23,24)) + 
  xlab("NMDS1") + ylab("NMDS2") + labs( fill = "Maturation Type", shape = "Locality") + theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23)))

annotations <- data.frame(
  xpos = c(-Inf,Inf,Inf),
  ypos =  c(Inf,Inf, Inf),
  annotateText = c(sprintf("NMDS Stress= %s", round(ordination_R1$stress, digits = 4))),
  hjustvar = c(-0.2) ,
  vjustvar = c(2)) #<- adjust

p_R1 <- p_1 + 
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) + ggtitle("Fungi NMDS Ordination - BAL + PAN sample")
p_R1 + facet_wrap(.~Treatment)

ggsave("NMDS1_BALePAN.pdf",device = "pdf")



######################
## Alpha Diversity ####
######################

# With the same logic as beta diversity, I directly explore data divided in Group1 and Group2

sample_data(phyloseq_obj_raref_R1)$Group <- sample_data(phyloseq_obj_raref_R1)$Locality_short
sample_data(phyloseq_obj_raref_R1)$Group[sample_data(phyloseq_obj_raref_R1)$Group == "BAL"] <- "Group2"
sample_data(phyloseq_obj_raref_R1)$Group[sample_data(phyloseq_obj_raref_R1)$Group == "PAN"] <- "Group2"
sample_data(phyloseq_obj_raref_R1)$Group[sample_data(phyloseq_obj_raref_R1)$Group == "CAL"] <- "Group1"
sample_data(phyloseq_obj_raref_R1)$Group[sample_data(phyloseq_obj_raref_R1)$Group == "VTA"] <- "Group1"

R1_group1_raref <- subset_samples(phyloseq_obj_raref_R1, Group == "Group1")  %>% prune_taxa(taxa_sums(.) > 0.00001, .)
R1_group2_raref <- subset_samples(phyloseq_obj_raref_R1, Group == "Group2")  %>% prune_taxa(taxa_sums(.) > 0.00001, .)

### GROUP 1

alphadiv <- microbiome::alpha(R1_group1_raref)
df <- as(sample_data(R1_group1_raref), "data.frame")
df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson
df$dom <- alphadiv$dominance_gini

p_rich <- ggboxplot(data = df, x = "Days_maturation", y = "rich", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_rich <- p_rich + ggtitle(label =   "Richness (n° of ASVs)")

p_eve <- ggboxplot(data = df, x = "Days_maturation", y = "eve", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_eve <- p_eve + ggtitle(label =   "Pielou's Evenness")

p_shan <- ggboxplot(data = df, x = "Days_maturation", y = "shan", add = "jitter", 
                   color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_shan <- p_shan + ggtitle(label =   "Shannon's Index")

p_rich + p_eve + p_shan

ggsave("Richness_CALeVAL.pdf",device = "pdf")

## 

### GROUP 2

alphadiv <- microbiome::alpha(R1_group2_raref)
df <- as(sample_data(R1_group2_raref), "data.frame")
df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson
df$dom <- alphadiv$dominance_gini

p_rich <- ggboxplot(data = df, x = "Days_maturation", y = "rich", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_rich <- p_rich + ggtitle(label =   "Richness (n° of ASVs)")

p_eve <- ggboxplot(data = df, x = "Days_maturation", y = "eve", add = "jitter", 
                   color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_eve <- p_eve + ggtitle(label =   "Pielou's Evenness")

p_shan <- ggboxplot(data = df, x = "Days_maturation", y = "shan", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_shan <- p_shan + ggtitle(label =   "Shannon's Index")

p_rich + p_eve + p_shan

ggsave("Richness_BALePAN.pdf",device = "pdf")


### GROUP 1

alphadiv <- microbiome::alpha(R1_group1_raref)
df <- as(sample_data(R1_group1_raref), "data.frame")
df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson
df$dom <- alphadiv$dominance_gini

p_rich <- ggboxplot(data = df, x = "Treatment", y = "rich", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_rich <- p_rich + ggtitle(label =   "Richness (n° of ASVs)")

p_eve <- ggboxplot(data = df, x = "Treatment", y = "eve", add = "jitter", 
                   color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_eve <- p_eve + ggtitle(label =   "Pielou's Evenness")

p_shan <- ggboxplot(data = df, x = "Treatment", y = "shan", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_shan <- p_shan + ggtitle(label =   "Shannon's Index")

p_rich + p_eve + p_shan

ggsave("Richness_CALeVAL_treat.pdf",device = "pdf")

## 

### GROUP 2

alphadiv <- microbiome::alpha(R1_group2_raref)
df <- as(sample_data(R1_group2_raref), "data.frame")
df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson
df$dom <- alphadiv$dominance_gini

p_rich <- ggboxplot(data = df, x = "Treatment", y = "rich", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_rich <- p_rich + ggtitle(label =   "Richness (n° of ASVs)")

p_eve <- ggboxplot(data = df, x = "Treatment", y = "eve", add = "jitter", 
                   color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_eve <- p_eve + ggtitle(label =   "Pielou's Evenness")

p_shan <- ggboxplot(data = df, x = "Treatment", y = "shan", add = "jitter", 
                    color = "Maturation_type", palette ="jama", rotate = F, dot.size = 4, facet.by = "Locality") + theme_cleveland() + stat_compare_means()
p_shan <- p_shan + ggtitle(label =   "Shannon's Index")

p_rich + p_eve + p_shan

ggsave("Richness_BALePAN_treat.pdf",device = "pdf")


#####################
## Composition  #####
#####################


otu_table_collapsed <- tax_glom(phyloseq_obj_raref_R1, taxrank="Species")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 0.1, TRUE)
tb <- psmelt(qiime_file_proportional_oneperc) %>%
  as_tibble
colnames(tb)
tb0 <- tb %>%
  group_by(OTU, Locality_short, Label, Maturation_type, Species) %>%
  add_tally() %>%
  dplyr::summarize(Mean = mean(Abundance), SD = sd(Abundance), n = mean(n)) %>%
  arrange(desc(Mean)) %>%
  ungroup

p_phylum <- ggplot(aes(x=Label, y=Mean, fill = Species), data = tb0)
p_phylum <- p_phylum + 
  geom_bar(stat="identity", position = "stack", color = "black") + 
  theme_linedraw() + 
  scale_fill_manual(values = as.character(polychrome(18))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab("Abundance (%)") + 
  xlab("") + ylim(c(0,100)) 
p_phylum + ggtitle("Fungi Species composition (> 1%)")

ggsave("Species_composition.pdf",device = "pdf")



