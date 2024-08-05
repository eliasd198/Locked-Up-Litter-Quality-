###############################################
#### Figure 5 ####
###############################################

# Set Working Directory
setwd()

# Load required packages 

devtools::install_github('smin95/smplot2', force = TRUE)
library(smplot2)

library(plyr)
library(dplyr)
library(tidyr)
library(car)
library(phyloseq)
library(ggplot2)
library(data.table)
library(olsrr)
library(relaimpo)
library(emmeans)
library(microbiome)
library(tibble)
library(vegan)
library(ggordiplots)
library(forcats)
library(patchwork)
library(ggpubr)
library(ggpattern)
library(see)

#### Import Data into Phyloseq ####
fun = read.csv("ITS_ASV.csv", header = TRUE, row.names=1, check.names=F)
fun.tax = read.csv("ITS_TAXA.csv", header = TRUE,row.names=1, check.names=F)
fun.env <- read.csv("16S_ITS_META.csv",header=T, row.names=2, check.names = F)

# Recode soil and litter levels for plotting
fun.env$Mineral <- fct_recode(fun.env$Mineral, "No Minerals"="Unamended")
fun.env$Litter <- fct_recode(fun.env$Litter, "Low Quality" = "Winter Wheat", "High Quality" = "White Clover")

# Create phyloseq object
fun.tax <- as.matrix(fun.tax)

OTU_its <- otu_table(fun, taxa_are_rows = T)   
TAX_its <- tax_table(fun.tax)
ENV_its <- sample_data(fun.env)

fun_its <- merge_phyloseq(OTU_its,TAX_its,ENV_its)
fun_its
ntaxa(fun_its)
nsamples(fun_its)

# Add reference sequences and rename ASV's to something more user-friendly
# for downstream analysis

sequences <- Biostrings::DNAStringSet(taxa_names(fun_its))
names(sequences) <- taxa_names(fun_its)
fun_its <- merge_phyloseq(fun_its, sequences)
fun_its

taxa_names(fun_its) <- paste0("ASV", seq(ntaxa(fun_its)))

##########################
#### Data Exploration ####
##########################

# Firstly remove any non Fungal sequences
fun_its <- subset_taxa(fun_its, Kingdom=="Fungi")
ntaxa(fun_its)

# Summarize sequencing depths with a Histogram
sdt = data.table(as(sample_data(fun_its), "data.frame"),
                 TotalReads = sample_sums(fun_its), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

# Read Counts
min(sample_sums(fun_its))
max(sample_sums(fun_its))
mean(sample_sums(fun_its))

# Sequencing depth by plot and land use
# Check that there is no systematic bias
pSeqDepth + facet_wrap(~Mineral)
pSeqDepth + facet_wrap(~Litter)

# Look at distribution of Taxa
tdt = data.table(tax_table(fun_its),
                 TotalCounts = taxa_sums(fun_its),
                 OTU = taxa_names(fun_its))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts") + xlim(0, 25000) + ylim(0,1000)

# How many zero sum OTUs?
tdt[(TotalCounts <= 0), .N]

# Remove these zero sum OTU's
fun_its <- prune_taxa(taxa_sums(fun_its) > 0, fun_its)
ntaxa(fun_its)

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Define the plot - Can use to assess filtering thresholds 
# fore rare taxa if required
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() + xlim(0,100) +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum

############################
#### Data Normalization ####
############################

# Rarefy to an even sequencing depth for ordinations of community composition. 

# Calculate total reads to inform at what depth to rarefy
sdt = data.table(as(sample_data(fun_its), "data.frame"),
                 TotalReads = sample_sums(fun_its), keep.rownames = TRUE)

min(sdt$TotalReads)
head(sort(sdt$TotalReads), 50)

# Min reads is 0. A few other samples with very low reads. Rarefy to 10984
# Set a seed so this step is reproducible
set.seed(3)

fun.r <- rarefy_even_depth(fun_its, sample.size = 10984,
                                 replace = TRUE, trimOTUs = TRUE,
                                 verbose = TRUE, rngseed = 3)

nsamples(fun.r) # Lost 6 samples due to low read count
ntaxa(fun.r) 

################
# Data Prep ####
################

# Subset soil and litter saprotrophs, transform to relative abundance and summarize data by treatment ####
SS_FUN_bar <- fun.r %>% 
  subset_taxa(primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="litter_saprotroph"|
                primary_lifestyle=="wood_saprotroph"|primary_lifestyle=="unspecified_saprotroph") %>%
  tax_glom(taxrank = "primary_lifestyle") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>% # Melt to long format
  arrange(primary_lifestyle)

# Gives you the sum of % contributions of fungal saprotroph groups by land use
DT <- data.table(SS_FUN_bar, key = c("primary_lifestyle", "Sample", "Timepoint", "Mineral", "Litter"))
DT <- DT[, sum(Abundance), by = key(DT)]
setnames(DT, "V1", "Abundance")

# Remove NA's and summarize data prior to plotting. 
# Make sure plyr is loaded before dplyr
DT_summary <- DT %>% drop_na() %>%
  group_by(primary_lifestyle, Timepoint, Mineral, Litter) %>%
  dplyr::summarise(
    mean = mean(Abundance),
    sd = sd(Abundance), 
    n = n(),
    se = sd / sqrt(n)
  )

# Write to CSV for future use if necessary
write.csv(DT, "Fungi_Relative_abundance.csv")

############################
# Figure 5 & statistics ####
############################

# Subset T126
T2 <- subset_samples(fun.r, Timepoint=="126 Days")

#Extract OTU table in correct format for vegan package 
f.otu <- (otu_table(T2))
f.otu <- otu_table(T2)
if (taxa_are_rows(f.otu)) {
  f.otu <- t(f.otu)
}
f.mat <- as(f.otu, "matrix")
f.df <- as.data.frame(f.mat)

# Create distance matrix using bray-curtis dissimlarities
f.dist <- vegdist(f.df, method="bray")

# Run NMDS
f.nmds <- metaMDS(f.dist, k = 2, trymax = 500, autotransform = T)

# Extract metadata
f.meta <- meta(T2)
f.meta$Mineral <- fct_relevel(f.meta$Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")

# Plot NMDS Ordination
f.plot <- gg_ordiplot(f.nmds, groups = f.meta$Mineral, ellipse = F, spiders = T, label = T, pt.size = 2)
f.plot2 <- gg_ordiplot(f.nmds, groups = f.meta$Litter, ellipse = F, spiders = T, label = T, pt.size = 2)

Figure_5A <- ggplot(data = f.plot$df_ord, aes(x = x, y = y)) +
  scale_shape_manual(name = "Plot Name", values=c(3,4,23,24)) +
  scale_color_manual(name = "Plot Name", values=c("#000000", "#0072B2", "#D55E00"), guide="none") +
  geom_point(data = f.plot$df_ord, aes(color = f.meta$Litter, shape = f.meta$Mineral), size = 3) + 
  xlab("NMDS1") + 
  ylab("NMDS2") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_segment(data = f.plot2$df_spiders,
               aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group), 
               show.legend = FALSE) +
  geom_label(data = f.plot2$df_mean.ord, aes(x = x, y = y, label = Group), show.legend = F, colour = "black", fontface = "bold", size = 5) +
  theme_bw() + theme(
    text=element_text(color = "black", size=21),
    axis.text = element_text(color = "black", size=21),
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
    plot.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.8,0.8)) +
  guides(size = "none", color = "none")

Figure_5A

# Test for differences between Mineral and Litter treatments using PERMANOVA (Table S6)
adonis2(f.dist ~ Mineral*Litter, data = f.meta)

#adonis2(formula = f.dist ~ Mineral * Litter, data = f.meta)
#Df SumOfSqs      R2       F Pr(>F)    
#Mineral         3   0.8260 0.13190  5.7418  0.001 ***
#  Litter          2   2.6699 0.42636 27.8393  0.001 ***
#  Mineral:Litter  6   0.6563 0.10481  2.2812  0.001 ***
#  Residual       44   2.1099 0.33693                   
#Total          55   6.2622 1.00000    

# Plot a Bar graph of the relative abundance soil and litter saprotrophs between treatments

# Add error bars to stacked bar plot in correct position
DT_summary$y_pos = NA
DT_summary$y_pos[DT_summary$primary_lifestyle == "soil_saprotroph"] = DT_summary$mean[DT_summary$primary_lifestyle == "soil_saprotroph"]
DT_summary$y_pos[DT_summary$primary_lifestyle == "litter_saprotroph"] = DT_summary$mean[DT_summary$primary_lifestyle == "soil_saprotroph"] + 
DT_summary$mean[DT_summary$primary_lifestyle == "litter_saprotroph"]
str(DT_summary)

# Overlay data points to stacked bar plot in correct positions

DT_filter <- DT %>% filter(primary_lifestyle=="litter_saprotroph"|primary_lifestyle=="soil_saprotroph")
DT_filter_126 <- DT_filter %>% filter(Timepoint=="126 Days")

DT_filter_126$y_pos = NA
DT_filter_126$y_pos[DT_filter_126$primary_lifestyle == "soil_saprotroph"] = DT_filter_126$Abundance[DT_filter_126$primary_lifestyle == "soil_saprotroph"]
DT_filter_126$y_pos[DT_filter_126$primary_lifestyle == "litter_saprotroph"] = DT_filter_126$Abundance[DT_filter_126$primary_lifestyle == "soil_saprotroph"] + 
  DT_filter_126$Abundance[DT_filter_126$primary_lifestyle == "litter_saprotroph"]

DT_filter_126$primary_lifestyle = fct_recode(DT_filter_126$primary_lifestyle, Litter = "litter_saprotroph", Soil = "soil_saprotroph")

Figure_5B <- DT_summary %>%  group_by(Timepoint, Mineral, Litter) %>%
  filter(Timepoint=="126 Days", primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="litter_saprotroph") %>%
  mutate(primary_lifestyle = fct_recode(primary_lifestyle, Litter = "litter_saprotroph", Soil = "soil_saprotroph"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality")) %>%
  ggplot(mapping=aes(x = Litter, y = mean, fill=Mineral, alpha=primary_lifestyle)) + 
  geom_bar(aes(), color="black", stat = "identity") +
  geom_errorbar(aes(ymin = y_pos-se, ymax = y_pos+se, width=.3)) +
  geom_jitter(DT_filter_126,  mapping = aes(x = Litter,y = y_pos)) +
  facet_grid(~Mineral) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  labs(alpha='', size = 12)  +
  ylab("Relative Abundance \n") +
  ylim(0,1) +
  guides(fill="none") +
  theme_bw() + theme(
    legend.position="right",
    text=element_text(color = "black", size=21),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 14),
    strip.background = element_rect(fill = NA))

Figure_5B

# Plot correlations between litter-C respired and relative abundance of soil saprotrophs

# Extract metadata
metadata <- meta(fun.r)
metadata <- rownames_to_column(metadata)

# Merge fungal abundance with metadata
DT_summary2 <- inner_join(metadata, DT, by=c('rowname'='Sample', 'Timepoint'='Timepoint', 'Mineral'='Mineral', 'Litter'='Litter'))
print(DT_summary2)

Figure_5C <- DT_summary2 %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="126 Days", !Litter=="No Litter", primary_lifestyle=="soil_saprotroph") %>%
  ggplot(DT_summary2, mapping=aes(x=Abundance, y=respired_litter_c_perc, linetype=Mineral)) + 
  geom_point(aes(shape=Litter, fill=Mineral), size=2) + 
  facet_grid(~Mineral, scales="fixed") + 
  stat_cor(label.x = 0.15, label.y = 1, p.accuracy = 0.001, r.accuracy = 0.1) +
  sm_statCorr(aes(group=Mineral), color="black", label_x = c(0.15), label_y = c(1), show_text=FALSE) +
  xlab("\n Relative abundance of soil saprotrophs") +
  ylab("% of litter-derived C respired") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_linetype_manual(values = c("No Minerals" = "solid", "Kaolinite" = "solid", "Goethite" = "blank", "Montmorillonite" = "solid")) +
  scale_shape_manual(values = c(24, 21)) +
  guides(fill = "none", linetype = "none") +
  theme_bw() + theme(
    legend.position="right",
    legend.text = element_text(size=16),
    text=element_text(color = "black", size=21),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    strip.background.y = element_rect(fill = NA))

Figure_5C

# Plot correlations between MAOM formation efficiency and relative abundance of soil saprotrophs
Figure_5D <- DT_summary2 %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="126 Days", !Litter=="No Litter", primary_lifestyle=="soil_saprotroph") %>%
  ggplot(DT_summary2, mapping=aes(x=Abundance, y=maom_formation_efficiency, fill=Mineral, linetype=Mineral)) + 
  geom_point(aes(shape=Litter), size=2) + 
  facet_grid(~Mineral, scales="fixed") + 
  stat_cor(label.x = 0.2, label.y = 0.2, p.accuracy = 0.001, r.accuracy = 0.1) +
  sm_statCorr(aes(group=Litter), color="black", label_x = c(0.2), label_y = c(0.2), show_text=FALSE) +
  xlab("\n Relative abundance of soil saprotrophs") +
  ylab("MAOM Formation Efficiency") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  #scale_fill_manual(values=c('white', "white")) +
  scale_linetype_manual(values = c("No Minerals" = "blank", "Kaolinite" = "solid", "Goethite" = "blank", "Montmorillonite" = "solid")) +
  scale_shape_manual(values = c(24, 21)) +
  #scale_alpha_manual(values=c("Winter Wheat"=0, "White Clover"=1)) +
  guides(fill = "none", linetype = "none") +
  theme_bw() + theme(
    legend.position="right",
    #legend.box.background = element_rect(color="black", linewidth=1),
    legend.text = element_text(size=16),
    #legend.box.margin = margin(2, 2, 2, 2),
    text=element_text(color = "black", size=21),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    #panel.grid = element_blank(),
    legend.title = element_blank(),
    strip.text.x = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    strip.background.y = element_rect(fill = NA))

Figure_5D

# Plot fungal community composition ordinations at T126 alongside saprotroph 
# relative abundance bar plots and correlations with litter-C respired 
# and MAOM formation efficiency
Figure_5 <- (Figure_5A | (Figure_5B/Figure_5C/Figure_5D + 
                              plot_layout(guides = 'auto'))) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1.25, 1)) &
  theme(
        legend.direction = 'vertical',
        legend.box = 'vertical',
        plot.tag = element_text(size = 16, face = 'bold'))

Figure_5

ggsave("Figure_5.tiff", Figure_5, height = 30, width = 50, units = "cm", compression="lzw")

#############################
# Figure S2 & Statistics ####
#############################

T1 <- subset_samples(fun.r, Timepoint=="15 Days")

#Extract OTU table in correct format for vegan package 
b.otu_T1 <- (otu_table(T1))
b.otu_T1 <- otu_table(T1)
if (taxa_are_rows(b.otu_T1)) {
  b.otu_T1 <- t(b.otu_T1)
}
b.mat_T1 <- as(b.otu_T1, "matrix")
b.df_T1 <- as.data.frame(b.mat_T1)

# Create distance matrix using bray-curtis dissimlarities
b.dist_T1 <- vegdist(b.df_T1, method="bray")

# Run NMDS
b.nmds_T1 <- metaMDS(b.dist_T1, k = 2, trymax = 500, autotransform = T)

# Extract metadata
b.meta_T1 <- meta(T1)
b.meta_T1$Mineral <- fct_relevel(b.meta_T1$Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")

# Plot NMDS Ordination
b.plot_T1 <- gg_ordiplot(b.nmds_T1, groups = b.meta_T1$Mineral, ellipse = F, spiders = T, label = T, pt.size = 2)
b.plot2_T1 <- gg_ordiplot(b.nmds_T1, groups = b.meta_T1$Litter, ellipse = F, spiders = T, label = T, pt.size = 2)

Figure_S2A <- ggplot(data = b.plot_T1$df_ord, aes(x = x, y = y)) +
  scale_shape_manual(name = "Plot Name", values=c(3,4,23,24)) +
  scale_color_manual(name = "Plot Name", values=c("#000000", "#0072B2", "#D55E00")) +
  geom_point(data = b.plot_T1$df_ord, aes(color = b.meta_T1$Litter, shape = b.meta_T1$Mineral), size = 3) + 
  xlab("NMDS1") + 
  ylab("NMDS2") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_segment(data = b.plot2_T1$df_spiders,
               aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group), 
               show.legend = FALSE) +
  geom_label(data = b.plot2_T1$df_mean.ord, aes(x = x, y = y, label = Group), show.legend = F, colour = "black", fontface = "bold", size = 3) +
  theme_bw() + theme(
    text=element_text(color = "black", size=21),
    axis.text = element_text(color = "black", size=21),
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
    plot.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.85,0.85)) +
  guides(size = "none", color = "none")

Figure_S2A

# Test for Differences between Mineral and Litter treatments
adonis2(b.dist_T1 ~ Mineral*Litter, data = b.meta_T1)

#adonis2(formula = b.dist_T1 ~ Mineral * Litter, data = b.meta_T1)
#Df SumOfSqs      R2       F Pr(>F)    
#Mineral         3   0.5221 0.08840  4.8535  0.001 ***
#  Litter          2   3.0630 0.51865 42.7135  0.001 ***
#  Mineral:Litter  6   0.7071 0.11973  3.2869  0.001 ***
#  Residual       45   1.6135 0.27321                   
#Total          56   5.9056 1.00000                   


## Plot a Bar chart of relative abundance of soil and litter saprotrophs at T15
DT_filter_15 <- DT_filter %>% filter(Timepoint=="15 Days")

DT_filter_15$y_pos = NA
DT_filter_15$y_pos[DT_filter_15$primary_lifestyle == "soil_saprotroph"] = DT_filter_15$Abundance[DT_filter_15$primary_lifestyle == "soil_saprotroph"]
DT_filter_15$y_pos[DT_filter_15$primary_lifestyle == "litter_saprotroph"] = DT_filter_15$Abundance[DT_filter_15$primary_lifestyle == "soil_saprotroph"] + 
  DT_filter_15$Abundance[DT_filter_15$primary_lifestyle == "litter_saprotroph"]

DT_filter_15$primary_lifestyle = fct_recode(DT_filter_15$primary_lifestyle, Litter = "litter_saprotroph", Soil = "soil_saprotroph")

Figure_S2B <- DT_summary %>%  group_by(Timepoint, Mineral, Litter) %>%
  filter(Timepoint=="15 Days", primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="litter_saprotroph") %>%
  mutate(primary_lifestyle = fct_recode(primary_lifestyle, Litter = "litter_saprotroph", Soil = "soil_saprotroph"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality")) %>%
  ggplot(mapping=aes(x = Litter, y = mean, fill=Mineral, alpha=primary_lifestyle)) + 
  geom_bar(aes(), color="black", stat = "identity") +
  geom_errorbar(aes(ymin = y_pos-se, ymax = y_pos+se, width=.3)) +
  geom_jitter(DT_filter_15,  mapping = aes(x = Litter,y = y_pos)) +
  facet_grid(~Mineral) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  labs(alpha='', size = 12)  +
  ylab("Relative Abundance \n") +
  ylim(0,1) +
  guides(fill="none") +
  theme_bw() + theme(
    legend.position="right",
    text=element_text(color = "black", size=21),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 14),
    strip.background = element_rect(fill = NA))

Figure_S2B

Figure_S2C <- DT_summary2 %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="15 Days", !Litter=="No Litter", primary_lifestyle=="soil_saprotroph") %>%
  ggplot(DT_summary2, mapping=aes(x=Abundance, y=respired_litter_c_perc)) + 
  geom_point(aes(shape=Litter, fill=Mineral), size=2) + 
  facet_grid(~Mineral, scales="fixed") + 
  sm_statCorr(aes(group=Mineral), color="black", label_x = c(0.15), label_y = c(1)) +
  xlab("\n Relative abundance of soil saprotrophs") +
  ylab("% of litter-derived C respired") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_linetype_manual(values = c("No Minerals" = "solid", "Kaolinite" = "solid", "Goethite" = "dashed", "Montmorillonite" = "solid")) +
  scale_shape_manual(values = c(24, 21)) +
  guides(fill="none", linetype = "none") +
  theme_bw() + theme(
    legend.position="right",
    legend.text = element_text(size=16),
    text=element_text(color = "black", size=21),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size=14),
    strip.text = element_text(face = "bold", hjust = 0),
    strip.background.x = element_rect(fill = NA),
    strip.background.y = element_rect(fill = NA))

Figure_S2C

# Figure S2: Plot fungal community composition at T15 alongside relative abundance barplots 
# of soil and litter saprotrophs and correlations with litter-C respired

Figure_S2 <- Figure_S2A + Figure_S2B / Figure_S2C + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1.25, 1)) &
  theme(
    legend.direction = 'vertical',
    legend.box = 'vertical',
    plot.tag = element_text(size = 16, face = 'bold'))

Figure_S2

ggsave("Figure_S2.tiff", Figure_S2, height = 25, width = 50, units = "cm", compression="lzw")

# Figure S3: Correlation between litter C respired at T15 and soil saprotrophs ####

Figure_S3 <- DT_summary2 %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="15 Days", !Litter=="No Litter", primary_lifestyle=="litter_saprotroph") %>%
  ggplot(DT_summary2, mapping=aes(x=Abundance, y=respired_litter_c_perc)) + 
  geom_point(aes(fill=Mineral), shape=21, size=2) + 
  facet_grid(~Litter, scales="fixed") + 
  sm_statCorr(aes(group=Mineral), color="black", label_x = c(0.15), label_y = c(1)) +
  xlab("\n Relative abundance of litter saprotrophs") +
  ylab("% of litter-derived C Respired") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_linetype_manual(values = c("No Minerals" = "solid", "Kaolinite" = "solid", "Goethite" = "dashed", "Montmorillonite" = "solid")) +
  scale_shape_manual(values = c(24, 21)) +
  guides(linetype = "none") +
  theme_bw() + theme(
    legend.position="right",
    legend.text = element_text(size=16),
    text=element_text(color = "black", size=21),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size=16),
    strip.text = element_text(face = "bold", hjust = 0),
    strip.background.x = element_rect(fill = NA),
    strip.background.y = element_rect(fill = NA))

Figure_S3

ggsave("Figure_S3.tiff", sap_soilresp_15, height = 15, width = 30, units = "cm", compression="lzw")

# Figure S4: Relative abundance of soil saprotrophic genera at T15 and T126 ####

soil_taxa <- fun.r %>%
  subset_taxa(primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="litter_saprotroph"|
                primary_lifestyle=="unspecified_saprotroph"| primary_lifestyle=="wood_saprotroph") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  subset_taxa(primary_lifestyle =="soil_saprotroph") %>%      
  psmelt() %>% # Melt to long format
  arrange(primary_lifestyle)

# Gives you the sum of % contributions of Phyla by land use
DT_soil <- data.table(soil_taxa, key = c("Genus", "primary_lifestyle", "Sample", "Timepoint", "Mineral", "Litter"))
DT_soil <- DT_soil[, sum(Abundance), by = key(DT_soil)]
setnames(DT_soil, "V1", "Abundance")

# Remove NA's and summarize data prior to plotting. 
# Make sure plyr is loaded before dplyr
DT_soil_summary <- DT_soil %>% drop_na() %>%
  group_by(Genus, primary_lifestyle, Timepoint, Mineral, Litter) %>%
  dplyr::summarise(
    mean = mean(Abundance),
    sd = sd(Abundance), 
    n = n(),
    se = sd / sqrt(n)
  )

DT_soil_summary$Genus <- as.factor(DT_soil_summary$Genus)
DT_soil_summary$Genus <- as.character(DT_soil_summary$Genus) # convert to character
DT_soil_summary$Genus[DT_soil_summary$mean < 0.02] <- "< 2% Abundance"

# Plot the Bar graph
Figure_S4 <- DT_soil_summary %>%  group_by(Timepoint, Mineral, Litter) %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days")) %>%
  ggplot(mapping=aes(x = Litter, y = mean, fill=Genus, color=Genus)) + 
  geom_bar(aes(), stat = "identity") +
  facet_grid(Timepoint~Mineral) +
  scale_fill_okabeito(palette = "black_first") +
  scale_color_okabeito(palette = "black_first") +
  labs(alpha='', size = 12)  +
  ylab("Relative Abundance \n") +
  ylim(0,1) +
  theme_bw() + theme(
    legend.position="right",
    text=element_text(color = "black", size=21),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 14),
    strip.background = element_rect(fill = NA))

Figure_S4

# Run statistics on relative abundance of soil saprotrophs at T15 and T126 (Table S6)

# Function to summarize ANOVA assumptions
mod_summary <- function(x){
  ols_plot_resid_fit(x)
  ols_plot_resid_qq(x)
  ols_plot_resid_hist(x)
  ols_plot_resid_box(x)
  ols_test_normality(x)
} 

DT_soil_T15 <- DT %>%
  filter(primary_lifestyle =="soil_saprotroph", 
         Timepoint=="15 Days")

mod_soil_15 <- lm(Abundance ~ Mineral*Litter, data=DT_soil_T15)
mod_summary(mod_soil_15)  

# Model summary and statistics #
summary(mod_soil_15)
Anova(mod_soil_15)

#Anova Table (Type II tests)

#Response: Abundance
#Sum Sq Df F value    Pr(>F)    
#Mineral        0.000905  3  0.1802 0.9092567    
#Litter         0.097052  2 28.9817 8.164e-09 ***
#  Mineral:Litter 0.051351  6  5.1115 0.0004471 ***
#  Residuals      0.075346 45                      

#Assess relative importance of predictors #
calc.relimp(mod_soil_15, rela=FALSE)

DT_soil_T126 <- DT %>%
  filter(primary_lifestyle =="soil_saprotroph", 
         Timepoint=="126 Days")

mod_soil_126 <- lm(Abundance ~ Mineral*Litter, data=DT_soil_T126)
mod_summary(mod_soil_126)  

# Model summary and statistics #
summary(mod_soil_126)
Anova(mod_soil_126)

#Anova Table (Type II tests)

#Response: Abundance
#Sum Sq Df F value    Pr(>F)    
#Mineral        0.07879  3  6.8653  0.000681 ***
#  Litter         0.32702  2 42.7431 4.865e-11 ***
#  Mineral:Litter 0.05820  6  2.5355  0.033969 *  
#  Residuals      0.16832 44                    

#Assess relative importance of predictors #
calc.relimp(mod_soil_126, rela=FALSE)


# Figure S5: relative abundance of litter saprotrophic genera at T15 and T126 ####

litter_taxa <- fun.r %>%
  subset_taxa(primary_lifestyle=="soil_saprotroph"|primary_lifestyle=="litter_saprotroph"|
                primary_lifestyle=="unspecified_saprotroph"| primary_lifestyle=="wood_saprotroph") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  subset_taxa(primary_lifestyle =="litter_saprotroph") %>%      
  psmelt() %>% # Melt to long format
  arrange(primary_lifestyle)

# Gives you the sum of % contributions of Phyla by land use

DT_litter <- data.table(litter_taxa, key = c("Genus", "primary_lifestyle", "Sample", "Timepoint", "Mineral", "Litter"))
DT_litter <- DT_litter[, sum(Abundance), by = key(DT_litter)]
setnames(DT_litter, "V1", "Abundance")

# Remove NA's and summarize data prior to plotting. 
# Make sure plyr is loaded before dplyr
DT_litter_summary <- DT_litter %>% drop_na() %>%
  group_by(Genus, primary_lifestyle, Timepoint, Mineral, Litter) %>%
  dplyr::summarise(
    mean = mean(Abundance),
    sd = sd(Abundance), 
    n = n(),
    se = sd / sqrt(n)
  )

DT_litter_summary$Genus <- as.factor(DT_litter_summary$Genus)
DT_litter_summary$Genus <- as.character(DT_litter_summary$Genus) # convert to character
DT_litter_summary$Genus[DT_litter_summary$mean < 0.02] <- "< 2% Abundance"

# Plot the Bar graph
Figure_S5 <- DT_litter_summary %>%  group_by(Timepoint, Mineral, Litter) %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days")) %>%
  ggplot(mapping=aes(x = Litter, y = mean, fill=Genus, color=Genus)) + 
  geom_bar(aes(), stat = "identity") +
  facet_grid(Timepoint~Mineral) +
  scale_fill_okabeito(palette = "black_first") +
  scale_color_okabeito(palette = "black_first") +
  labs(alpha='', size = 12)  +
  ylab("Relative Abundance \n") +
  ylim(0,1) +
  theme_bw() + theme(
    legend.position="right",
    text=element_text(color = "black", size=21),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 14),
    strip.background = element_rect(fill = NA))

Figure_S5

# Run statistics on relative abundance of litter saprotrophs at T15 and T126 (Table S6)

DT_litter_T15 <- DT %>%
  filter(primary_lifestyle =="litter_saprotroph", 
         Timepoint=="15 Days")

mod_litter_15 <- lm(Abundance ~ Mineral*Litter, data=DT_litter_T15)
mod_summary(mod_litter_15)  

# Model summary and statistics #
summary(mod_litter_15)
Anova(mod_litter_15)

#Anova Table (Type II tests)

#Response: Abundance
#Sum Sq Df F value    Pr(>F)    
#Mineral        0.102326  3 16.3654 2.481e-07 ***
#  Litter         0.173547  2 41.6340 5.816e-11 ***
#  Mineral:Litter 0.016961  6  1.3563    0.2529    
#Residuals      0.093789 45                   

#Assess relative importance of predictors #
calc.relimp(mod_litter_15, rela=FALSE)

#Proportion of variance explained by model: 76.25%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.2696202
#Litter         0.4499549
#Mineral:Litter 0.0429463

DT_litter_T126 <- DT %>%
  filter(primary_lifestyle =="litter_saprotroph", 
         Timepoint=="126 Days")

mod_litter_126 <- lm(Abundance ~ Mineral*Litter, data=DT_litter_T126)
mod_summary(mod_litter_126)  

# Model summary and statistics #
summary(mod_litter_126)
Anova(mod_litter_126)

#Anova Table (Type II tests)

#Response: Abundance
#Sum Sq Df F value    Pr(>F)    
#Mineral        0.049342  3  4.3254  0.009322 ** 
#  Litter         0.242276  2 31.8575 2.792e-09 ***
#  Mineral:Litter 0.017116  6  0.7502  0.612555    
#Residuals      0.167310 44                 

#Assess relative importance of predictors #
calc.relimp(mod_litter_126, rela=FALSE)

#Proportion of variance explained by model: 65%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.10530176
#Litter         0.50889982
#Mineral:Litter 0.03580411

# Save plots
ggsave("Figure_S4.tiff", Figure_S4, height = 20, width = 30, units = "cm", compression="lzw")
ggsave("Figure_S5.tiff", Figure_S5, height = 20, width = 30, units = "cm", compression="lzw")


# Richness and diversity Estimates on raw samples ####
# Estimate Fungal Richness for Table 2 ####
f.rich <- estimate_richness(fun.r)
f.rich <- data.frame(f.rich, sample_data(fun.r))
mean(f.rich$Observed)
sd(f.rich$Observed)

f.rich$Mineral <- as.factor(f.rich$Mineral)
f.rich$Litter <- as.factor(f.rich$Litter)

# Summarize the data by soil and litter treatments
my.summary = function(x) list(mean = mean(x), SE = parameters::standard_error(x))

f.rich <- data.table(f.rich)

f.rich_summary <- f.rich[, sapply(.SD, my.summary), by = c("Mineral", "Litter", "Timepoint"), .SDcols = c('Observed',"Shannon")]
setnames(f.rich_summary, 4:7, c("Observed", "Observed_SE", "Shannon", "Shannon_SE"))

f.rich_summary

write.csv(f.rich_summary, "Fungal_Richness_Summary.csv")

# Model effects of soil and litter on fungal richness at T15 and T126 (Table S4)
f.rich_T1 <- subset(f.rich, Timepoint=="15 Days")
f.rich_T2 <- subset(f.rich, Timepoint=="126 Days")

mod1 <- lm(Observed ~ Mineral*Litter, data=f.rich_T1)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod1)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod1)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod1)
# 4. Create a Boxplot II
ols_plot_resid_box(mod1)
# 5. Perform a Normality Test
ols_test_normality(mod1)

# Model summary and statistics #
summary(mod1)
Anova(mod1)

#Anova Table (Type II tests)

#Response: Observed
#Sum Sq Df F value Pr(>F)    
#Mineral         14919  3  2.1744 0.1042    
#Litter         440112  2 96.2199 <2e-16 ***
#  Mineral:Litter  19515  6  1.4221 0.2273    
#Residuals      102916 45   

#Assess relative importance of predictors #
calc.relimp(mod1, rela=FALSE)

#Proportion of variance explained by model: 82.07%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.02302476
#Litter         0.76370388
#Mineral:Litter 0.03399444

# Post-Hoc comparisons using emmeans #
emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")

mod2 <- lm(Observed ~ Mineral*Litter, data=f.rich_T2)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod2)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod2)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod2)
# 4. Create a Boxplot II
ols_plot_resid_box(mod2)
# 5. Perform a Normality Test
ols_test_normality(mod2)

# Model summary and statistics #
summary(mod2)
Anova(mod2)

#Anova Table (Type II tests)

#Response: Observed
#Sum Sq Df  F value    Pr(>F)    
#Mineral          4217  3   1.6514  0.191298    
#Litter         242337  2 142.3343 < 2.2e-16 ***
#  Mineral:Litter  19306  6   3.7797  0.004029 ** 
#  Residuals       37457 44           

#Assess relative importance of predictors #
calc.relimp(mod2, rela=FALSE)

#Proportion of variance explained by model: 87.71%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.01631615
#Litter         0.79747258
#Mineral:Litter 0.06333285

# Post-Hoc comparisons using emmeans #
emmeans(mod2, specs = pairwise ~ Litter|Mineral, type="response")

# END #