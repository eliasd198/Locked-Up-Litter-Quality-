# 16S ####

# Set Working Directory for 16S
setwd()

# Load required packages
library(ggordiplots)
library(plyr)
library(dplyr)
library(phyloseq)
library(microbiome)
library(ggpubr)
library(ggplot2)
library(forcats)
library(olsrr)
library(car)
library(relaimpo)
library(emmeans)
library(multcomp)
library(multcompView)
library(data.table)
library(FD)
library(patchwork)
library(ggnewscale)
library(vegan)
library(devtools)
library(Biostrings)

#devtools::install_github('smin95/smplot2', force = TRUE)
library(smplot2)

###################################
#### Import Data into Phyloseq ####
###################################

bac = read.csv("16S_ASV.csv", header = TRUE,row.names=1, check.names = F)
bac.tax = read.csv("16S_TAX.csv", header = TRUE,row.names=1)
bac.env <- read.csv("16S_ITS_META.csv",header=T, row.names = 1)
bac.traits <- read.csv("16S_COPY_NUMBER.csv",header=T, row.names = 1)

# Recode soil and litter levels for plotting
bac.env$Mineral <- fct_recode(bac.env$Mineral, "No Minerals"="Unamended")
bac.env$Litter <- fct_recode(bac.env$Litter, "Low Quality" = "Winter Wheat", "High Quality" = "White Clover")

# Create Phyloseq object
bac.tax <- as.matrix(bac.tax, rownames.force = 1)

OTU_16s <- otu_table(bac, taxa_are_rows = T)   
TAX_16s <- tax_table(bac.tax)
ENV_16s <- sample_data(bac.env)

bac_16S <- merge_phyloseq(OTU_16s,TAX_16s,ENV_16s)
bac_16S
ntaxa(bac_16S)
nsamples(bac_16S)

# Add reference sequences and rename ASV's to something more user-friendly
# for downstream analysis
#sequences <- Biostrings::DNAStringSet(taxa_names(bac_16S))
#names(sequences) <- taxa_names(bac_16S)
#bac_16S <- merge_phyloseq(bac_16S, sequences)
#bac_16S

#taxa_names(bac_16S) <- paste0("ASV", seq(ntaxa(bac_16S)))

##########################
#### Data Exploration ####
##########################

# Firstly remove any non bacterial sequences that may have been amplified: Chromista and Plantae)
bac_16S <- subset_taxa(bac_16S, Kingdom=="Bacteria")

# And remove Thermus spike used for normalisation
bac_16S <- subset_taxa(bac_16S, !Genus=="Thermus")

# Summarize sequencing depths with a Histogram
sdt = data.table(as(sample_data(bac_16S), "data.frame"),
                 TotalReads = sample_sums(bac_16S), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

# Read Counts
min(sample_sums(bac_16S))
max(sample_sums(bac_16S))
mean(sample_sums(bac_16S))

# Sequencing depth by treatments
# Check that there is no systematic bias by treatment

pSeqDepth + facet_wrap(~Mineral)
pSeqDepth + facet_wrap(~Litter)

# Look at distribution of Taxa

tdt = data.table(tax_table(bac_16S),
                 TotalCounts = taxa_sums(bac_16S),
                 OTU = taxa_names(bac_16S))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts") + xlim(0, 25000) + ylim(0,250)

# How many zero sum OTUs?
tdt[(TotalCounts <= 0), .N]

# Remove these zero sum OTU's
bac_16S <- prune_taxa(taxa_sums(bac_16S) > 0, bac_16S)

# taxa cumulative sum. 
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
sdt = data.table(as(sample_data(bac_16S), "data.frame"),
                 TotalReads = sample_sums(bac_16S), keep.rownames = TRUE)

min(sdt$TotalReads)
head(sort(sdt$TotalReads), 10)

# Min reads is 129. V. low so a bad sample. Rarefy to 2nd lowest depth (2888)
# Set a seed so this step is reproducible
set.seed(3)

bac.r <- rarefy_even_depth(bac_16S, sample.size = 2888,
                           replace = TRUE, trimOTUs = TRUE,
                           verbose = TRUE, rngseed = 3)
nsamples(bac.r)
ntaxa(bac.r)

# Transform community data to relative abundance to calculate abundance weighted 16s copy number  

# Make sure we use functions from correct package
transform <- microbiome::transform
bac_16S_RA <- transform(bac_16S, "compositional")

# Extract corresponding ASV table and metadata and Tax
spc <- as.data.frame(otu_table(bac_16S_RA, "matrix"))
metadata <- meta(bac_16S_RA)
tax <- as.data.frame(tax_table(bac_16S_RA))

# Extract predicted copy number per ASV
ASV <- rownames(tax)
copy_no <- as.data.frame(tax[,7])
names(copy_no)[names(copy_no) == "tax[, 7]"] <- "copy_no"
rownames(copy_no) <- ASV
copy_no$copy_no <- as.numeric(copy_no$copy_no)

# Calculate abundance weighted mean 16S rRNA gene copy number per sample
# This gives more weight to more abundant ASV's and vica versa
copy_no_weighted <- FD::functcomp(copy_no, t(spc))

# Add abundance weighted mean 16S rRNA gene copy number to metadata 
metadata <- merge(copy_no_weighted,metadata, by="row.names")
metadata$Timepoint <- as.factor(metadata$Timepoint)
metadata$Mineral <- as.factor(metadata$Mineral)
metadata$Litter <- as.factor(metadata$Litter)
str(metadata)

# Summarize data for plotting. Make sure plyr is loaded before dplyr
bac_env_summary <- metadata %>% 
  group_by(Timepoint, Mineral, Litter) %>%
  dplyr::summarise(
    mean = mean(copy_no),
    sd = sd(copy_no), 
    n = n(),
    se = sd / sqrt(n)
  )

# Save summary data for inclusion in Table 2
write.csv(bac_env_summary, "16S_Copy_Number.csv")

############################
# Figure 4 & statistics ####
############################

# Subset T126
bac_T2 <- subset_samples(bac.r, Timepoint=="126 Days")

#Extract OTU table in correct format for vegan package 
b.otu_T2 <- (otu_table(bac_T2))
b.otu_T2 <- otu_table(bac_T2)
if (taxa_are_rows(b.otu_T2)) {
  b.otu_T2 <- t(b.otu_T2)
}
b.mat_T2 <- as(b.otu_T2, "matrix")
b.df_T2 <- as.data.frame(b.mat_T2)

# Create distance matrix using bray-curtis dissimlarities
b.dist_T2 <- vegdist(b.df_T2, method="bray")

# Run NMDS
b.nmds_T2 <- metaMDS(b.dist_T2, k = 2, trymax = 500, autotransform = T)

# Extract metadata
b.meta_T2 <- meta(bac_T2)

# Plot NMDS Ordination
b.plot_T2 <- gg_ordiplot(b.nmds_T2, groups = b.meta_T2$Mineral, ellipse = F, spiders = T, label = T, pt.size = 2)
b.plot2_T2 <- gg_ordiplot(b.nmds_T2, groups = b.meta_T2$Litter, ellipse = F, spiders = T, label = T, pt.size = 2)

Figure_4A <- ggplot(data = b.plot_T2$df_ord, aes(x = x, y = y)) +
  scale_shape_manual(name = "Plot Name", values=c(3,4,23,24)) +
  scale_color_manual(name = "Plot Name", values=c("#000000", "#0072B2", "#D55E00"), guide="none") +
  geom_point(data = b.plot_T2$df_ord, aes(color = b.meta_T2$Litter, shape = b.meta_T2$Mineral), size = 3) + 
  xlab("NMDS1") + 
  ylab("NMDS2") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_segment(data = b.plot2_T2$df_spiders,
               aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group), 
               show.legend = FALSE) +
  geom_label(data = b.plot2_T2$df_mean.ord, aes(x = x, y = y, label = Group), show.legend = F, colour = "black", fontface = "bold", size = 5) +
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
  guides(size = FALSE)

Figure_4A

# Test for differences between Mineral and Litter treatments using PERMANOVA (Table S5)
adonis2(b.dist_T2 ~ Litter*Mineral, data = b.meta_T2)

#adonis2(formula = b.dist_T2 ~ Soil * Litter, data = b.meta_T2)
#Df SumOfSqs      R2      F Pr(>F)    
#Soil         3   0.9230 0.13016 3.8308  0.001 ***
#  Litter       2   1.4695 0.20724 9.1489  0.001 ***
#  Soil:Litter  6   0.8434 0.11895 1.7504  0.001 ***
#  Residual    48   3.8548 0.54365                  
#Total       59   7.0907 1.00000    

# Model the effect of Mineral and Litter treatments on copy number #
# This allows us to add lowercase letters on plots to denote significant differences

# Subset T126
metadata_126 <- subset(metadata, Timepoint=="126 Days")

# T126 modelling
mod2 <- lm((copy_no) ~ Mineral*Litter, data=metadata_126)

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

# Model summary and statistics (Table S5) #
summary(mod2)
Anova(mod2)

#Anova Table (Type II tests)

#Response: (copy_no)
#Sum Sq Df F value    Pr(>F)    
#Soil        0.23057  3  25.649 4.759e-10 ***
#  Litter      0.31788  2  53.044 6.972e-13 ***
#  Soil:Litter 0.20181  6  11.225 8.523e-08 ***
#  Residuals   0.14383 48          

#Assess relative importance of predictors #
calc.relimp(mod2, rela=FALSE)

#Relative importance metrics: 
  
#  lmg
#Soil        0.2578819
#Litter      0.3555389
#Soil:Litter 0.2257141

# Post-Hoc comparisons using emmeans #
emm2 = emmeans(mod2, specs = pairwise ~ Litter|Mineral, type="response")
emm2

# Assign letters to denote significant differences at p<0.05 #
cld2 <- cld(emm2, Letters=letters, sort=T)

# Subset columns of interest
cld2_sub <- cld2[,c(1,2,8)] 

# Bind with summary data for plotting
bac_env_summary_126 <- subset(bac_env_summary, Timepoint=="126 Days")
bac_env_summary_126 <- inner_join(cld2_sub, bac_env_summary_126, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(bac_env_summary_126)

# Plot bar chart of 16s rRNA operon copy number

Figure_4B <- bac_env_summary_126 %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  ggplot(bac_env_summary, mapping=aes(x=Litter, y=mean, fill=Mineral)) + 
  geom_bar(aes(alpha=Litter), stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
  geom_jitter(metadata_126, mapping = aes(x=Litter, y=copy_no)) +
  geom_text(aes(label = .group, y = mean + se), hjust = 0.65, vjust = -0.5) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"), guide="none") +
  facet_grid(Timepoint~Mineral) +
  ylim(0,3) +
  ylab("Abundance weighted mean \n rRNA gene copy number") +
  theme_bw() + theme(
    legend.position="none",
    text=element_text(color = "black", size=21),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    strip.text.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 14),
    strip.background = element_rect(fill = NA))

Figure_4B

# Plot relationship between litter C respired and copy number for T126

Figure_4C <- metadata %>% 
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="126 Days", !Litter=="No Litter") %>%
  ggplot(mapping=aes(x=copy_no, y=respired_litter_c_perc, fill=Mineral, linetype=Mineral)) + 
  geom_point(aes(shape=Litter), size=2) + 
  facet_grid(~Mineral, scales="fixed") + 
  stat_cor(label.y = 1.5, p.accuracy = 0.01, r.accuracy = 0.1) +
  sm_statCorr(aes(group=Mineral), color="black", label_y = 1.5, show_text=FALSE) +
  xlab("\n 16S rRNA Gene Copy Number") +
  ylab(" % of litter-derived C respired") +
  scale_linetype_manual(values = c("Unamended" = "dashed", "Kaolinite" = "solid", "Goethite" = "solid", "Montmorillonite" = "solid")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_shape_manual(values = c(24, 21)) +
  scale_x_continuous(breaks=c(seq(2,3,0.2))) +
  guides(fill = FALSE, linetype=FALSE) +
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

Figure_4C

# Plot Relationship between maom formation efficiency and copy number for T126

Figure_4D <- metadata %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="126 Days", !Litter=="No Litter") %>%
  ggplot(mapping=aes(x=copy_no, y=maom_formation_efficiency, fill=Mineral, linetype=Mineral)) + 
  geom_point(aes(shape=Litter), size=2) + 
  facet_grid(~Mineral, scales="fixed") + 
  sm_statCorr(aes(group=Mineral), color="black", label_y = 0.2) +
  xlab("\n 16S rRNA Gene Copy Number") +
  ylab("MAOM Formation Efficiency") +
  scale_linetype_manual(values = c("Unamended" = "blank", "Kaolinite" = "solid", "Goethite" = "solid", "Montmorillonite" = "solid")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_shape_manual(values = c(24, 21)) +
  scale_x_continuous(breaks=c(seq(2,3,0.2))) +
  guides(fill = FALSE, linetype=FALSE) +
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

Figure_4D

# Plot ordination of bacterial community composition at T126 alongside 16S copy number bar graphs
# and correlations with litter-C respired and MAOM formation efficiency

Figure_4 <- (Figure_4A | (Figure_4B/Figure_4C/Figure_4D)) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1.25, 1)) &
  theme(
        legend.direction = 'vertical',
        legend.box = 'vertical',
        plot.tag = element_text(size = 16, face = 'bold'))

Figure_4

ggsave("Figure_4.tiff", Figure_4, height = 30, width = 50, units = "cm", compression="lzw")

#############################
# Figure S1 & statistics ####
#############################

# Subset T15
bac_T1 <- subset_samples(bac.r, Timepoint=="15 Days")

#Extract OTU table in correct format for vegan package 
b.otu_T1 <- (otu_table(bac_T1))
b.otu_T1 <- otu_table(bac_T1)
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
b.meta_T1 <- microbiome::meta(bac_T1)

# Plot NMDS Ordination
b.plot_T1 <- gg_ordiplot(b.nmds_T1, groups = b.meta_T1$Mineral, ellipse = F, spiders = T, label = T, pt.size = 2)
b.plot2_T1 <- gg_ordiplot(b.nmds_T1, groups = b.meta_T1$Litter, ellipse = F, spiders = T, label = T, pt.size = 2)

Figure_S1A <- ggplot(data = b.plot_T1$df_ord, aes(x = x, y = y)) +
  scale_shape_manual(name = "Plot Name", values=c(3,4,23,24)) +
  scale_color_manual(name = "Plot Name", values=c("#000000", "#0072B2", "#D55E00"), guide="none") +
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
  guides(size = FALSE)

Figure_S1A

# Test for difference in bacteria communities between Mineral and Litter treatments using PERMANOVA (Table S5)
adonis2(b.dist_T1 ~ Mineral*Litter, data = b.meta_T1)

#adonis2(formula = b.dist_T1 ~ Soil * Litter, data = b.meta_T1)
#Df SumOfSqs      R2       F Pr(>F)    
#Soil         3   0.9391 0.10699  3.4585  0.001 ***
#  Litter       2   2.7506 0.31339 15.1953  0.001 ***
#  Soil:Litter  6   0.9238 0.10525  1.7011  0.002 ** 
#  Residual    46   4.1634 0.47436                   
#Total       57   8.7769 1.00000   

# Subset T15
metadata_15 <- subset(metadata, Timepoint=="15 Days")

# T15 modelling
mod1 <- lm((copy_no) ~ Mineral*Litter, data=metadata_15)

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

# Model summary and statistics (Table S5) #
summary(mod1)
Anova(mod1)

#Anova Table (Type II tests)

#Response: (copy_no)
#Sum Sq Df F value    Pr(>F)    
#Soil        1.2510  3 10.6674 1.823e-05 ***
#  Litter      7.6718  2 98.1267 < 2.2e-16 ***
#  Soil:Litter 1.1548  6  4.9234 0.0005581 ***
#  Residuals   1.8373 47    

#Assess relative importance of predictors #
calc.relimp(mod1, rela=FALSE)

#Relative importance metrics: 

#  lmg
#Soil        0.10690469
#Litter      0.64318874
#Soil:Litter 0.09645036

# Post-Hoc comparisons using emmeans #
emm1 = emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")

# Assign letters to denote significant differences at p<0.05 #
cld1 <- cld(emm1, Letters=letters, sort=T)

# Subset columns of interest
cld1_sub <- cld1[,c(1,2,8)] 

# bind with summary data for plotting
bac_env_summary_15 <- subset(bac_env_summary, Timepoint=="15 Days")
bac_env_summary_15 <- inner_join(cld1_sub, bac_env_summary_15, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(bac_env_summary_15)

# Plot Results
Figure_S1B <- bac_env_summary_15 %>%
  mutate(Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  ggplot(bac_env_summary, mapping=aes(x=Litter, y=mean, fill=Mineral)) + 
  geom_bar(aes(alpha=Litter), stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
  geom_jitter(metadata_15, mapping = aes(x=Litter, y=copy_no)) +
  geom_text(aes(label = .group, y = mean + se), hjust = 0.65, vjust = -0.5) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"), guide="none") +
  facet_grid(Timepoint~Mineral) +
  ylim(0,4) +
  ylab("Abundance weighted mean \n rRNA gene copy number") +
  theme_bw() + theme(
    legend.position="none",
    text=element_text(color = "black", size=21),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    strip.text.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    strip.background = element_rect(fill = NA))

Figure_S1B

Figure_S1C <- metadata %>% 
  mutate(Timepoint = fct_relevel(Timepoint, "15 Days", "126 Days"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  filter(Timepoint=="15 Days", !Litter=="No Litter") %>%
  ggplot(mapping=aes(x=copy_no, y=respired_litter_c_perc, fill=Mineral, linetype=Mineral)) + 
  geom_point(aes(shape=Litter), size=2) + 
  facet_grid(~Mineral, scales="fixed") + 
  sm_statCorr(aes(group=Mineral), color="black", label_y = 0.8) +
  xlab("\n 16S rRNA Gene Copy Number") +
  ylab("% of litter-derived C respired") +
  xlim(2.5,4.3) +
  scale_linetype_manual(values = c("No Minerals" = "blank", "Kaolinite" = "blank", "Goethite" = "blank", "Montmorillonite" = "solid")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_shape_manual(values = c(24, 21)) +
  guides(fill = FALSE, linetype=FALSE) +
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

Figure_S1C

# Plot ordination of bacterial community composition at T15 alongside 16S copy number bar graphs
# and correlations with litter-C respired

Figure_S1 <- (Figure_S1A | (Figure_S1B/Figure_S1C)) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1.25, 1)) &
  theme(
    legend.direction = 'vertical',
    legend.box = 'vertical',
    plot.tag = element_text(size = 16, face = 'bold'))

Figure_S1

ggsave("Figure_S1.tiff", Figure_S1, height = 30, width = 55, units = "cm", compression="lzw")

# Calculate Bacterial Richness for Table 2 ####

# Richness and diversity Estimates on rarefied samples
b.rich <- estimate_richness(bac.r)
# Add metadata to dataframe 
b.rich <- data.frame(b.rich, sample_data(bac.r))

# Summarize the data by soil and litter treatments: Define function
my.summary = function(x) list(mean = mean(x), SE = parameters::standard_error(x))

b.rich <- data.table(b.rich)

# Select observed (number of ASV's) and Shannon diversity
b.rich_summary <- b.rich[, sapply(.SD, my.summary), by = c("Mineral", "Litter", "Timepoint"), .SDcols = c('Observed',"Shannon")]
setnames(b.rich_summary, 4:7, c("Observed", "Observed_SE", "Shannon", "Shannon_SE"))

b.rich_summary

write.csv(b.rich_summary, "Bacterial_Richness_Summary.csv")

# Model effects of soil and litter on bacterial richness at T15 and T126 (Table S4)

b.rich_T1 <- subset(b.rich, Timepoint=="15 Days")
b.rich_T2 <- subset(b.rich, Timepoint=="126 Days")

mod1 <- lm(Observed ~ Mineral*Litter, data=b.rich_T1)

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
#Mineral          5427  3  0.5081 0.6786
#Litter          13307  2  1.8690 0.1658
#Mineral:Litter  13519  6  0.6329 0.7031
#Residuals      163755 46    

#Assess relative importance of predictors #
calc.relimp(mod1, rela=FALSE)

#Proportion of variance explained by model: 16.29%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.02675177
#Litter         0.06703601
#Mineral:Litter 0.06910631

# Post-Hoc comparisons using emmeans #
emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")

mod2 <- lm(Observed ~ Mineral*Litter, data=b.rich_T2)

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
#Sum Sq Df F value Pr(>F)
#Mineral          4984  3  0.5059 0.6800
#Litter           3867  2  0.5888 0.5590
#Mineral:Litter   5692  6  0.2889 0.9394
#Residuals      157626 48           

#Assess relative importance of predictors #
calc.relimp(mod2, rela=FALSE)

#Proportion of variance explained by model: 8.45%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.02894946
#Litter         0.02245933
#Mineral:Litter 0.03306286

# Post-Hoc comparisons using emmeans #
emmeans(mod2, specs = pairwise ~ Litter|Mineral, type="response")

# END ####