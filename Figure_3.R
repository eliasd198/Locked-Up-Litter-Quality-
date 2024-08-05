###########################################################################################
# Figure 3 - MAOM Formation efficiency and total mineral-associated organic matter stocks #
###########################################################################################

# Set working directory #
setwd()

#remotes::install_github("coolbutuseless/ggpattern")
library("ggpattern")
# Load Required Packages #
library(olsrr)
library(data.table)
library(dplyr)
library(forcats)
library(ggplot2)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(relaimpo)
library(patchwork)

# Define any functions #
my.summary = function(x) list(mean = mean(x), SE = parameters::standard_error(x))

# Read Data #
c13 <- read.csv("litter_c_data.csv", header = T)

# Set treatments to factors
c13$Mineral <- as.factor(c13$Mineral)
c13$Litter <- as.factor(c13$Litter)
c13$Mineral <- fct_recode(c13$Mineral, "No Minerals" = "Unamended")
c13$Litter <- fct_recode(c13$Litter, "Low Quality" = "Winter Wheat", "High Quality" = "White Clover")

c13_sub <- subset(c13, !Litter=="No Litter") # Remove No Litter controls

# Change file list to Data Table
DT_13c <- data.table(c13_sub, key = c("Sample", "Mineral", "Litter"))
DT_13c

# Model the effect of Mineral and Litter treatments on C transfer efficiency #
mod_fe <- lm(sqrt(maom_formation_efficiency) ~ Mineral*Litter, data=DT_13c)

# Post-Hoc comparisons using emmeans #
emm1 = emmeans(mod_fe, specs = pairwise ~ Litter|Mineral, type="response")

# Assign letters to denote significant differences at p<0.05 #
cld <- cld(emm1, Letters=letters, sort=T)

# Plotting Results #

# summarize results by Mineral and Litter Treatments #
DT_13c_summary <- DT_13c[, sapply(.SD, my.summary), by = c("Mineral", "Litter"), .SDcols = c("maom_formation_efficiency")]
setnames(DT_13c_summary, 3:4, c("maom_fe", "maom_fe_SE"))

# add letters to denote significance to dt
cld_sub <- cld[,c(1,2,8)] # Subset columns of interest
DT_13c_summary <- inner_join(cld_sub, DT_13c_summary, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(DT_13c_summary)

# Plot the bar graph with raw data overlaid

# Define facet label names for Litter variable
litter.labs <- c("A. High Quality","B. Low Quality")
names(litter.labs) <- c("White Clover", "Winter Wheat")

# Plot
Figure_3A <- DT_13c_summary %>%
  mutate(Litter = fct_relevel(Litter, "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  ggplot(DT_13c_summary, mapping=aes(x = Litter,y = maom_fe, fill = Mineral)) +
  geom_bar(stat="identity", color="black", position = position_dodge()) +
  geom_jitter(DT_13c,  mapping = aes(x = Litter,y = maom_formation_efficiency)) +
  geom_errorbar(aes(ymin = maom_fe-maom_fe_SE, ymax = maom_fe+maom_fe_SE), width = 0.3, position = position_dodge(0.9)) +
  geom_text(aes(label = .group, y = maom_fe + maom_fe_SE), hjust = 0.65, vjust = -0.5) +
  ylab(expression(paste("MAOM formation efficiency")))+
  scale_y_continuous(limits=c(0,0.6)) +
  facet_grid(~Mineral, labeller = labeller(Litter = litter.labs)) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=14),
                     axis.text = element_text(color = "black"),
                     axis.text.x = element_text(angle=45, hjust=1),
                     panel.grid = element_blank(),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_text(face = "bold", hjust = 0),
                     strip.background.x = element_rect(color = NA, fill = NA))

Figure_3A

# reshape data to plot MAOM-C content partitioned into litter and SOM derived C
maom_part <- melt(setDT(c13[,c(4,5,16:18)]), measure.vars = c("litter_maom_c_content_mg_g","som_maom_c_content_mg_g", "total_maom_c_content_mg_g"), 
             value.name = c("maom_mg_g"), value.factor = TRUE, variable.name = "source")

maom_part$source <- fct_recode(maom_part$source, "Litter" = "litter_maom_c_content_mg_g", "SOM" = "som_maom_c_content_mg_g", "Total" = "total_maom_c_content_mg_g")
str(maom_part)

# Summarize data for plotting #
maom_part_summary <- maom_part[, sapply(.SD, my.summary), by = c("Mineral", "Litter", "source"), .SDcols = c('maom_mg_g')]
setnames(maom_part_summary, 4:5, c("maom", "maom_SE"))

# Figure 3B ####

# First model total MAOM stocks to add subscript letters to plot to indicate sig diffs between litters
maom_mod <- maom_part %>% filter(source=="Total", !Litter=="No Litter")
mod1 <- lm((maom_mg_g) ~ Litter*Mineral, data=maom_mod)

# Post-Hoc comparisons using emmeans #
emm1 = emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")
emm1

# Assign letters to denote significant differences at p<0.05 #
cld <- cld(emm1, Letters=letters, sort=T, alpha=0.05)

# add letters to denote significance to dt
cld_sub <- cld[,c(1,2,8)] # Subset columns of interest
dt <- inner_join(cld_sub, maom_part_summary, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(dt)

dt <- dt %>%
  filter(source =="Total") %>%
  mutate(source = fct_recode(source, "Litter derived"="Total"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality"))

 # Add asterisk to denote significance at p<0.1. Just do this manually
fix(dt)

# Add error bars to stacked bar plot in correct positions
maom_part_summary$y_pos = NA
maom_part_summary$y_pos[maom_part_summary$source == "SOM"] = maom_part_summary$maom[maom_part_summary$source == "SOM"]
maom_part_summary$y_pos[maom_part_summary$source == "Litter"] = maom_part_summary$maom[maom_part_summary$source == "SOM"] + 
  maom_part_summary$maom[maom_part_summary$source == "Litter"]

# Overlay data points to stacked bar plot in correct positions
maom_part_filter <- maom_part %>% filter(!source=="Total")

maom_part_filter$y_pos = NA
maom_part_filter$y_pos[maom_part_filter$source == "SOM"] = maom_part_filter$maom[maom_part_filter$source == "SOM"]
maom_part_filter$y_pos[maom_part_filter$source == "Litter"] = maom_part_filter$maom[maom_part_filter$source == "SOM"] + 
  maom_part_filter$maom[maom_part_filter$source == "Litter"]

# Remove Litter as a source from No Litter samples
maom_part_filter <- maom_part_filter %>% filter(!source=="Litter" | !Litter=="No Litter")
maom_part_filter$source = fct_recode(maom_part_filter$source, "Litter derived"="Litter", "SOM derived"="SOM")

# Define facet label names for Litter variable
litter.labs <- c("No Litter","Wheat", "Clover")
names(litter.labs) <- c("No Litter", "Winter Wheat", "White Clover")

list <- c("Litter derived" = "stripe", "SOM derived" = "none")

# Plot
Figure_3B <- maom_part_summary %>%
  filter(!source =="Total") %>%
  mutate(Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         source = fct_relevel(source, "Litter", "SOM"),
         source = fct_recode(source, "Litter derived"="Litter", "SOM derived"="SOM"),
         Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality")) %>%
  ggplot(maom_part_summary, mapping = aes(x = Litter,y = maom, fill = Mineral, pattern = source)) +
  ggpattern::geom_bar_pattern(stat="identity", colour  = 'black', show.legend = TRUE, pattern_fill = "white") +
  geom_jitter(maom_part_filter,  mapping = aes(x = Litter,y = y_pos, shape = source)) +
  geom_errorbar(aes(ymax = y_pos + maom_SE, ymin = y_pos - maom_SE), position = "identity", width = 0.5) +
  facet_grid(~Mineral) +
  geom_text(data = dt, aes(label = .group, y = maom + 1), hjust = 0.5, vjust = -0.5) +
  ylab(expression(paste("MAOM-C content ( mg C",~g^-1,"dry soil)"))) +
  scale_y_continuous(limits=c(0,20)) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_pattern_manual(values = list)  +
  guides(fill="none") +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     legend.key.size = unit(1.5, 'cm'),
                     legend.key=element_rect(fill="white"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=14),
                     axis.text = element_text(color = "black"),
                     axis.text.x = element_text(angle=45, hjust=1),
                     panel.grid = element_blank(),
                     legend.title = element_blank(),
                     strip.text = element_blank(),
                     strip.background.x = element_rect(color = NA, fill = NA))

Figure_3B

# Create Figure 3 #
figure_3 <- Figure_3A / Figure_3B +
  plot_annotation(tag_levels = 'A')

figure_3

ggsave("Figure_3.tiff", figure_3, height = 25, width = 35, units = "cm", compression="lzw")

# Paper Statistics #

# Statistics for Paper

########################################
# Supplementary Information - Table S3 #
########################################

# MAOM formation efficiency #

# Model the effect of Mineral and Litter treatments on maom formation efficiency #
mod_fe <- lm((maom_formation_efficiency) ~ Mineral*Litter, data=DT_13c)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_fe)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_fe)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_fe)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_fe)
# 5. Perform a Normality Test
ols_test_normality(mod_fe)

# Model summary and statistics #
summary(mod_fe)
Anova(mod_fe)

#Anova Table (Type II tests)

#Response: (maom_formation_efficiency)
#Sum Sq Df F value    Pr(>F)    
#Mineral        0.129312  3 183.942 < 2.2e-16 ***
#  Litter         0.054917  1 234.353  2.80e-16 ***
#  Mineral:Litter 0.038891  3  55.321  9.24e-13 ***
#  Residuals      0.007499 32                           


#Assess relative importance of predictors #
calc.relimp(mod_fe, rela=FALSE, type="lmg")

#Proportion of variance explained by model: 96.75%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.5607181
#Mineral:Litter 0.1686376
#Litter         0.2381287

# Plot Soil x Litter Interaction
interaction.plot(x.factor     = DT_13c$Litter,
                 trace.factor = DT_13c$Mineral,
                 response     = DT_13c$maom_formation_efficiency,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

# Post-Hoc comparisons using emmeans #
emmeans(mod_fe, specs = pairwise ~ Litter|Mineral, type="response")

# Total MAOM Stocks

mod_maom <- lm((maom_mg_g) ~ Litter*Mineral, data=maom_mod)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_maom)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_maom)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_maom)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_maom)
# 5. Perform a Normality Test
ols_test_normality(mod_maom)

# Model summary and statistics #
summary(mod_maom)
Anova(mod_maom)

# Anova Table (Type II tests)

# Response: (maom_mg_g)
# Sum Sq Df F value    Pr(>F)    
# Litter         0.0754  1  0.5538   0.46221    
# Mineral        4.2599  3 10.4322 6.087e-05 ***
# Litter:Mineral 0.9647  3  2.3624   0.08971 .  
# Residuals      4.3557 32   

#Assess relative importance of predictors #
calc.relimp(mod_maom, rela=FALSE, type="lmg")

# Proportion of variance explained by model: 54.89%
# Metrics are not normalized (rela=FALSE). 

# Relative importance metrics: 
  
#  lmg
# Mineral        0.441184922
# Litter:Mineral 0.099906177
# Litter         0.007806506

# Plot Soil x Litter Interaction
interaction.plot(x.factor     = maom_mod$Litter,
                 trace.factor = maom_mod$Mineral,
                 response     = maom_mod$maom_mg_g,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

# Post-Hoc comparisons using emmeans #
emmeans(mod_maom, specs = pairwise ~ Litter|Mineral, type="response")

# END #
