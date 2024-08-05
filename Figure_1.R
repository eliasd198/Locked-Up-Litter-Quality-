###########################################################################################################
# Figure 1 - Litter and mineralogy treatment effects on litter and soil organic matter derived C respired.#
###########################################################################################################

# Set working Directory #
setwd("N:\\LockedUp\\Zenodo Upload")

# Load Required Packages #
library(olsrr)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(patchwork)
library(multcomp)
library(multcompView)
library(relaimpo)
library(emmeans)
library(ggpattern)

# Define any functions #
my.summary = function(x) list(mean = mean(x), SE = parameters::standard_error(x))

# Read Soil Respiration Data #
resp <- read.csv("13CO2_Summary_Cumulatives_mixingmodel.csv", header = T)
str(resp)

# Set data structure
resp$Date <- as.POSIXct(as.character(resp$sampling_date), format = "%d/%m/%Y")
resp$Mineral <- as.factor(resp$Mineral)
resp$Litter <- as.factor(resp$Litter)

# Change factor levels for plotting
resp$Mineral <- fct_recode(resp$Mineral, "No Minerals" = "Unamended")
resp$Litter <- fct_recode(resp$Litter, "Low Quality" = "Winter Wheat", "High Quality" = "White Clover")

# Change file to Data Table
DT <- data.table(resp, key = c("Treatment_ID", "DaysSinceAddn", "Jar_ID"))
DT

# Plotting cumulative litter-derived C respired as a % of initial litter-C addition ####

# And Summarize by treatment for plotting
cumPerc13C.plot <- DT[, sapply(.SD, my.summary), by = c("Treatment_ID", "DaysSinceAddn", "Mineral", "Litter"), .SDcols = c('cum_perc13c', 'cum_total_respired_c', 'cum_soil_derived_c', 'cum_primed_c')]
setnames(cumPerc13C.plot, 5:12, c("Cum.percent.13C", "SE", "total_soil_derived_c", "total_soil_SE", "soil_derived_c", "soil_SE", "primed_c", "primed_SE"))

# Plotting cumulative timeseries of litter-C respired ####
# Use a colour blind friendly palette for all manuscript figures #

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Define facet label names for Litter variable
litter.labs_soil <- c("A. No Litter","B. Wheat", "C. Clover")
names(litter.labs_soil) <- c("No Litter", "Low Quality", "High Quality")

Figure_1A <- cumPerc13C.plot %>%
  filter(!Litter =="No Litter") %>%
  mutate(Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Litter = fct_relevel(Litter, "No Litter", "High Quality", "Low Quality")) %>%
  ggplot(aes(x = DaysSinceAddn,y = Cum.percent.13C, color=Mineral, fill = Mineral))+
  geom_point(aes(shape=Litter), colour="black", size=3)+
  geom_line(aes(group=Litter), linewidth=0.8)+
  ylab(expression(paste("% of litter-derived C respired"))) +
  xlab("Days since litter addition") +
  scale_y_continuous(limits=c(0,100)) +
  facet_grid(~Mineral, labeller = labeller(Litter = litter.labs_soil)) +
  geom_errorbar(aes(ymin = Cum.percent.13C-SE, ymax = Cum.percent.13C+SE), color="black") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"), guide="none") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"), guide="none") +
  scale_shape_manual(values = c(21, 24)) +
  labs(title = NULL) +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     legend.key.size = unit(1, 'cm'),
                     legend.key=element_rect(fill="white"),
                     axis.title.y = element_text(size=14),
                     axis.title.x = element_blank(),
                     legend.position = "right",
                     axis.text = element_text(color = "black"),
                     panel.grid = element_blank(),
                     legend.title = element_blank(),
                     strip.text = element_text(face = "bold", hjust = 0),
                     strip.background.x = element_rect(color = NA, fill = NA))

Figure_1A

# Priming of native SOM ####

# read soils data #
c13 <- read.csv("litter_c_distribution.csv", header = T)
c13 <- subset(c13, !Litter=="No Litter") # Remove No Litter controls
str(c13)

# Set data structure
c13$Mineral <- as.factor(c13$Mineral)
c13$Litter <- as.factor(c13$Litter)

# Change factor levels for plotting
c13$Mineral <- fct_recode(c13$Mineral, "No Minerals" = "Unamended")
c13$Litter <- fct_recode(c13$Litter, "Low Quality" = "Winter Wheat", "High Quality" = "White Clover")


# Change file list to Data Table
DT_13c <- data.table(c13, key = c("Sample", "Mineral", "Litter"))
DT_13c

# Primed C #

# Model the effect of Mineral and Litter treatments on primed C in Soil and add letters to plot signifying differences between treatments  #
mod1 <- lm(primed_c_mug_g_T126 ~ Mineral*Litter, data=DT_13c)

# Pairwise Comparisons
emm1 = emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")

# Assign letters to denote significant differences at p<0.05 #
cld <- cld(emm1, Letters=letters, sort=T)

# Plotting Results #

# And Summarize by treatment for plotting
DT_13c_summary <- DT_13c[, sapply(.SD, my.summary), by = c("Mineral", "Litter"), .SDcols = c('primed_c_mug_g_T126')]
setnames(DT_13c_summary, 3:4, c("soil_13c", "SE"))

# add letters to denote significance to dt
cld_sub <- cld[,c(1,2,8)] # Subset columns of interest
DT_13c_summary <- inner_join(cld_sub, DT_13c_summary, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(DT_13c_summary)

# Plotting Boxplots
primed_list <- c("Winter Wheat" = "none", "White Clover" = "none")

# Define facet label names for Mineral variable
labs_primed <- c("A. No Minerals","B. Kaolinite", "C. Goethite", "D. Montmorillonite")
names(labs_primed) <- c("No Minerals", "Kaolinite", "Goethite", "Montmorillonite")

Figure_1B <- DT_13c %>%
  mutate(Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite"),
         Litter = fct_relevel(Litter, "High Quality", "Low Quality")) %>%
  ggplot(aes(x = Litter, y = primed_c_mug_g_T126, fill = Mineral)) +
  geom_boxplot(aes(fill=Mineral), outlier.shape = 21, outlier.fill = "white") +
  stat_summary(fun=mean, geom="point", shape=4, size=6, color="black") +
  ylab(expression(atop("Total primed C from SOC", paste("(", mu, g, " CO"[2],"-C",~g^-1," dry soil)")))) +
  geom_text(data = filter(DT_13c_summary, Litter =="Low Quality"), aes(label = .group, y = soil_13c + 50), hjust = 0.5, vjust = -0.5) +
  geom_text(data = filter(DT_13c_summary, Litter =="High Quality"), aes(label = .group, y = soil_13c + 65), hjust = 0.5, vjust = -0.5) +
  facet_grid(~Mineral) +
  scale_y_continuous(limits=c(-100,200)) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  labs(title = NULL) +
  scale_pattern_manual(values=primed_list) +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     panel.grid = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=14),
                     legend.title = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_text(angle = 45, hjust=1),
                     axis.text = element_text(color = "black"),
                     strip.text = element_blank(),
                     strip.background.x = element_rect(color = NA, fill = NA))

Figure_1B

# Create Figure 1

Figure_1 <- Figure_1A/Figure_1B +
  plot_layout(nrow=2,ncol=1, guides="collect") + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = 'bold'))

Figure_1

# Save Plot
setwd("N:\\LockedUp\\Figures\\Revised_Figures")

ggsave("Figure_1.tiff", Figure_1, height = 20, width = 30, units = "cm", compression="lzw")

# Paper Statistics ####

# Statistics for Paper

########################################
# Supplementary Information - Table S1 #
########################################

# Model the effect of Mineral and Litter treatments on litter-C respired #
mod_resp <- lm(sqrt(respired_litter_c_perc_T126) ~ Mineral*Litter, data=DT_13c)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_resp)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_resp)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_resp)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_resp)
# 5. Perform a Normality Test
ols_test_normality(mod_resp)

# Model summary and statistics #
summary(mod_resp)
car::Anova(mod_resp)

#Anova Table (Type II tests)

#Response: sqrt(respired_litter_c_perc_T126)
#Sum Sq Df F value    Pr(>F)    
#Mineral        12.5569  3 268.656 < 2.2e-16 ***
#  Litter          6.2230  1 399.426 < 2.2e-16 ***
#  Mineral:Litter  2.6831  3  57.406 5.633e-13 ***
#  Residuals       0.4986 32                   

# Plot Interaction
interaction.plot(
  x.factor = DT_13c$Litter,
  trace.factor = DT_13c$Mineral,
  response = DT_13c$respired_litter_c_perc_T126,
  fun = mean,
  ylab = "% Litter derived C respired",
  xlab = "Litter type",
  trace.label = "Mineral",
  col = c("#0198f9", "#f95801"),
  lyt = 1,
  lwd = 3
)
#Assess relative importance of predictors #
calc.relimp(mod_resp, rela=FALSE)

#Proportion of variance explained by model: 97.73%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.5717655
#Mineral:Litter 0.1221744
#Litter         0.2833589

# Post-Hoc comparisons using emmeans #
emmeans(mod_resp, specs = pairwise ~ Litter|Mineral, type="response")

# Model the effect of Mineral and Litter treatments on primed C in Soil  #
mod_prime <- lm(primed_c_mug_g_T126 ~ Mineral*Litter, data=DT_13c)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_prime)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_prime)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_prime)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_prime)
# 5. Perform a Normality Test
ols_test_normality(mod_prime)

# Model summary and statistics #
summary(mod_prime)
car::Anova(mod_prime)

#Anova Table (Type II tests)

#Response: primed_c_mug_g_T126
#Sum Sq Df F value    Pr(>F)    
#Mineral         22336  3  7.2701 0.0007492 ***
#  Litter          29977  1 29.2714  6.02e-06 ***
#  Mineral:Litter  15428  3  5.0218 0.0057719 ** 
#  Residuals       32771 32       

# Plot Interaction
interaction.plot(
  x.factor = DT_13c$Litter,
  trace.factor = DT_13c$Mineral,
  response = DT_13c$primed_c_mug_g_T126,
  fun = mean,
  ylab = "Primed C",
  xlab = "Litter type",
  trace.label = "Mineral",
  col = c("#0198f9", "#f95801"),
  lyt = 1,
  lwd = 3
)

#Assess relative importance of predictors #
calc.relimp(mod_prime, rela=FALSE, type = "lmg")

#Proportion of variance explained by model: 67.4%
#Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
#Mineral        0.2222212
#Mineral:Litter 0.1534977
#Litter         0.2982398

emmeans(mod_prime, specs = pairwise ~ Litter|Mineral, type="response")

#######
# END #
#######
