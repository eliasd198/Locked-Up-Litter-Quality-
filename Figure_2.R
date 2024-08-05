############################################################
# Figure 2 - Fate of carbon across experimental treatments #
############################################################

# Set working Directory #
setwd()

# Load Required Packages #
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(patchwork)
library(ggpubr)
library(olsrr)
library(emmeans)
library(relaimpo)
library(car)
library(multcomp)
library(tidyr)

# Define any functions #
my.summary = function(x) list(mean = mean(x), SE = parameters::standard_error(x))

# read soils data #
c13 <- read.csv("litter_c_data.csv", header = T)
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

# Model the effect of Mineral and Litter treatments on % 13C in MAOM and #
# indicate significant differences on boxplots with lowercase letters #
mod1 <- lm((maom_litter_c_perc_T126) ~ Mineral*Litter, data=DT_13c)

emm1 = emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")

# Assign letters to denote significant differences at p<0.05 #
cld1 <- cld(emm1, Letters=letters, sort=T)

# Plotting Results as boxplots #

# And Summarize by treatment for plotting
DT_13c_summary_maom <- DT_13c[, sapply(.SD, my.summary), by = c("Mineral", "Litter"), .SDcols = c('maom_litter_c_perc_T126')]
setnames(DT_13c_summary_maom, 3:4, c("maom_13c", "SE"))

# add letters to denote significance
cld1_sub <- cld1[,c(1,2,8)] # Subset columns of interest
DT_13c_summary_maom <- inner_join(cld1_sub, DT_13c_summary_maom, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(DT_13c_summary_maom)

################
# Figure 2A ####
################

# Colour blind friendly palette #

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Figure_2A <- DT_13c %>%
  mutate(Litter = fct_relevel(Litter, "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  ggplot(aes(x = Litter, y = maom_litter_c_perc_T126, fill = Mineral)) +
  geom_boxplot(alpha=1, outlier.shape = 21, outlier.fill = "white") +
  stat_summary(fun=mean, geom="point", shape=4, size=6, color="black") +
  ylab(expression(paste("% of litter-derived C recovered"))) +
  geom_text(data = DT_13c_summary_maom, aes(label = .group, y = maom_13c + 5), hjust = 0.5, vjust = -0.5) +
  scale_y_continuous(limits=c(20,55)) +
  facet_grid(~Mineral) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  labs(title = "A. MAOM") +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     panel.grid = element_blank(),
                     plot.title = element_text(size=14, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=14),
                     legend.title = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_blank(),
                     axis.text = element_text(color = "black"),
                     strip.text = element_text(size=14, face = "bold", hjust = 0),
                     strip.background.x = element_rect(color = "black", fill = NA))

Figure_2A

# Model the effect of Mineral and Litter treatments on % 13C in POM  #
mod2 <- lm(log(totalpom_litter_c_perc_T126) ~ Mineral*Litter, data=DT_13c)
emm2 = emmeans(mod2, specs = pairwise ~ Litter|Mineral, type="response")

# Assign letters to denote significant differences at p<0.05 #
cld2 <- cld(emm2, Letters=letters, sort=T)

# Plotting Results as boxplots #

# And Summarize by treatment for plotting
DT_13c_summary_pom <- DT_13c[, sapply(.SD, my.summary), by = c("Mineral", "Litter"), .SDcols = c('totalpom_litter_c_perc_T126')]
setnames(DT_13c_summary_pom, 3:4, c("pom_13c", "SE"))

# add letters to denote significance to dt
cld2_sub <- cld2[,c(1,2,8)] # Subset columns of interest
DT_13c_summary_pom <- inner_join(cld2_sub, DT_13c_summary_pom, by=c('Mineral'='Mineral', 'Litter'='Litter'))
print(DT_13c_summary_pom)

################
# Figure 2B ####
################

Figure_2B <- DT_13c %>%
  mutate(Litter = fct_relevel(Litter, "High Quality", "Low Quality"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Kaolinite", "Goethite", "Montmorillonite")) %>%
  ggplot(aes(x = Litter, y = totalpom_litter_c_perc_T126, fill = Mineral)) +
  geom_boxplot(alpha=1, outlier.shape = 21, outlier.fill = "white") +
  stat_summary(fun=mean, geom="point", shape=4, size=6, color="black") +
  ylab(expression(paste("% of litter-derived C recovered"))) +
  geom_text(data = DT_13c_summary_pom, aes(label = .group, y = pom_13c + 3), hjust = 0.5, vjust = -0.5) +
  scale_y_continuous(limits=c(0,15)) +
  facet_grid(~Mineral) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  labs(title = "B. POM") +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     panel.grid = element_blank(),
                     plot.title = element_text(size=14, face = "bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=14),
                     legend.title = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_text(angle = 45, hjust=1),
                     axis.text = element_text(color = "black"),
                     strip.text = element_text(size=14, face = "bold", hjust = 0),
                     strip.background.x = element_rect(color = "black", fill = NA))

Figure_2B

##################
# Figure 2C-D ####
##################

# Reshape relevant data from wide to long format. Express respired C as minus for plotting to represent loss from the microcosms
litter_c_data <- c13 %>% pivot_longer(cols = matches(c("respired_litter_c_perc_T126", "maom_litter_c_perc_T126", 
                                                       "totalpom_litter_c_perc_T126", "mbc_litter_c_perc_T126")), 
                                      names_to = "pool", values_to = "perc_13c") %>%
  dplyr::select(!6:14) %>% mutate(pool = fct_recode(pool, Respired="respired_litter_c_perc_T126", MAOM="maom_litter_c_perc_T126", 
                           POM="totalpom_litter_c_perc_T126", `Microbial Biomass`="mbc_litter_c_perc_T126"),
                           perc_13c = case_when(pool == "Respired" ~ -perc_13c,
                                                       pool == "MAOM" ~ perc_13c,
                                                       pool == "POM" ~ perc_13c,
                                                       pool == "Microbial Biomass" ~ perc_13c))
         
  
print(litter_c_data)

# Change dataframe to data.table
litter_c_data <- data.table(litter_c_data)

# Plotting a stacked bar chart ####
# Summarize data for plotting #
data_summary <- litter_c_data[, sapply(.SD, my.summary), by = c("Mineral", "Litter", "pool", "Timepoint"), .SDcols = c("perc_13c")]
setnames(data_summary, 5:6, c("perc_13c", "perc_SE"))

# Stacked Bar Plotting ####

# Add error bars to stacked bar plot in correct positions
data_summary$y_pos = NA
data_summary$y_pos[data_summary$pool == "Respired"] =  data_summary$perc_13c[data_summary$pool == "Respired"]
data_summary$y_pos[data_summary$pool == "Microbial Biomass"] = data_summary$perc_13c[data_summary$pool == "Microbial Biomass"]
data_summary$y_pos[data_summary$pool == "POM"] = data_summary$perc_13c[data_summary$pool == "Microbial Biomass"] + 
  data_summary$perc_13c[data_summary$pool == "POM"]
data_summary$y_pos[data_summary$pool == "MAOM"] = data_summary$y_pos[data_summary$pool == "POM"] + 
  data_summary$perc_13c[data_summary$pool == "MAOM"]

# Overlay data points to stacked bar plot in correct positions
litter_c_data$y_pos = NA
litter_c_data$y_pos[litter_c_data$pool == "Respired"] =  litter_c_data$perc_13c[litter_c_data$pool == "Respired"]
litter_c_data$y_pos[litter_c_data$pool == "Microbial Biomass"] = litter_c_data$perc_13c[litter_c_data$pool == "Microbial Biomass"]
litter_c_data$y_pos[litter_c_data$pool == "POM"] = litter_c_data$perc_13c[litter_c_data$pool == "Microbial Biomass"] + 
  litter_c_data$perc_13c[litter_c_data$pool == "POM"]
litter_c_data$y_pos[litter_c_data$pool == "MAOM"] = litter_c_data$y_pos[litter_c_data$pool == "POM"] + 
  litter_c_data$perc_13c[litter_c_data$pool == "MAOM"]

# Define facet label names for Litter variable
litter.labs <- c("C. High Quality","D. Low Quality")
names(litter.labs) <- c("White Clover", "Winter Wheat")

# Plot
Figure_2C_D <- data_summary %>%
  mutate(pool = fct_relevel(pool, "Respired", "MAOM", "POM", "Microbial Biomass"),
         Mineral = fct_relevel(Mineral, "No Minerals", "Goethite", "Kaolinite", "Montmorillonite")) %>%
  ggplot(data_summary, mapping = aes(x = Mineral,y = perc_13c, fill = pool)) +
  geom_bar(stat="identity", color="black") +
  geom_jitter(litter_c_data,  mapping = aes(x = Mineral,y = y_pos)) +
  geom_errorbar(aes(ymax = y_pos + perc_SE, ymin = y_pos - perc_SE), position = "identity", width = 0.5) +
  ylab(expression(paste("% of litter-derived C")))+
  #scale_y_continuous(limits=c(0,105)) +
  facet_grid(~Litter, labeller = labeller(Litter = litter.labs)) +
  scale_fill_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  scale_color_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme_bw() + theme(text=element_text(color = "black", size=21),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=14),
                     axis.text = element_text(color = "black"),
                     axis.text.x = element_text(angle=45, hjust=1),
                     panel.grid = element_blank(),
                     legend.title = element_blank(),
                     strip.text = element_text(face = "bold", hjust = 0),
                     strip.background.x = element_rect(color = NA, fill = NA))

Figure_2C_D

#############################
# Create composite Figure 2 #
#############################

figure_2 <- ((Figure_2A / Figure_2B + plot_layout(guides = 'auto')) | Figure_2C_D) + plot_layout(guides = 'collect', widths = c(2,1))
figure_2

ggsave("Figure_2.tiff", figure_2, height = 20, width = 38, units = "cm", compression="lzw")

# Paper Statistics #

# Statistics for Paper

########################################
# Supplementary Information - Table S2 #
########################################

# % litter-C recovery in MAOM #
mod_maom <- lm(maom_litter_c_perc_T126 ~ Mineral*Litter, data=DT_13c)

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

# Response: (maom_litter_c_perc)
# Sum Sq Df F value    Pr(>F)    
#  Mineral        324.57  3  44.352 1.663e-11 ***
#  Litter         131.59  1  53.944 2.370e-08 ***
#  Mineral:Litter  97.63  3  13.342 8.118e-06 ***
#  Residuals       78.06 32    

# Plot Soil x Litter Interaction
interaction.plot(x.factor     = DT_13c$Litter,
                 trace.factor = DT_13c$Mineral,
                 response     = DT_13c$maom_litter_c_perc_T126,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

#Assess relative importance of predictors #
calc.relimp(mod_maom, rela=FALSE, type = "lmg")

# Proportion of variance explained by model: 87.65%
# Metrics are not normalized (rela=FALSE). 

# Relative importance metrics: 
  
#  lmg
# Mineral        0.5136812
# Mineral:Litter 0.1545212
# Litter         0.2082581

# Post-Hoc comparisons using emmeans #
emmeans(mod_maom, specs = pairwise ~ Litter|Mineral, type="response")


# % litter-C recovery in POM #
mod_pom <- lm(log(totalpom_litter_c_perc_T126) ~ Mineral*Litter, data=DT_13c)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_pom)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_pom)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_pom)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_pom)
# 5. Perform a Normality Test
ols_test_normality(mod_pom)

# Model summary and statistics #
summary(mod_pom)
Anova(mod_pom)

# Anova Table (Type II tests)

# Response: log(totalpom_litter_c_perc)
# Sum Sq Df F value    Pr(>F)    
#  Mineral        2.5028  3 11.9458 2.069e-05 ***
#  Litter         6.8638  1 98.2805 2.802e-11 ***
#  Mineral:Litter 0.2197  3  1.0485    0.3846    
#Residuals      2.2348 32     

# Plot Soil x Litter Interaction
interaction.plot(x.factor     = DT_13c$Litter,
                 trace.factor = DT_13c$Mineral,
                 response     = DT_13c$totalpom_litter_c_perc_T126,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

#Assess relative importance of predictors #
calc.relimp(mod_pom, rela=FALSE, type = "lmg")

# Proportion of variance explained by model: 81.09%
# Metrics are not normalized (rela=FALSE). 

#Relative importance metrics: 
  
#  lmg
# Mineral        0.21172564
# Mineral:Litter 0.01858322
# Litter         0.58063658

# Post-Hoc comparisons using emmeans #
emmeans(mod1, specs = pairwise ~ Litter|Mineral, type="response")

# % litter-C recovery in MBC #
mod_mbc <- lm(mbc_litter_c_perc_T126 ~ Mineral*Litter, data=DT_13c)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_mbc)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_mbc)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_mbc)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_mbc)
# 5. Perform a Normality Test
ols_test_normality(mod_mbc)

# Model summary and statistics #
summary(mod_mbc)
Anova(mod_mbc)

# Anova Table (Type II tests)

# Response: mbc_litter_c_perc
# Sum Sq Df F value Pr(>F)
# Mineral        0.08913  3  1.2706 0.3010
# Litter         0.01862  1  0.7964 0.3788
# Mineral:Litter 0.09545  3  1.3608 0.2724
# Residuals      0.74822 32 

# Plot Soil x Litter Interaction
interaction.plot(x.factor     = DT_13c$Litter,
                 trace.factor = DT_13c$Mineral,
                 response     = DT_13c$mbc_litter_c_perc_T126,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

#Assess relative importance of predictors #
calc.relimp(mod_mbc, rela=FALSE, type = "lmg")

# Proportion of variance explained by model: 21.36%
# Metrics are not normalized (rela=FALSE). 

# Relative importance metrics: 
  
#  lmg
# Mineral        0.09367868
# Mineral:Litter 0.10032902
# Litter         0.01957089

# Post-Hoc comparisons using emmeans #
emmeans(mod_mbc, specs = pairwise ~ Litter|Mineral, type="response")

# % litter-C recovery in MBC at early-stage decomposition - 15 days #

mod_mbc_T15 <- lm(sqrt(mbc_litter_c_perc_T15) ~ Mineral*Litter, data=DT_13c)

# 1. Check the Residuals vs Fitted Plot
ols_plot_resid_fit(mod_mbc_T15)
# 2. Check the Normal Q-Q Plot 
ols_plot_resid_qq(mod_mbc_T15)
# 3. Create a Histogram of the Residuals
ols_plot_resid_hist(mod_mbc_T15)
# 4. Create a Boxplot II
ols_plot_resid_box(mod_mbc_T15)
# 5. Perform a Normality Test
ols_test_normality(mod_mbc_T15)

# Model summary and statistics #
summary(mod_mbc_T15)
Anova(mod_mbc_T15)

#Anova Table (Type II tests)

#Response: sqrt(mbc_litter_c_perc_T15)
#Sum Sq Df F value   Pr(>F)   
#Mineral        0.40236  3  5.1233 0.005759 **
#  Litter         0.23469  1  8.9651 0.005580 **
#  Mineral:Litter 0.20940  3  2.6663 0.066323 . 
#Residuals      0.75917 29   

# Plot Soil x Litter Interaction

interaction.plot(x.factor     = DT_13c$Litter,
                 trace.factor = DT_13c$Mineral,
                 response     = DT_13c$mbc_litter_c_perc_T15,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

#Assess relative importance of predictors #
calc.relimp(mod_mbc_T15, rela=FALSE, type = "lmg")

# Proportion of variance explained by model: 54.19%
# Metrics are not normalized (rela=FALSE). 

# Relative importance metrics: 
  
#  lmg
# Mineral        0.2583691
# Mineral:Litter 0.1263508
# Litter         0.1571979

# Post-Hoc comparisons using emmeans #
emmeans(mod1, specs = pairwise ~ Mineral|Litter, type="response")

#######
# END #
#######
