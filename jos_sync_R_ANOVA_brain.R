# Brain data analysis for jos_sync2
# using this tutorial: https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#:~:text=The%20repeated%2Dmeasures%20ANOVA%20is,are%20measured%20more%20than%20once.&text=three%2Dway%20repeated%20measures%20ANOVA,on%20a%20continuous%20outcome%20variable.

# Beta (central cluster; 16-32Hz)-abs value
# 1 way repeated measure ANOVA between three conditions over the full gait cycle 
# Then 2 follow-up corrected t-tests between natural and sync and natural and control, controlling for false discovery rate (Benjamini and Yekutieli, 2001; Wagner et al., 2016)- expect main effect

# Alpha/mu (right and left parietal clusters; 7.5-12.5Hz)
# Data is not normally distributed sot he wilcoxon signed rank test is used between natural and sync and natural and control
#  controlling for false discovery rate (Benjamini and Yekutieli, 2001; Wagner et al., 2016)- expect main effect


library(tidyverse)
library(ggpubr)
library(rstatix)
library(MASS)

########################################################
# Beta

## Data preparation
# insert data here
# csv files
my_beta <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/data_measures/betacentr.csv")
head(my_beta,3)

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_beta <- my_beta %>%
  gather(key = "wcond", value = "power", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_beta, 3)

# Compute some summary statistics 
# (time): mean and sd (standard deviation)
my_beta %>%
  group_by(wcond) %>%
  get_summary_stats(power, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_beta, x = "wcond", y = "power", add = "point")
bxp

# look for outliers
my_beta %>%
  group_by(wcond) %>%
  identify_outliers(power)

# normality test: p>0.05 denotes normal distribution
my_beta %>%
  group_by(wcond) %>%
  shapiro_test(power)
  
# qq plot
  ggqqplot(my_beta, "power", facet.by = "wcond")
  
  res.aov <- anova_test(data = my_beta, dv = power, wid = id, within = wcond, effect.size = "pes")
  get_anova_table(res.aov, correction = "GG")

# post-hoc tests
#paired t-tests
# pairwise comparisons
  pwc <- my_beta %>%
    pairwise_t_test(
      power ~ wcond, paired = TRUE,
      p.adjust.method = "fdr"
    )
  pwc
  
# Visualization: box plots with p-values
  pwc <- my_beta %>% add_xy_position(x = "wcond")
  bxp + 
    stat_pvalue_manual(pwc) +
  labs(
     subtitle = get_test_label(res.aov, detailed = TRUE),
     caption = get_pwc_label(pwc)
   )
  
 # effect size
  my_beta  %>% cohens_d(power ~ wcond, paired = TRUE)
  ########################################################
# Part 2: 2 way Anova
# time = wcond; score = power; treatment = hemi
# alpha-mu
  
# load and show one row of data
# Wide format
  set.seed(123)
  my_alpha <- read.csv("C:/Users/ebmi2273/jos_sync/data/ana3_processed/data_measures/alphalat2.csv")
  my_alpha %>% sample_n_by(hemi, size = 1)
  head(my_alpha,3)  
  
  my_alpha_wil<- read.csv("C:/Users/ebmi2273/jos_sync/data/ana3_processed/data_measures/alphalat1.csv")
  head(my_alpha_wil,3) 
  
# Gather the columns t1, t2 and t3 into long format.
# Convert id and time into factor variables
  my_alpha <- my_alpha %>%
    gather(key = "wcond", value = "power", a_natural, b_control, c_sync) %>%
    convert_as_factor(id, wcond)
# Inspect some random rows of the data by groups
  set.seed(123)
  my_alpha %>% sample_n_by(hemi, wcond, size = 1)
  
# Group the data by treatment and time, and then 
# compute some summary statistics of the score variable: 
# mean and sd (standard deviation).
  my_alpha %>%
    group_by(hemi, wcond) %>%
    get_summary_stats(power, type = "mean_sd")
  
# Create box plots of the score colored by treatment groups:
  bxp <- ggboxplot(
    my_alpha, x = "wcond", y = "power",
    color = "hemi", palette = "jco"
  )
  bxp
  
# identify outliers
  my_alpha %>%
    group_by(hemi, wcond) %>%
  identify_outliers(power)
  
# Compute Shapiro-Wilk test for each combinations of factor levels:
# p>0.05 denotes normal distribution
  my_alpha %>%
    group_by(hemi, wcond) %>%
    shapiro_test(power)
  
# Create QQ plot for each cell of design:
  ggqqplot(my_alpha, "power", ggtheme = theme_bw()) +
    facet_grid(wcond ~ hemi, labeller = "label_both")

# cond effects
  wilcox.test( my_alpha_wil$a_natural_right,my_alpha_wil$b_control_right, paired=TRUE)
  wilcox.test( my_alpha_wil$a_natural_right,my_alpha_wil$c_sync_right, paired=TRUE)

  wilcox.test( my_alpha_wil$d_natural_left,my_alpha_wil$e_control_left, paired=TRUE)
  wilcox.test( my_alpha_wil$d_natural_left,my_alpha_wil$f_sync_left, paired=TRUE)

  # effect sizes
  my_alpha_wil <- my_alpha_wil %>%
    gather(key = "wcond", value = "power",  a_natural_right, b_control_right, c_sync_right, d_natural_left, e_control_left, f_sync_left)
  
wilcox_effsize(my_alpha_wil, power ~ wcond,  'paired'= TRUE)
