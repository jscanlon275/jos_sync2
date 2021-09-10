## jos_sync_R_ANOVA_beh.R
# cadence, stride time, frequency and phase locking
# analyze the data with anova, t-tests, wilcoxon signed rank test,
# find effect sizes

library(tidyverse)
library(ggpubr)
library(rstatix)
library(MASS)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Step speed - EXP
  # csv files
  my_estepspeed <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/EXP_steptime.csv")
head(my_estepspeed,3)

my_estepspeed_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/EXP_steptime.csv")
head(my_estepspeed_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_estepspeed <- my_estepspeed %>%
  gather(key = "wcond", value = "time", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_estepspeed, 3)

# Compute some summary statistics of the self-esteem score by groups 
# (time): mean and sd (standard deviation)
my_estepspeed %>%
  group_by(wcond) %>%
  get_summary_stats(time, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_estepspeed, x = "wcond", y = "time", add = "point")
bxp

# look for outliers
my_estepspeed %>%
  group_by(wcond) %>%
  identify_outliers(time)

# is the effect mostly due to outliers??
outliers <- boxplot(my_estepspeed$time, plot=FALSE)$out
nooutlier<-my_estepspeed 
nooutlier<- nooutlier[-which(nooutlier$time %in% outliers),]

#shapiro wilk normality test: p>0.05 denotes normal distribution
my_estepspeed %>%
  group_by(wcond) %>%
  shapiro_test(time)

# qq plot
ggqqplot(my_estepspeed, "time", facet.by = "wcond")


#cond effects
wilcox.test( my_estepspeed_wil$a_natural,my_estepspeed_wil$b_control, paired=TRUE)
wilcox.test( my_estepspeed_wil$a_natural,my_estepspeed_wil$c_sync, paired=TRUE)

# now do the effect sizes
wilcox_effsize(my_estepspeed, time ~ wcond, NULL, NULL, 'paired'= TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Step speed - PAR
  # csv files
  my_pstepspeed <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/PAR_steptime.csv")
head(my_pstepspeed,3)

my_pstepspeed_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/PAR_steptime.csv")
head(my_pstepspeed_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_pstepspeed <- my_pstepspeed %>%
  gather(key = "wcond", value = "time", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_pstepspeed, 3)

# Compute some summary statistics 
# (time): mean and sd (standard deviation)
my_pstepspeed %>%
  group_by(wcond) %>%
  get_summary_stats(time, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_pstepspeed, x = "wcond", y = "time", add = "point")
bxp

# look for outliers
my_pstepspeed %>%
  group_by(wcond) %>%
  identify_outliers(time)

#shapiro wilk normality test: p>0.05 denotes normal distribution
my_pstepspeed %>%
  group_by(wcond) %>%
  shapiro_test(time)

# qq plot
ggqqplot(my_pstepspeed, "time", facet.by = "wcond")

#cond effects
wilcox.test( my_pstepspeed_wil$a_natural,my_pstepspeed_wil$b_control, paired=TRUE)
wilcox.test( my_pstepspeed_wil$a_natural,my_pstepspeed_wil$c_sync, paired=TRUE)

# now do the effect sizes!
wilcox_effsize(my_pstepspeed, time ~ wcond, NULL, NULL, 'paired'= TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Step cadence - EXP
  # csv files
my_esteppermin <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/EXP_steppermin.csv")
head(my_esteppermin,3)

my_esteppermin_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/EXP_steppermin.csv")
head(my_esteppermin_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_esteppermin <- my_esteppermin %>%
  gather(key = "wcond", value = "time", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_esteppermin, 3)

# Compute some summary statistics 
# (time): mean and sd (standard deviation)
my_esteppermin %>%
  group_by(wcond) %>%
  get_summary_stats(time, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_esteppermin, x = "wcond", y = "time", add = "point")
bxp

# look for outliers
my_esteppermin %>%
  group_by(wcond) %>%
  identify_outliers(time)

#shapiro wilk normality test: p>0.05 denotes normal distribution
my_esteppermin %>%
  group_by(wcond) %>%
  shapiro_test(time)

# qq plot
ggqqplot(my_esteppermin, "time", facet.by = "wcond")

#cond effects
wilcox.test( my_esteppermin_wil$a_natural,my_esteppermin_wil$b_control, paired=TRUE)
wilcox.test( my_esteppermin_wil$a_natural,my_esteppermin_wil$c_sync, paired=TRUE)

# now do the effect sizes!
wilcox_effsize(my_esteppermin, time ~ wcond, NULL, NULL, 'paired'= TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Step cadence - PAR
# csv files
my_psteppermin <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/PAR_steppermin.csv")
head(my_psteppermin,3)

my_psteppermin_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/PAR_steppermin.csv")
head(my_psteppermin_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_psteppermin <- my_psteppermin %>%
  gather(key = "wcond", value = "time", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_psteppermin, 3)

# Compute some summary statistics 
# (time): mean and sd (standard deviation)
my_psteppermin %>%
  group_by(wcond) %>%
  get_summary_stats(time, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_psteppermin, x = "wcond", y = "time", add = "point")
bxp

# look for outliers
my_psteppermin %>%
  group_by(wcond) %>%
  identify_outliers(time)

#shapiro wilk normality test: p>0.05 denotes normal distribution
my_psteppermin %>%
  group_by(wcond) %>%
  shapiro_test(time)

# qq plot
ggqqplot(my_psteppermin, "time", facet.by = "wcond")

#cond effects
wilcox.test( my_psteppermin_wil$a_natural,my_psteppermin_wil$b_control, paired=TRUE)
wilcox.test( my_psteppermin_wil$a_natural,my_psteppermin_wil$c_sync, paired=TRUE)

# now do the effect sizes!
wilcox_effsize(my_psteppermin, time ~ wcond, NULL, NULL, 'paired'= TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## frequency locking percentage
# csv files
my_beh <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/perc_freqlocked.csv")
head(my_beh,3)

my_beh_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/perc_freqlocked.csv")
head(my_beh_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_beh <- my_beh %>%
  gather(key = "wcond", value = "time", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_beh, 3)

# Compute some summary statistics 
# (time): mean and sd (standard deviation)
my_beh %>%
  group_by(wcond) %>%
  get_summary_stats(time, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_beh, x = "wcond", y = "time", add = "point")
bxp

# look for outliers
my_beh %>%
  group_by(wcond) %>%
  identify_outliers(time)

#normality test: p>0.05 denotes normal distribution
my_beh %>%
  group_by(wcond) %>%
  shapiro_test(time)

# qq plot
ggqqplot(my_beh, "time", facet.by = "wcond")

#Now for the ANOVA
res.aov <- anova_test(data = my_beh, dv = time, wid = id, within = wcond, effect.size = "pes")
get_anova_table(res.aov, correction = "GG")

# post-hoc tests
#paired t-tests
# pairwise comparisons
pwc <- my_beh %>%
  pairwise_t_test(
    time ~ wcond, paired = TRUE,
    p.adjust.method = "fdr"
  )
pwc

# Visualization: box plots with p-values
pwc <- my_beh %>% add_xy_position(x = "wcond")
bxp + 
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

my_beh  %>% cohens_d(time ~ wcond, paired = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## in-phase locking percentage
# insert data here
# csv files
my_beh <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/perc_inphase_lock.csv")
head(my_beh,3)

my_beh_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/perc_inphase_lock.csv")
head(my_beh_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_beh <- my_beh %>%
  gather(key = "wcond", value = "stepperc", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_beh, 3)

# Compute some summary statistics of the self-esteem score by groups 
# (time): mean and sd (standard deviation)
my_beh %>%
  group_by(wcond) %>%
  get_summary_stats(stepperc, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_beh, x = "wcond", y = "stepperc", add = "point")
bxp

# look for outliers
my_beh %>%
  group_by(wcond) %>%
  identify_outliers(stepperc)

#normality test: p>0.05 denotes normal distribution
my_beh %>%
  group_by(wcond) %>%
  shapiro_test(stepperc)

# qq plot
ggqqplot(my_beh, "stepperc", facet.by = "wcond")

#cond effects
wilcox.test( my_beh_wil$a_natural,my_beh_wil$b_control, paired=TRUE)
wilcox.test( my_beh_wil$a_natural,my_beh_wil$c_sync, paired=TRUE)

# now do the effect sizes!
wilcox_effsize(my_beh, stepperc ~ wcond, NULL, NULL, 'paired'= TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  (in and anti-)phase locking percentage
# csv files
my_beh <- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/perc_inorantiphase_lock.csv")
head(my_beh,3)

my_beh_wil<- read.csv("C:/Users/ebmi2273/jos_sync2/data/ana3_processed/perc_inorantiphase_lock.csv")
head(my_beh_wil,3) 

# Gather columns t1, t2 and t3 into long format
# Convert id and time into factor variables
my_beh <- my_beh %>%
  gather(key = "wcond", value = "stepperc", a_natural, b_control, c_sync) %>%
  convert_as_factor(id, wcond)
head(my_beh, 3)


# Compute some summary statistics of the self-esteem score by groups 
# (time): mean and sd (standard deviation)
my_beh %>%
  group_by(wcond) %>%
  get_summary_stats(stepperc, type = "mean_sd")

# Create a box plot and add points corresponding to individual values:
bxp <- ggboxplot(my_beh, x = "wcond", y = "stepperc", add = "point")
bxp

# look for outliers
my_beh %>%
  group_by(wcond) %>%
  identify_outliers(stepperc)

#normality test: p>0.05 denotes normal distribution
my_beh %>%
  group_by(wcond) %>%
  shapiro_test(stepperc)

# qq plot
ggqqplot(my_beh, "stepperc", facet.by = "wcond")

#cond effects
wilcox.test( my_beh_wil$a_natural,my_beh_wil$b_control, paired=TRUE)
wilcox.test( my_beh_wil$a_natural,my_beh_wil$c_sync, paired=TRUE)

# now do the effect sizes!
wilcox_effsize(my_beh, stepperc ~ wcond, NULL, NULL, 'paired'= TRUE)
