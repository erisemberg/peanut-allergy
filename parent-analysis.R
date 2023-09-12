### This script analyzes data from inbred CC027 and C3H mice (parents of the 
### CC027xC3H backcross)

library(pzfx)
library(MASS)
library(lme4)
library(rstatix)
library(RColorBrewer)
source("code-dependencies/qtl_functions.R")

load_themes() # load ggplot themes 
ensure_directory('results')
ensure_directory('figs')

# load parent data 
parentdata <- readr::read_csv("source_data/parent_challenge_data.csv") # temp 
parentdata <- as_tibble(parentdata)
parentige <- readr::read_csv("source_data/parent_PNsIgE.csv")[,c(1:3)] # PNsIgE

# pivot data longer for plotting
dat <- parentdata %>% 
  pivot_longer(cols = c('0', '15', '30', '45', '60'), names_to = 'Time', values_to = 'Temp')

# dark2 color palette 
dark2cols <- brewer.pal(n = 5, name = 'Dark2')
# replacing #E58700 (orange) and blue with Dark2 colors 4 and 5 (cluster using 1-3)

pparent <- ggplot(dat, aes(x = Time, y = Temp)) +
  geom_line(aes(group = mouse_ID), color = c(rep(dark2cols[4], 50), rep(dark2cols[5], 45)), alpha=0.25, size=0.2) +
  stat_summary(aes(group = Strain, color = paste("mean", Strain)), fun.y = mean, geom="line", show.legend = TRUE) +
  stat_summary(aes(group = Strain, color = paste("mean", Strain)), fun.data = mean_se, geom = "errorbar", width=0.15) +
  scale_colour_manual('Strain', values = c(dark2cols[4:5]), labels = c('C3H/HeJ', 'CC027')) +
  labs(x = 'Minutes post-challenge', y = 'Body Temperature') +
  scale_x_discrete(expand = c(0.03,0.03)) +
  ylim(32.5, 40)

print(pparent) # print 

# save 
png('figs/parent-temps.png', width = 600, height = 400)
pparent + big_theme
dev.off()

# save R object (for loading into main script / part of Fig 1)
ensure_directory("derived_data/other_Robjects")
save(pparent, file = 'derived_data/other_Robjects/pparent.Rdata')

## compare temperature trajectories in CC0027 vs C3H 
parentdata$Strain <- factor(parentdata$Strain) # convert strain to factor 

# find mean/sd of temps at each time point (for both CC027 and C3H mice)
dat %>% 
  group_by(Time) %>%
  get_summary_stats(Temp, type = 'mean_sd')

# identify outliers
dat %>% group_by(Time) %>%
  identify_outliers(Temp) 

# no extreme outliers
# C3H mouse who has a significant drop is an outlier but not extreme.

## Calculate derived measures    
# AAC statistic 
parentdata$Temp_aac <- rep(NA, nrow(parentdata))
x <- c(0,15,30,45,60)
for (i in 1:nrow(parentdata)){
  line <- parentdata[i,c('0','15','30','45','60')]
  auc <- auc(x = x, y = line)
  aac <- line[1]*tail(x, n=1) - auc
  parentdata[i, 'Temp_aac'] <- as.numeric(aac)
}

# minimum temperature 
parentdata$Min_Temp <- rep(NA, nrow(parentdata))
parentdata$Min_Temp <- apply(parentdata[,c('0','15','30','45','60')], 1, FUN=min)

# time of minimum temperature 
parentdata$Min_Temp_Time <- rep(NA, nrow(parentdata))
for (i in 1:nrow(parentdata)){
  parentdata$Min_Temp_Time[i] <- as.numeric(names(which.min(parentdata[i,c('0', '15', '30', '45', '60')])))
}

# doesn't make sense to cluster parent data because we only see the tri-modal distribution in the backcross mice

# save raw (unnormalized) data with derived measures 
write_csv(parentdata, file = 'derived_data/parentdata.csv') 


## Means/ranges in parent strains 
phenos <- c('0', '15', '30', '45', '60', 'Temp_aac', 'Min_Temp', 'Min_Temp_Time',
            'PNsIgE')
n <- length(phenos)

datasum <- data.frame(phenotype = phenos,
                      CC027mean = rep(NA, n),
                      CC027rangeL = rep(NA, n),
                      CC027rangeU = rep(NA, n),
                      C3Hmean = rep(NA, n),
                      C3HrangeL = rep(NA, n),
                      C3HrangeU = rep(NA, n))

for (i in 1:(n-1)){
  pheno <- phenos[i]
  dat <- parentdata[[pheno]]
  # CC027 inbreds
  datasum[i,'CC027mean'] <- mean(dat[which(parentdata$Strain == 'CC027')])
  datasum[i,'CC027rangeL'] <- range(dat[which(parentdata$Strain == 'CC027')])[1]
  datasum[i,'CC027rangeU'] <- range(dat[which(parentdata$Strain == 'CC027')])[2]
  # C3H inbreds
  datasum[i,'C3Hmean'] <- mean(dat[which(parentdata$Strain == 'C3H/HeJ')])
  datasum[i,'C3HrangeL'] <- range(dat[which(parentdata$Strain == 'C3H/HeJ')])[1]
  datasum[i,'C3HrangeU'] <- range(dat[which(parentdata$Strain == 'C3H/HeJ')])[2]
}

# Means/ranges for IgE (comes from different dataframe)
postparentige <- parentige %>% filter(`Pre/post` == 'Post-sensitization')
dat <- postparentige$Value

# CC027 inbreds
datasum[n,'CC027mean'] <- mean(dat[which(postparentige$Strain == 'CC027')])
datasum[n,'CC027rangeL'] <- range(dat[which(postparentige$Strain == 'CC027')])[1]
datasum[n,'CC027rangeU'] <- range(dat[which(postparentige$Strain == 'CC027')])[2]
# C3H inbreds
datasum[n,'C3Hmean'] <- mean(dat[which(postparentige$Strain == 'C3H/HeJ')])
datasum[n,'C3HrangeL'] <- range(dat[which(postparentige$Strain == 'C3H/HeJ')])[1]
datasum[n,'C3HrangeU'] <- range(dat[which(postparentige$Strain == 'C3H/HeJ')])[2]

# format for saving 
datasum$phenotype[1:5] <- c('T0', 'T15', 'T30', 'T45', 'T60')
# save
write_csv(datasum, file = 'results/inbred-parent-data-summary.csv')



## Normalize data    
# Test for normality with Shapiro-Wilks test. 
# Null hypothesis = data are normally distributed, 
# so p < 0.05 means the data are not normally distributed.
apply(parentdata[,c(4:11)], 2, FUN = shapiro.test)

# Only baseline temp is normally distributed with high confidence. 

# Normalize other phenotypes (except for min temp time which is categorical)
rawparentdata <- parentdata # save raw phenotypes
parentdata[,c(4:11)] <- apply(parentdata[,c(4:11)], MARGIN = 2, FUN = trint)

# Perform one-way repeated-measures ANOVA on the transformed temperature data 
# (just temp @ baseline, 15 min, 30 min, 45 min and 60 min)
dat <- parentdata %>% 
  pivot_longer(cols = c('0', '15', '30', '45', '60'), names_to = 'Time', values_to = 'Temp')

res.aov <- aov(Temp ~ Strain + Time + Error(mouse_ID), data = dat)
summary(res.aov)


## Calculate heritability  
# Rename time columns
names(parentdata)[4:8] <- c('T0', 'T15', 'T30', 'T45', 'T60')

# Calculate h^2 for all normal phenotypes using FSS/TSS
heritability <- data.frame(pheno = names(parentdata)[4:11], 
                           h2 = rep(NA, 8))

for (i in 1:nrow(heritability)){
  pheno <- heritability$pheno[i]
  
  if ((pheno == 'T0') | (pheno == 'Min_Temp_Time')){
    form <- paste(pheno, ' ~ Strain')
  } else {
    form <- paste0(pheno, ' ~ T0 + Strain')
  }
  
  aov <- anova(lm(form, parentdata))
  
  if ((pheno == 'T0') | (pheno == 'Min_Temp_Time')){
    heritability$h2[i] <- aov$`Sum Sq`[1] / sum(aov$`Sum Sq`) # 1st row = strain
  } else {
    heritability$h2[i] <- aov$`Sum Sq`[2] / sum(aov$`Sum Sq`) # 1st row = T0, 2nd row = strain
  }
}

# calculate heritability for IgE:
aov <- anova(lm(Value ~ Strain, data = postparentige))
h2ige <- aov$`Sum Sq`[1] / sum(aov$`Sum Sq`) # 1st row = strain
heritability[9,] <- c('PNsIgE', h2ige)

# format for saving 
heritability$h2 <- as.numeric(heritability$h2)
heritability$h2 <- sprintf(heritability$h2, fmt = '%#.2f')

# save 
write_csv(heritability, file = 'results/h2-from-inbred-parents.csv')

