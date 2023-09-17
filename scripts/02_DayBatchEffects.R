# Linnea Sandell
# April 2021
# This script does two things:
#   it reports the variance attributable to day and machine in four datasets
#       - sep2016
#       - nov2019
#       - sep2019
#       - jan2020
#   it reports the correlation of identical lines measured across different datasets 
#       - haploid sep2016 to no knockout nov2019 (exact same lines, 3 years freezer)
#       - haploid knockout lines sep2019 and nov2019 (exact same lines, two months freezer)
#   it reports the correlation (or lack thereof) of the true replicates in sep2019 and jan2020.

library(performance)
library(tidyverse)


sep2016 <- read.csv("data/bioscreens/sep2016.csv") %>% 
  select(-X) %>% filter(Ploidy == "hap") %>% 
  mutate_if(is.integer, as.factor)

plot(maxSlope ~ initOD, sep2016)
abline(v=0.2)
abline(v=0.3)
sep2016 <- sep2016 %>% filter(initOD > 0.2, initOD < 0.3)
# Remove five lines we decided not to transform
sep2016 <- sep2016 %>% filter(!Line %in% c("54", "201", "29", "110", "189"))  

# Variance attributable to day and machine: 3.5005e-05
m1_2016 <- lmer(maxSlope ~ Group + RDH54 + initOD + (1|machine) + (1|day) + (1|Line), sep2016, REML = F)
summary(m1_2016) # Line ID variance 1.034e-05, day variance is 2.648e-05, machine 8.575e-06
icc(m1_2016, by_group = TRUE) # 0.039


nov2019 <- read.csv("data/bioscreens/nov2019_full.csv") %>% 
  filter(!Aberrant) %>%  # 198 cases excluded
  mutate_if(is.integer, as.factor) %>% 
  rename(machine = Machine, day = Day)

plot(maxSlope ~ initOD, nov2019)
abline(v=0.15)
abline(v=0.31)
nov2019 <- nov2019 %>% filter(initOD > 0.15, initOD < 0.31) # 7 points with outlier values for initial OD removed

# Variance attributable to day and machine: 3.2403e-04 (much higher)
m1_nov2019 <- lmer(maxSlope ~ MA + KO + RDH + initOD + (1|machine) + (1|day) + (1|Line), nov2019, REML = F)
summary(m1_nov2019) # Line ID variance 1.938e-04, day variance is 8.953e-05, machine 2.345e-04
icc(m1_nov2019, by_group = TRUE) # 0.177

# To make a comparison to the sep2016 data, we need to remove the knockout lines
noKO2019 <- nov2019 %>% filter(!KO)
# Variance attributable to day and machine: 2.82e-04 (much higher)
m2_nov2019 <- lmer(maxSlope ~ MA + RDH + initOD + (1|machine) + (1|day) + (1|Line), noKO2019, REML = F)
summary(m2_nov2019) # Line ID variance 5.503e-05, day variance is 6.841e-05, machine 2.136e-04
icc(m2_nov2019, by_group = TRUE) # 0.054


# Let's look to see if day effect is separate from batch effect when comparing sep2016 to nov2019
colnames(sep2016)
colnames(noKO2019)

noKO2019 <- noKO2019 %>% 
  select(-Line) %>% 
  mutate(batch = "nov2019", day = paste0("nov2019_",day), Group = if_else(MA,"MA","control")) %>% 
  rename(Line = anc_ID,
         RDH54 = RDH) %>% 
  select(Group, Line, RDH54, batch, day, machine, initOD, OD_maxSlope, finalOD, maxSlope, time_maxSlope, finalTime)

sep2016 <- sep2016 %>% mutate(batch = "sep2016", day = paste0("sep2016_", day)) %>% 
  select(Group, Line, RDH54, batch, day, machine, initOD, OD_maxSlope, finalOD, maxSlope, time_maxSlope, finalTime)

# First look at day vs batch effect
haploid_beforeKO <- rbind(noKO2019,sep2016) %>% mutate_if(is.character, as.factor)
m1_hap_beforeKO <- lmer(maxSlope ~ Group + RDH54 + initOD + (1|machine) + (1|day) + (1|Line), haploid_beforeKO, REML = F)
summary(m1_hap_beforeKO) # Line ID variance 2.069e-05, day variance is 1.553e-04, machine 2.926e-05
icc(m1_hap_beforeKO, by_group = TRUE) # 0.030

m2_hap_beforeKO <- lmer(maxSlope ~ Group + RDH54 + initOD + batch + (1|machine) + (1|day) + (1|Line), haploid_beforeKO, REML = F)
summary(m2_hap_beforeKO) # Line ID variance 2.062e-05, day variance is 3.432e-05, machine 2.476e-05
icc(m2_hap_beforeKO, by_group = TRUE) # 0.037


plot(maxSlope ~ batch, haploid_beforeKO)

# Is the model with batch effect better than without batch effect?
anova(m1_hap_beforeKO, m2_hap_beforeKO) # m2 significantly better

# September 2019 --------------------------------------------------------
sep2019 <- read.csv("data/bioscreens/sep2019.csv") 
plot(maxSlope ~ initOD, sep2019)
abline(v=0.27)
abline(v=0.185)
sep2019 <- sep2019 %>% filter(initOD < 0.27, initOD > 0.185) # 4 instances where the initial OD was very low or very high
blanks <- read.csv("data/bioscreens/sep2019_blanks.csv") %>% group_by(day, machine) %>% summarise(mIni = mean(init_OD))
sep2019 <- left_join(sep2019, blanks, by = c("machine", "day")) %>% mutate(machine = as.factor(machine), day = as.factor(day), mateA = as.factor(mateA), mateB = as.factor(mateB), genotype = as.factor(genotype))
sep2019$mIni[is.na(sep2019$mIni)] <- blanks$mIni[blanks$day == 6]
sep2019$petite <- factor(sep2019$petite, levels = c("FALSE", "TRUE"), labels = c("not petite","petite"))
sep2019$MA <- factor(sep2019$MA, levels = c("FALSE", "TRUE"), labels = c("control", "MA"))
sep2019$Line <- paste0(sep2019$mateA,sep2019$genotype)
  
# Variance attributable to day and machine: 3.2403e-04 (much higher)
m1_sep2019 <- lmer(maxSlope ~ petite + MA + RDH + ploidy + initOD + (1|machine) + (1|day) + (1|Line), sep2019, REML = F)
summary(m1_sep2019) # Line ID variance 1.728e-04, day variance is 2.427e-06, machine 1.648e-04
icc(m1_sep2019, by_group = TRUE) # 0.286

# Let's check nov2019 KO haploids to these haploids
# Need to first rbind the dataframes

KO2019 <- nov2019 %>% filter(KO)
