# Linnea Sandell

# Functions for analysing growth curves
source("scripts/bioscreen_functions.R")
library(tidyverse)
library(lmerTest)


# Knockout of MAT and STE4 lead to increases in fitness -------------------

growth <- read.csv("data/bioscreens/nov2019_full.csv", stringsAsFactors = F, colClasses = c(rep("factor",4),
                                                                                            "integer",
                                                                                            rep("factor",5),
                                                                                            "logical",
                                                                                            "integer",
                                                                                            rep("factor",7),
                                                                                            "logical", 
                                                                                            "logical",
                                                                                            "factor",
                                                                                            "logical",
                                                                                            "logical",
                                                                                            rep("numeric",6)))
growth$petite <- factor(growth$petite, levels = c("FALSE", "TRUE"), labels = c("not petite", "petite"))
growth$MA <- factor(growth$MA, levels = c("FALSE", "TRUE"), labels = c("control", "MA"))
growth$KO <- factor(growth$KO, levels = c("FALSE", "TRUE"), labels = c("before KO", "after KO"))
growth <- growth %>% filter(!Aberrant) # 198 cases excluded

# Remove two lines that we don't trust
growth <- growth %>% filter(!anc_ID %in% c("81", "158"))

beforeKO <- growth %>% filter(KO == "before KO")
afterKO <- growth %>% filter(KO == "after KO")

bf_full <- lmer(maxSlope ~ MA + RDH + initOD + (1|Machine) + (1|Day), beforeKO, REML = F)
beforeKO$residuals <- resid(bf_full) 

beforeKO_sum <- beforeKO %>% group_by(KO, MA, Line2, anc_ID) %>% 
  summarise(mean_resid = mean(residuals), sem_resid =sd(residuals)/sqrt(n()))
var.test(beforeKO_sum$mean_resid[beforeKO_sum$MA == "MA"], beforeKO_sum$mean_resid[beforeKO_sum$MA == "control"])

ak_full <- lmer(maxSlope ~ MA + RDH + initOD + (1|Machine) + (1|Day), afterKO, REML = F)
afterKO$residuals <- resid(ak_full)
afterKO_sum <- afterKO %>% group_by(KO, MA, petite, Line2, anc_ID) %>% 
  summarise(mean_resid = mean(residuals), sem_resid =sd(residuals)/sqrt(n()))
var.test(afterKO_sum$mean_resid[afterKO_sum$MA == "MA"], afterKO_sum$mean_resid[afterKO_sum$MA == "control"])
var.test(afterKO_sum$mean_resid[afterKO_sum$MA == "MA" & afterKO_sum$petite == "not petite"], afterKO_sum$mean_resid[afterKO_sum$MA == "control" & afterKO_sum$petite == "not petite"])

before_after <- left_join(beforeKO_sum[, c("anc_ID", "MA", "mean_resid", "sem_resid")], afterKO_sum[,c("anc_ID","MA", "mean_resid", "sem_resid")], by = c("anc_ID", "MA"))
long_before_after <- bind_rows(beforeKO_sum, afterKO_sum)
cor.test(before_after$mean_resid.x, before_after$mean_resid.y)


before_after$change <- before_after$mean_resid.x - before_after$mean_resid.y
before_after$zscore <- sapply(before_after$change, function(x){(x-mean(before_after$change))/sd(before_after$change)})
hist(before_after$zscore, breaks = 25)


# Find outliers -----------------------------------------------------------
tmp <- before_after %>% filter(abs(zscore) > 2.5) # Much lower growth rates after knockout.
growth_sub <- growth %>% filter(!anc_ID %in% tmp$anc_ID)
full <- lmer(maxSlope ~ MA + KO + RDH + initOD + (1|Machine) + (1|Day) + (1|Line), growth_sub, REML = F)
summary(full)
full_noKO <- lmer(maxSlope ~ MA + RDH + initOD + (1|Machine) + (1|Day) + (1|Line), growth_sub, REML = F)
anova(full, full_noKO)

# Figure 2 ----------------------------------------------------------------
mylines <- before_after
mylines <- mylines %>% rename(beforeKO = mean_resid.x, afterKO = mean_resid.y, sembefore = sem_resid.x, semafter = sem_resid.y)
mylines_sub <- mylines %>% filter(!anc_ID %in% tmp$anc_ID, MA == "MA")
afterKO %>% filter(anc_ID %in% tmp$anc_ID) %>% summarise(unique(Line2)) # K2, 106, 122, 127
shapiro.test(mylines_sub$beforeKO)
shapiro.test(mylines_sub$afterKO) # not ok
ggpubr::ggqqplot(mylines$beforeKO)
ggpubr::ggqqplot(mylines$afterKO)


cor.test(mylines_sub$beforeKO, mylines_sub$afterKO) # P = 0.17, cor = 0.14
cor.test(mylines_sub$beforeKO, mylines_sub$afterKO, method = "kendall") # P = 0.217, cor = 0.087 NO CORRELATION

mylines$lose <- mylines$anc_ID %in% tmp$anc_ID
mylines$MAlose <- factor(paste(mylines$MA, mylines$lose), 
                         levels = c("control FALSE", "MA FALSE", "MA TRUE"), labels = c("control", "MA", "exclude"))
(p1 <- ggplot(mylines, aes(x=beforeKO, y = afterKO))+
    #geom_smooth(method = lm, show.legend = FALSE, linetype = "dashed", alpha = 0.5) +
    geom_smooth(data = mylines_sub, method = lm) +
    geom_errorbar(aes(ymin = afterKO - semafter, ymax = afterKO + semafter), alpha = 0.8) +
    geom_errorbarh(aes(xmin = beforeKO - sembefore, xmax = beforeKO + sembefore), alpha = 0.8) +
    geom_point(aes(fill = MAlose, shape = MAlose)) +
    ggpubr::stat_cor(data = mylines_sub, label.y.npc = 0.35, method = "kendall", cor.coef.name = "tau", show.legend = FALSE)+
    cowplot::theme_cowplot() +
    scale_shape_manual(values = c(21,21, 24)) +
    theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = c(0.05, 0.95), legend.background = element_rect(fill=alpha("gray", 0.2)))+
    scale_fill_manual(values = c(alpha("#F8766D", 0.75), alpha("#00BFC4", 0.75), "white")) +
    labs(y = "resid. Max slope after KO", x = "resid. Max slope before KO"))


mrev <- growth %>% group_by(MA, KO, Line, Line2, anc_ID) %>% summarise(stdev = sd(maxSlope), mean_max = mean(maxSlope), sem = sd(maxSlope)/sqrt(n())) 
mrev$lose <- mrev$anc_ID %in% tmp$anc_ID
mrev$MAlose <- factor(paste(mrev$MA, mrev$lose), levels = c("control FALSE", "MA FALSE", "MA TRUE"), labels = c("control", "MA", "exclude"))

grand <- mrev %>% filter(!lose) %>% group_by(KO, MA) %>% summarise(grand_mean = mean(mean_max), stdev = sd(mean_max), samplesize = n(), grand_sem = sd(mean_max)/sqrt(n()))

pd <- position_jitterdodge(jitter.width = 1)
(p3 <- ggplot(mrev, aes(x = KO, y = mean_max, fill = MAlose, group = KO))+
    #ggbeeswarm::geom_quasirandom(dodge.width = 1, shape = 21) +
    geom_linerange(aes(ymin = mean_max - sem, ymax = mean_max + sem), position = pd, alpha = 0.7) +
    geom_point(alpha = 0.75, position = pd, aes(fill = MAlose, shape = MAlose), show.legend = FALSE) +
    #geom_pointrange(aes(ymin = mean_max - sem, ymax = mean_max + sem), shape = 21, size = 0.4, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.5)) +
    geom_errorbar(data = grand, aes(y = grand_mean, ymin = grand_mean - grand_sem, ymax = grand_mean + grand_sem, fill = MA), position = position_dodge(width = 1))+
    geom_text(data = grand, aes(x=KO, y = 0.32, label = samplesize, fill = MA)) +
    cowplot::theme_cowplot()+
    facet_grid(~MA, )+
    scale_shape_manual(values = c(21,21,24)) +
    scale_fill_manual(values = c(alpha("#F8766D", 0.75), alpha("#00BFC4", 0.75), "white")) +
    ylim(0.225,0.32)+
    #theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = c(0.05, 0.99), legend.background = element_rect(fill=alpha("gray", 0.2)))+
    theme(axis.text.x = element_text(angle = 30, vjust=0.75), 
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    labs(x = element_blank(), y = "Max slope")
)


cowplot::plot_grid(p3, p1, labels = "AUTO", align = "h", axis = "bt")

ggsave("figures/MA_KO_comparison_meanMax.pdf", width = 8, height = 3.5)
ggsave("figures/MA_KO_comparison_meanMax.png", width = 8, height = 3.8)


# Increase in fitness -----
full <- lmer(maxSlope ~ KO + MA + RDH + initOD + (1|Machine) + (1|Day) + (1|anc_ID), growth_sub, REML = F)
summary(full) # t = 15.87

noKO <- lmer(maxSlope ~ MA + RDH + initOD + (1|Machine) + (1|Day) + (1|anc_ID), growth_sub, REML = F)
anova(full,noKO) # chisq = 241.36

# Knockout of MAT and STE4 lead to new variation in growth rates ----------
# Variance among lines in MA and control before and after KO
bf_MA <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day), beforeKO %>% filter(MA == "MA", !anc_ID %in% tmp$anc_ID), REML = F)
bf_MA_line <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day) + (1|Line), beforeKO %>% filter(MA == "MA", !anc_ID %in% tmp$anc_ID), REML = F)
anova(bf_MA, bf_MA_line)

bf_cont <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day), beforeKO %>% filter(MA != "MA", !anc_ID %in% tmp$anc_ID), REML = F)
bf_cont_line <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day) + (1|Line), beforeKO %>% filter(MA != "MA", !anc_ID %in% tmp$anc_ID), REML = F)
anova(bf_cont, bf_cont_line)

af_MA <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day), afterKO %>% filter(MA == "MA", !anc_ID %in% tmp$anc_ID), REML = F)
af_MA_line <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day) + (1|Line), afterKO %>% filter(MA == "MA", !anc_ID %in% tmp$anc_ID), REML = F)
anova(af_MA, af_MA_line)

af_cont <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day), afterKO %>% filter(MA != "MA", !anc_ID %in% tmp$anc_ID), REML = F)
af_cont_line <- lmer(maxSlope ~ RDH + initOD + (1|Machine) + (1|Day) + (1|Line), afterKO %>% filter(MA != "MA", !anc_ID %in% tmp$anc_ID), REML = F)
anova(af_cont, af_cont_line)



# Figure 3 ----------------------------------------------------------------

# September 2019
sep2019 <- read.csv("data/bioscreens/sep2019.csv") %>% filter(initOD < 0.27, initOD > 0.185) # 4 instances where the initial OD was very low or very high
blanks <- read.csv("data/bioscreens/sep2019_blanks.csv") %>% group_by(day, machine) %>% summarise(mIni = mean(init_OD))
sep2019 <- left_join(sep2019, blanks, by = c("machine", "day")) %>% mutate(machine = as.factor(machine), day = as.factor(day), mateA = as.factor(mateA), mateB = as.factor(mateB), genotype = as.factor(genotype))
sep2019$mIni[is.na(sep2019$mIni)] <- blanks$mIni[blanks$day == 6]
sep2019$petite <- factor(sep2019$petite, levels = c("FALSE", "TRUE"), labels = c("not petite","petite"))
sep2019$MA <- factor(sep2019$MA, levels = c("FALSE", "TRUE"), labels = c("control", "MA"))

#sep2019 <- sep2019 %>% filter(!mateA %in% c("2", "106", "122", "127")) # These are the previously identified outliers in terms of effect of KO

interaction.plot(response = sep2019$maxSlope, x.factor = sep2019$genotype, 
                 trace.factor = sep2019$mateA, legend = FALSE, lty = 1, xlab = "Genotype",
                 ylab = "Max slope", type = "b", pch = 16, las = 1, cex = 1.5)
ds1_noline <- lmerTest::lmer(maxSlope ~ petite + MA*genotype + RDH + initOD + (1|machine) + (1|day), sep2019, REML = F)
sep2019$resid_noline <- resid(ds1_noline)

ds1_summary <- sep2019 %>% group_by(Line.ID, mateA, genotype, petite, MA) %>% 
  summarise(meanSlope = mean(maxSlope), meanResid = mean(resid_noline), 
            samplesize = n(), se = sd(maxSlope)/sqrt(samplesize), resid_se = sd(resid_noline)/sqrt(samplesize), petite = unique(petite))
ds1_largesum <- ds1_summary %>% group_by(genotype, MA, petite) %>% summarise(grandMean = mean(meanSlope), grandResidMean = mean(meanResid), linesize = n(), grandSE = sd(meanSlope)/sqrt(linesize), grandResidSE = sd(meanResid)/sqrt(linesize))

ds1_summary$ds <- "DS1"
ds1_largesum$ds <- "DS1"

ds1_largesum <- droplevels(ds1_largesum)
ds1_summary <- droplevels(ds1_summary)



# January 2020
jan2020 <- read.csv("data/bioscreens/jan2020.csv")
blanks <- read.csv("data/bioscreens/jan2020_blanks.csv") %>% group_by(day) %>% summarise(mIni = mean(init_OD))
jan2020 <- left_join(jan2020, blanks, by = c("day")) %>% mutate(day = as.factor(day), mateA = as.factor(mateA), mateB = as.factor(mateB), genotype = as.factor(genotype))
jan2020$petite <- factor(jan2020$petite, levels = c("FALSE", "TRUE"), labels = c("not petite","petite"))
jan2020$MA <- factor(jan2020$MA, levels = c("FALSE", "TRUE"), labels = c("control", "MA"))

ds2_noline <- lmer(maxSlope ~ petite + MA*genotype + RDH + initOD + (1|day), jan2020, REML = F)
jan2020$resid_noline <- resid(ds2_noline)

ds2_summary <- jan2020 %>% group_by(Line.ID, mateA, RDH, genotype, MA, petite) %>% summarise(meanSlope = mean(maxSlope), meanResid = mean(resid_noline), samplesize = n(), se = sd(maxSlope)/sqrt(samplesize), resid_se =sd(resid_noline)/sqrt(samplesize), petite = unique(petite))
ds2_largesum <- ds2_summary %>% group_by(genotype, MA, petite) %>% summarise(grandMean = mean(meanSlope), grandResidMean = mean(meanResid), linesize = n(), grandSE = sd(meanSlope)/sqrt(linesize), grandResidSE = sd(meanResid)/sqrt(linesize))

ds2_summary$ds <- "DS2"
ds2_largesum$ds <- "DS2"
ds2_largesum <- droplevels(ds2_largesum)
ds2_summary <- droplevels(ds2_summary)



summary <- rbind(ds1_summary,ds2_summary)
summary$genotype <- as.character(summary$genotype)
summary$genotype[summary$genotype == "homozygote"] <- "homo."
summary$genotype[summary$genotype == "heterozygote"] <- "hetero."
#summary$genotype[summary$genotype == "haploid"] <- "haploid"
summary$genotype <- as.factor(summary$genotype)

largesum <- rbind(ds1_largesum, ds2_largesum)
largesum$genotype <- as.character(largesum$genotype)
largesum$genotype[largesum$genotype == "homozygote"] <- "homo."
largesum$genotype[largesum$genotype == "heterozygote"] <- "hetero."
#largesum$genotype[largesum$genotype == "haploid"] <- "hap"
largesum$genotype <- as.factor(largesum$genotype)

pd <- position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5)
summary %>% 
  mutate(MApetite = factor(paste0(MA, petite), levels = c("controlnot petite", "MAnot petite","controlpetite", "MApetite"), labels = c("control", "MA",  "control petite", "MA petite"))) %>% 
  ggplot(aes(x = genotype, y = meanSlope)) +
  geom_linerange(aes(ymin=meanSlope-se, ymax=meanSlope+se, col = MApetite, fill = MApetite),
                 position = pd)+
  geom_point(aes(fill = MApetite, color = MApetite), position = pd, shape = 21)+
  facet_grid(ds~MA) +
  geom_point(data=largesum, aes(y = grandMean, group = petite), position = position_dodge(width = 0.75)) +
  geom_errorbar(data=largesum, aes(y = grandMean, ymin = grandMean-grandSE, ymax = grandMean + grandSE, group = petite), width = 0.75, position = position_dodge(width = 0.75))+
  ylab("Maximum slope") +
  scale_fill_manual(values = c(alpha("#F8766D", 0.75), alpha("#00BFC4", 0.75), alpha("white",0.5), alpha("white",0.5)))+
  scale_color_manual(values = c("black", "black", "#F8766D", "#00BFC4")) +
  cowplot::theme_minimal_hgrid() +
  theme(axis.title.x = element_blank(), legend.title = element_blank(), legend.position = "bottom")
ggsave("figures/summary_petites.pdf", width = 6, height = 4)

# RDH54-deletion protects against petiteness ------------------------------

contable <- sep2019 %>% 
  filter(ploidy == "diploid") %>% 
  select(MA, RDH, petite, Line.ID) %>% 
  distinct()

chisq.test(contable$RDH, contable$petite)


contable <- jan2020 %>% 
  filter(ploidy == "diploid") %>% 
  select(MA, RDH, petite, Line.ID) %>% 
  distinct()

chisq.test(contable$RDH, contable$petite)



# No correlation in growth rates across data sets -------------------------

petite_homozygotes_2019 <- ds1_summary %>% filter(genotype == "homozygote", petite == "petite")
petite_heterozygotes_2019 <- ds1_summary %>% filter(genotype == "heterozygote", petite == "petite") 
grande_haploids_2019 <- ds1_summary %>% filter(genotype == "haploid", petite == "not petite")
grande_homozygotes_2019 <- ds1_summary %>% filter(genotype == "homozygote", petite == "not petite")
grande_heterozygotes_2019 <- ds1_summary %>% filter(genotype == "heterozygote", petite == "not petite")


petite_homozygotes_2020 <- ds2_summary %>% filter(genotype == "homozygote", petite == "petite")
petite_heterozygotes_2020 <- ds2_summary %>% filter(genotype == "heterozygote", petite == "petite") 
grande_haploids_2020 <- ds2_summary %>% filter(genotype == "haploid", petite == "not petite")
grande_homozygotes_2020 <- ds2_summary %>% filter(genotype == "homozygote", petite == "not petite")
grande_heterozygotes_2020 <- ds2_summary %>% filter(genotype == "heterozygote", petite == "not petite")


petite_homozygotes <- petite_homozygotes_2019 %>% filter(mateA %in% intersect(petite_homozygotes_2019$mateA, petite_homozygotes_2020$mateA)) %>% left_join(petite_homozygotes_2020, by = c("mateA", "genotype", "petite", "MA"))
petite_heterozygotes <- petite_heterozygotes_2019 %>% filter(mateA %in% intersect(petite_heterozygotes_2019$mateA, petite_heterozygotes_2020$mateA)) %>% left_join(petite_heterozygotes_2020, by = c("mateA", "genotype", "petite", "MA"))
grande_haploids <- grande_haploids_2019 %>% filter(mateA %in% intersect(grande_haploids_2019$mateA, grande_haploids_2020$mateA)) %>% left_join(grande_haploids_2020, by = c("mateA", "genotype", "petite", "MA"))
grande_homozygotes <- grande_homozygotes_2019 %>% filter(mateA %in% intersect(grande_homozygotes_2019$mateA, grande_homozygotes_2020$mateA)) %>% left_join(grande_homozygotes_2020, by = c("mateA", "genotype", "petite", "MA"))
grande_heterozygotes <- grande_heterozygotes_2019 %>% filter(mateA %in% intersect(grande_heterozygotes_2019$mateA, grande_heterozygotes_2020$mateA)) %>% left_join(grande_heterozygotes_2020, by = c("mateA", "genotype", "petite", "MA"))


cor.test(grande_haploids$meanResid.x[grande_haploids$MA == "MA"], grande_haploids$meanResid.y[grande_haploids$MA == "MA"], method = "kendall") # tau =0.14, p = 0.19
cor.test(petite_heterozygotes$meanResid.x[petite_heterozygotes$MA == "MA"], petite_heterozygotes$meanResid.y[petite_heterozygotes$MA == "MA"], method = "kendall") # tau = -0.16, p = 0.15
cor.test(grande_heterozygotes$meanResid.x[grande_heterozygotes$MA == "MA"], grande_heterozygotes$meanResid.y[grande_heterozygotes$MA == "MA"], method = "kendall") # tau =-0.03, p = 0.89
cor.test(petite_homozygotes$meanResid.x[petite_homozygotes$MA == "MA"], petite_homozygotes$meanResid.y[petite_homozygotes$MA == "MA"], method = "kendall") # tau = 0.04, p = 0.69
cor.test(grande_homozygotes$meanResid.x[grande_homozygotes$MA == "MA"], grande_homozygotes$meanResid.y[grande_homozygotes$MA == "MA"], method = "kendall") # tau =-0.08, p = 0.69



# Inferences of mutational effects across haploid, heterozygote an --------
sep_full <- lmer(maxSlope ~ petite + MA*genotype + RDH + initOD + (1|machine) + (1|day) + (1|mateA), sep2019, REML = F)
summary(sep_full)

jan_full <- lmer(maxSlope ~ petite + MA*genotype + RDH + initOD + (1|day) + (1|mateA), jan2020, REML = F)
summary(jan_full)

# Without haplotype id
sep_noHaplotypeID <- lmer(maxSlope ~ petite + MA*genotype + RDH + initOD + (1|machine) + (1|day) + (1|Line.ID), sep2019, REML = F)
anova(sep_full,sep_noHaplotypeID) # 

jan_noHaplotypeID <- lmer(maxSlope ~ petite + MA*genotype + RDH + initOD + (1|day) + (1|Line.ID), jan2020, REML = F)
anova(jan_full, jan_noHaplotypeID)


# Figure 4 ----------------------------------------------------------------


gdips_ds1 <- bind_rows(grande_homozygotes_2019, grande_heterozygotes_2019) %>% 
  ungroup %>% 
  select(mateA, genotype, petite, MA, meanSlope, se) %>% 
  pivot_wider(names_from = genotype, values_from = c(meanSlope, se)) %>% 
  mutate(ds = "DS1")

gdips_ds2 <- bind_rows(grande_homozygotes_2020, grande_heterozygotes_2020) %>% 
  ungroup %>% 
  select(mateA, genotype, petite, MA, meanSlope, se) %>% 
  pivot_wider(names_from = genotype, values_from = c(meanSlope, se)) %>% 
  mutate(ds = "DS2")

pdips_ds1 <- bind_rows(petite_homozygotes_2019, petite_heterozygotes_2019) %>% 
  ungroup %>% 
  select(mateA, genotype, petite, MA, meanSlope, se) %>% 
  pivot_wider(names_from = genotype, values_from = c(meanSlope, se)) %>% 
  mutate(ds = "DS1")

pdips_ds2 <- bind_rows(petite_homozygotes_2020, petite_heterozygotes_2020) %>% 
  ungroup %>% 
  select(mateA, genotype, petite, MA, meanSlope, se) %>% 
  pivot_wider(names_from = genotype, values_from = c(meanSlope, se))%>% 
  mutate(ds = "DS2")

bigdf <- bind_rows(gdips_ds1,gdips_ds2,pdips_ds1,pdips_ds2) %>% 
  mutate(MApetite = factor(paste0(MA, petite), levels = c("controlnot petite", "MAnot petite","controlpetite", "MApetite"), labels = c("control", "MA",  "control petite", "MA petite"))) %>% 
  drop_na() %>% mutate(MAlose = mateA %in% c("2", "106", "122", "127"))

ggplot(bigdf, aes(x=meanSlope_homozygote, y = meanSlope_heterozygote, group = MApetite, fill = MApetite, color = MApetite))+
  geom_smooth(method = lm, show.legend = FALSE, fill = "gray") +
  ggpubr::stat_cor(label.y = c(0.14, 0.14, 0.11, 0.11), method = "kendall", cor.coef.name = "tau", show.legend = FALSE)+
  facet_grid(ds~MA, scales = "free") +
  geom_errorbar(aes(ymin = meanSlope_heterozygote - se_heterozygote, ymax = meanSlope_heterozygote + se_heterozygote)) +
  geom_errorbarh(aes(xmin = meanSlope_homozygote - se_homozygote, xmax = meanSlope_homozygote + se_homozygote)) +
  geom_point(shape = 21) +
  labs(y = "Max slope heterozygote", x = "Max slope homozygote")+
  scale_fill_manual(values = c(alpha("#F8766D", 0.75), alpha("#00BFC4", 0.75), alpha("white",0.5), alpha("white",0.5)))+
  scale_color_manual(values = c("black", "black", "#F8766D", "#00BFC4")) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.title = element_blank(), legend.position = "bottom")
ggsave("figures/hetero_homo_correlation.pdf", height = 5, width =6)



# Figure S1: Correlation in growth between haploid KO two datasets -------------------
# Haploid KO lines that were not plasmid transformed were run once in September 2019 (DS1)
# and once in November 2019 (in the MA to KO comparison).
# Is there a correlation between their growth rates?

test <- afterKO_sum %>% 
  separate(Line2, into = c("extra", "Line.ID"), sep = 1) %>% 
  left_join(ds1_summary, by = c("MA", "petite", "Line.ID")) %>% 
  rename(meanSep = meanResid, meanNov = mean_resid,
         seSep = resid_se, seNov = sem_resid)

test$lose <- test$anc_ID %in% tmp$anc_ID
test$MAlose <- factor(paste(test$MA, test$petite, test$lose), 
                      levels = c("control not petite FALSE", 
                                 "MA not petite FALSE", 
                                 "MA petite FALSE", 
                                 "MA petite TRUE"), labels = c("control", "MA", "petite", "exclude"))
    
testsub <- filter(test, petite == "not petite")
cor.test(testsub$meanSep, testsub$meanNov, method = "kendall") # tau = 0.306, P < 10^-3

ggplot(test, aes(x=meanSep, y = meanNov))+
    #geom_smooth(method = lm, show.legend = FALSE, linetype = "dashed", alpha = 0.5) +
    geom_smooth(data=testsub, method = lm) +
    geom_errorbar(aes(ymin = meanNov - seNov, ymax = meanNov + seNov), alpha = 0.8) +
    geom_errorbarh(aes(xmin = meanSep - seSep, xmax = meanSep + seSep), alpha = 0.8) +
  geom_point(aes(fill = MAlose, shape = MAlose)) +
  ggpubr::stat_cor(data = testsub, label.y.npc = 0.15, label.x.npc = 0.65, 
                     method = "kendall", cor.coef.name = "tau", show.legend = FALSE)+
    cowplot::theme_cowplot() +
    scale_shape_manual(values = c(21,21, 21, 24)) +
  scale_fill_manual(values = c(alpha("#F8766D", 0.75), alpha("#00BFC4", 0.75), "white", "white")) +
    theme(legend.title = element_blank(), legend.direction = "horizontal", 
          legend.position = c(0.05, 0.98), 
          legend.background = element_rect(fill=alpha("gray", 0.2)))+
    labs(y = "resid. Max slope November", x = "resid. Max slope September")
ggsave("figures/FigureS1.pdf", height = 5, width =6)



# Difference between plasmids ---------------------------------------------
library(readxl)
knockout <- readxl::read_xlsx("raw-data/knockout_records.xlsx", sheet = "KO records") %>% 
  mutate(Line.ID =as.character(KO_code))

haploid_plsasmids <- jan2020 %>% 
  filter(ploidy == "haploid") %>% 
  select(MA, RDH, petite, Line.ID) %>% 
  distinct() %>% 
  left_join(select(knockout, Line.ID, haploid_treatment), by = "Line.ID")

haploid_plasmids %>% group_by(haploid_treatment, petite) %>% tally()

chisq.test(haploid_plasmids$haploid_treatment, haploid_plasmids$petite)
chisq.test(haploid_plasmids$haploid_treatment, haploid_plasmids$RDH)

haploid_plasmids_2 <- left_join(haploid_plasmids, ds2_summary, by = c("MA", "RDH", "petite", "Line.ID"))

ggplot(haploid_plasmids_2, aes(x=haploid_treatment, y =))