source("scripts/bioscreen_functions.R")
library(tidyverse)

# September 2016 ----------------------------------------------------------
codes <- read.csv("raw-data/bioscreens/sep2016/codes.csv")
locations <- read.delim("raw-data/bioscreens/sep2016/locations.txt", sep = "\t")
locations$Well.ID <- 100*locations$plateLR+locations$well
locations <- locations %>% mutate(file = paste0("bioscreen", machine, "_day", day, ".csv"))
files <- list.files("raw-data/bioscreens/sep2016", pattern = ".csv")
locations <- locations[locations$file %in% files,] #

codes <- left_join(codes, locations, by = c("Line" = "line"))
codes <- codes[order(codes$day),]
df <- sapply(files, function(x){read.csv(paste("raw-data/bioscreens/sep2016/", x, sep = ""), stringsAsFactors = F)})

pdf("figures/bioscreen_sep2016.pdf", onefile = TRUE)
par(mfrow = c(10,4))
par(mar = c(0,0,0,0))
par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
for(r in 1:nrow(codes)){
  temp <- df[codes$file[r]][[1]][paste0("Well.", codes$Well.ID[r])][[1]]
  time <- seq(0.25, length(temp)/4, by = 0.25)
  test <- NA
  test[1] <- FALSE
  for(i in 2:length(temp)){test[i] <- ifelse(abs(temp[i] - temp[max(which(!(test)))]) >= 0.15, T, F)}
  codes$keep[r] <- T
  if(sum(test) >=1){
    temp <- temp[-which(test)]
    time <- time[-which(test)]
    }
  growth <- spline.slope(time, temp)
  t <- spline(time, temp)
  if(is.na(t)){
    codes$keep[r] <- F
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- NA
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$day[r], ", Sample ", codes$Line[r], "\nr = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
    
  }else{
    if(sum(test[1:t]) >=1 | sum(test[1:75]) >= 5){
      codes$keep[r] <- F
    }
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- temp[which(time == t)]
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$day[r], ", Sample ", codes$Line[r], "\nr = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
    abline(v=t)
    }
}
dev.off()

sum(!codes$keep) # 342 growth curves were discarded
codes <- codes[codes$keep,]
write.csv(codes, "data/bioscreens/sep2016.csv", row.names = F)


# September 2019 -------------------------------------------------------------
codes <- read.csv("raw-data/bioscreens/sep2019/codes.csv")

locations <- read.csv("raw-data/bioscreens/sep2019/locations.csv")[,1:7]
locations$Well.ID <- 100*locations$plate+locations$well.ID
locations <- locations %>% mutate(file = paste0("dominance_day", day, "_machine", machine, ".csv"))
files <- list.files("raw-data/bioscreens/sep2019", pattern = "dominance")
locations <- locations[locations$file %in% files,] #

codes <- left_join(codes, locations)
codes <- codes[order(codes$day),]
df <- sapply(files, function(x){read.csv(paste0("raw-data/bioscreens/sep2019/", x), stringsAsFactors = F)})
dim(df)

# First subtract the OD of the blanks from all days.
blanks <- codes[which(codes$sample == "blank"),]
blanks$mean_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(mean(temp)) # Take the mean value
})
blanks$min_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(min(temp)) # Take the min value
})
blanks$init_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(mean(temp[3:10])) # Take the mean starting value
})

blanks <- select(blanks, day, machine, plate, Well.ID, mean_OD, min_OD, init_OD)

write.csv(blanks, "data/bioscreens/sep2019_blanks.csv", row.names = F)

codes <- filter(codes, sample != "blank")

pdf("figures/bioscreen_sep2019.pdf", onefile = TRUE)
par(mfrow = c(10,4))
par(mar = c(0,0,0,0))
par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
for(r in 1:nrow(codes)){
  temp <- df[paste0("Well.", codes$Well.ID[r]), codes$file[r]][[1]]
  time <- seq(0.25, length(temp)/4, by = 0.25)
  test <- NA
  test[1] <- FALSE
  for(i in 2:length(temp)){test[i] <- ifelse(abs(temp[i] - temp[max(which(!(test)))]) >= 0.15, T, F)}
  codes$keep[r] <- T
  if(sum(test) >=1){
    temp <- temp[-which(test)]
    time <- time[-which(test)]
  }
  growth <- spline.slope(time, temp)
  t <- spline(time, temp)
  if(is.na(t)){
    codes$keep[r] <- F
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- NA
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$day[r], ", Sample ", codes$sample[r], "\n r = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
  }else{
    if(sum(test[1:t]) >=1 | sum(test[1:75]) >= 5){
      codes$keep[r] <- F
    }
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- temp[which(time == t)]
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$day[r], ", Sample ", codes$sample[r], "\n r = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
    abline(v=t)
  }
}
dev.off()

sum(!codes$keep) # 179 growth curves were discarded
codes <- codes[codes$keep,]
codes$ploidy <- "diploid"
codes$ploidy[codes$type == "haploid"] <- "haploid"
codes$plasmidTransformed <- T
codes$plasmidTransformed[codes$type == "haploid"] <- F

codes <- codes %>% rename(Line.ID = sample, genotype = type) %>% 
  mutate(dataset = "Sep2019", MA = n_bott != "0", KO = T) %>% 
  select(Line.ID, dataset, machine, day, plate, 
         initOD, finalOD, maxSlope, OD_maxSlope, time_maxSlope, finalTime, 
         mateA, mateB, genotype, ploidy, RDH, petite, MA, KO, plasmidTransformed)

write.csv(codes, "data/bioscreens/sep2019.csv", row.names = F)


# November 2019 -------------------------------------------------------------
codes <- read.csv("raw-data/bioscreens/nov2019/codes.csv") # MA lines not petite

locations <- read.csv("raw-data/bioscreens/nov2019/locations.csv")
locations$Well.ID <- 100*locations$Plate+locations$Well

locations <- locations %>% mutate(file = paste0("KOMA_machine", Machine, "_day", Day, ".csv"))
files <- list.files("raw-data/bioscreens/nov2019", pattern = "KOMA_machine")
locations <- locations[locations$file %in% files,] #
locations$Line2 <- locations$Line
locations$Line2[grep("[.]", locations$Line)] <- sapply(locations$Line[grep("[.]", locations$Line)], 
                                                       function(x){
  substr(x, 1, nchar(x)-2)
})


codes <- left_join(locations, codes, by = c("Line2" = "ID"))
codes <- codes[order(codes$Day),]

# I would definitely just disregard day 1-3
codes <- codes[!codes$Day %in% c(1, 2, 3),]

df <- sapply(files, function(x){read.csv(paste0("raw-data/bioscreens/nov2019/", x), stringsAsFactors = F)})
dim(df)

# First subtract the OD of the blanks from all days.
blanks <- codes[which(codes$Line == "blank"),]
blanks$mean_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(mean(temp)) # Take the mean value
})
blanks$min_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(min(temp)) # Take the min value
})
blanks$init_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(mean(temp[3:10])) # Take the min value
})
unique(blanks$Volume)
blanks <- select(blanks, Day, Machine, Plate, Well.ID, mean_OD, min_OD, init_OD)

write.csv(blanks, "data/bioscreens/nov2019_blanks.csv", row.names = F)


codes <- filter(codes, Line != "blank")

pdf("figures/bioscreen_nov2019.pdf", onefile = TRUE)
par(mfrow = c(10,4))
par(mar = c(0,0,0,0))
par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
for(r in 1:nrow(codes)){
  temp <- df[paste0("Well.", codes$Well.ID[r]), codes$file[r]][[1]]
  time <- seq(0.25, length(temp)/4, by = 0.25)
  test <- NA
  test[1] <- FALSE
  for(i in 2:length(temp)){test[i] <- ifelse(abs(temp[i] - temp[max(which(!(test)))]) >= 0.15, T, F)}
  codes$keep[r] <- T
  if(sum(test) >=1){
    temp <- temp[-which(test)]
    time <- time[-which(test)]
  }
  growth <- spline.slope(time, temp)
  t <- spline(time, temp)
  if(is.na(t)){
    codes$keep[r] <- F
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- NA
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$Day[r], ", Sample ", codes$Line[r], "\n r = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
  }else{
    if(sum(test[1:t]) >=1 | sum(test[1:70]) >= 5){
      codes$keep[r] <- F
    }
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- temp[which(time == t)]
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$Day[r], ", Sample ", codes$Line[r], "\n r = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
    abline(v=t)
  }
}
dev.off()

sum(!codes$keep) # 77 growth curves were discarded
codes <- codes[codes$keep,]
write.csv(codes, "data/bioscreens/nov2019_full.csv", row.names = F)

codes$ploidy <- "haploid"
codes <- codes %>% rename(Line.ID = Line) %>% 
  mutate(dataset = "Nov2019", genotype = "haploid", mateB = NA, plasmidTransformed = F) %>% 
  select(Line.ID, dataset, Machine, Day, Plate, 
         initOD, finalOD, maxSlope, OD_maxSlope, time_maxSlope, finalTime, 
         mateA, mateB, genotype, ploidy, RDH, petite, MA, KO, plasmidTransformed)

write.csv(codes, "data/bioscreens/nov2019.csv", row.names = F)



# January 2020 -------------------------------------------------------------
codes <- read.csv("raw-data/bioscreens/jan2020/codes.csv")

locations <- read.csv("raw-data/bioscreens/jan2020/locations.csv")
locations$Well.ID <- 100+locations$Well.ID

locations <- locations %>% mutate(file = paste0("day", day, ".csv"))
files <- list.files("raw-data/bioscreens/jan2020", pattern = "day")
locations <- locations[locations$file %in% files,] #

codes <- left_join(codes, locations)
codes <- codes[order(codes$day),]

# We need to exclude the replicates that got too much or less media in it
codes <- codes[!(codes$day == 5 & codes$Well.ID == 278),]
codes <- codes[!(codes$day == 5 & codes$Well.ID == 288),]
codes <- codes[!(codes$day == 9 & codes$Well.ID == 139),]
codes <- codes[!(codes$day == 9 & codes$Well.ID == 140),]
codes <- codes[!(codes$day == 9 & codes$Well.ID == 132),]
codes <- codes[!(codes$day == 9 & codes$Well.ID == 133),]
codes <- codes[!(codes$day == 12 & codes$Well.ID == 176),]
codes <- codes[!(codes$day == 12 & codes$Well.ID == 187),]


df <- sapply(files, function(x){read.csv(paste0("raw-data/bioscreens/jan2020/", x), stringsAsFactors = F)})
dim(df)

# First subtract the OD of the blanks from all days.
blanks <- codes[which(codes$line == "blank"),]
blanks$mean_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(mean(temp)) # Take the mean value
})
blanks$min_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(min(temp)) # Take the min value
})
blanks$init_OD <- sapply(1:nrow(blanks),function(r){
  temp <- df[paste0("Well.", blanks$Well.ID[r]), blanks$file[r]][[1]]
  return(mean(temp[3:10])) # Take the min value
})

blanks$Plate <- 2
blanks$Plate[blanks$Well.ID <= 200] <- 1

blanks <- select(blanks, day, machine, Plate, Well.ID, mean_OD, min_OD, init_OD)

write.csv(blanks, "data/bioscreens/jan2020_blanks.csv", row.names = F)


codes <- filter(codes, line != "blank")
pdf("figures/bioscreen_jan2020.pdf", onefile = TRUE)
par(mfrow = c(10,4))
par(mar = c(0,0,0,0))
par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
for(r in 1:nrow(codes)){
  temp <- df[paste0("Well.", codes$Well.ID[r]), codes$file[r]][[1]]
  time <- seq(0.25, length(temp)/4, by = 0.25)
  test <- NA
  test[1] <- FALSE
  for(i in 2:length(temp)){test[i] <- ifelse(abs(temp[i] - temp[max(which(!(test)))]) >= 0.15, T, F)}
  codes$keep[r] <- T
  if(sum(test) >=1){
    temp <- temp[-which(test)]
    time <- time[-which(test)]
  }
  growth <- spline.slope(time, temp)
  t <- spline(time, temp)
  if(is.na(t)){
    codes$keep[r] <- F
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- NA
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$day[r], ", Sample ", codes$line[r], "\n r = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
  }else{
    if(sum(test[1:t]) >=1 | sum(test[1:75]) >= 5){
      codes$keep[r] <- F
    }
    codes$initOD[r] <- mean(temp[3:8]) # Arbitrary choice, first between 45 min and 2.5 hours
    codes$OD_maxSlope[r] <- temp[which(time == t)]
    codes$finalOD[r] <- mean(temp[(length(temp)-5):length(temp)])
    codes$maxSlope[r] <- growth
    codes$time_maxSlope[r] <- t
    codes$finalTime[r] <- length(temp)/4
    plot(temp ~ time, xlim = c(0, 24), ylim = c(0, 1.5), xlab = "Time (h)", ylab = "OD")
    text(5, 1, paste0("Day ", codes$day[r], ", Sample ", codes$line[r], "\n r = ", round(growth, 3), "\ntime points\nlost: ", sum(test)), cex = 0.8)
    if(!codes$keep[r]){text(16, 0.5, "Discard", col = "red")}
    abline(v=t)
  }
}
dev.off()

sum(!codes$keep) # 4 growth curves were discarded
codes <- codes[codes$keep,]

codes$Plate <- 2
codes$Plate[codes$Well.ID <= 200] <- 1

codes <- codes %>% rename(Line.ID = line) %>% 
  mutate(dataset = "Jan2020", plasmidTransformed = T, MA = !control, KO = T) %>% 
  select(Line.ID, dataset, machine, day, Plate, 
         initOD, finalOD, maxSlope, OD_maxSlope, time_maxSlope, finalTime, 
         mateA, mateB, genotype, ploidy, RDH, petite, MA, KO, plasmidTransformed)



write.csv(codes, "data/bioscreens/jan2020.csv", row.names = F)
