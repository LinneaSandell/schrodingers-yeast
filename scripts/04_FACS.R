# FACS data of selfed MA lines
# July 31st 2018
# Linnea Sandell

library(flowCore)
#========= Functions ============
peak1<-function(x){ # To get ploidy
  dens<-density(x$BL1.H) # We look at the density curve of green fluorescence
  peak.x<-vector() # BLI-H value
  peak.y<-vector() # Height (density)
  for(i in c(2:(length(dens$y)-1))){ # Walk along the density curve 
    if(dens$y[i]>dens$y[i-1] & dens$y[i]>dens$y[i+1]){ # If the point is higher than the previous and the next (local peak)
      peak.x<-c(peak.x, dens$x[i]) # Save the x-coordinate of the peak
      peak.y<-c(peak.y, dens$y[i]) # Save the y-coordinate of the peak
    }
  }
  y.true <- peak.y > max(peak.y)/4
  return(min(peak.x[y.true])) # Take the first peak that has a density at least higher than half of the max
}






# DS1 lines ---------------------------------------------------------------

#Screening function:
get.data <- function(file, trim=T, ratio_limit = 0.1){
  log<-capture.output({
    dx<-read.FCS(paste("raw-data/FACS_august/", file,sep=""), emptyValue=F, alter.names=T)})
  d<-as.data.frame(dx@exprs)
  if(trim){
    # get rid of events with BL1.A or FSC.A less than or equal to 20000 or equal to the max.
    d<-d[d$BL1.H>20000 & d$BL1.H<max(d$BL1.H) & 
           d$BL1.A<max(d$BL1.A) & d$FSC.H>0 & d$FSC.H<max(d$FSC.H) & d$FSC.A<max(d$FSC.A),]
    d$flag <- sapply(1:nrow(d), function(r){if(d$FSC.A[r] > d$FSC.H[r]*(1+ratio_limit) | d$FSC.A[r] < d$FSC.H[r]*(1-ratio_limit) |
                                               d$BL1.A[r] > d$BL1.H[r]*(1+ratio_limit) | d$BL1.A[r] < d$BL1.H[r]*(1-ratio_limit)){T}else{F}})
  }
  return(d)
}

codes<-read.csv("raw-data/FACS_august/FACS_codes_July.csv", stringsAsFactors = F)
codes<- codes[1:687,1:6]
codes$sample <- codes$Sample
files <- list.files("raw-data/FACS_august/", pattern = ".fcs")
codes$file <- paste("set",codes$Set, "_plate", codes$plate, "_", codes$FACS_well, ".fcs", sep = "")
#codes$bad <- FALSE
##codes[codes$Sample %in% codes[codes$Set == "3_rerun",]$Sample,]$bad <- TRUE
codes <- codes[codes$file %in% files,] #
codes <- codes[4:687,] #Because the unstained controls don't have any fluorescence

codes$ploidy_exp <- "haploid"
codes$ploidy_exp[grep("x", codes$sample)] <- "diploid"
codes$ploidy_exp[codes$sample == "dip-2"] <- "diploid"


codes$mateA <- sapply(codes$sample, function(x){strsplit(x, split = "x")[[1]][1]})
codes$mateB <- sapply(codes$sample, function(x){strsplit(x, split = "x")[[1]][2]})
codes$genotype <- ifelse(is.na(codes$mateB), "haploid", ifelse(codes$mateA == codes$mateB, "homozygote", "heterozygote"))
codes$genotype <- factor(codes$genotype, levels = c("haploid", "heterozygote", "homozygote"))

codes <- codes[order(codes$mateA, codes$genotype),]


# Stats on raw data files
codes$all_events <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r], trim = F)
  return(nrow(d))
})

codes$gate_events <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  return(nrow(d))
})

codes$flagged_events <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  return(sum(d$flag))
})

codes$sd_size <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  #d <- d[!d$flag,]
  if(nrow(d) == 0){return("na")}else{
    return(round(sd(d$FSC.H),1))}
})

codes$sd_ploidy <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  #d <- d[!d$flag,]
  if(nrow(d) == 0){return("na")}else{
    return(round(sd(d$BL1.H),1))}
})

codes$ploidy<-sapply(c(1:nrow(codes)),function(r){
  d <- get.data(codes$file[r])
  #d <- d[!d$flag,]
  p1<-peak1(d)
  return(p1)
}
)

# Change the mateA for the control lines so that they fall on top
codes$mateA[codes$sample == "dip-2"] <- "0"
codes$mateA[codes$sample == "a-2"] <- "0"


# Order everything based on mateA
codes$mateA <- as.numeric(codes$mateA)

# There's a case of 83B x 46A which is supposed to be the 83 heterozygote
# I guess we need to relabel it so that the mateA is 83 and mateB is 46
codes$sample[codes$sample == "46x83"] <- "83x46"
codes$mateA[codes$sample == "83x46"] <- 83
codes$mateB[codes$sample == "83x46"] <- 46
codes <- codes[order(codes$mateA, codes$genotype),]

plotorder <- rep(unique(codes$mateA), each = 6)


# Good.
pdf("figures/facs_august2019.pdf", onefile = TRUE)
par(mfrow = c(10,3))
par(mar = c(0,0,0,0))
par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
#par(fig=c(0.1,0.1,0.3,0.1))
r <- 1
i <- 1
while(i <= nrow(codes)){
    d <- get.data(codes$file[i])
    d <- d[!d$flag,]
    p1<-peak1(d)
    plot(density(d$BL1.H), axes = F, xlab = NA, ylab = NA, xlim = c(0,5e+05), ylim = c(0, 1.5e-05), main = codes$sample[i])
    abline(v = p1)
    text(x = 4e+05, y = 5e-06, labels = paste(nrow(d)," events", sep = ""))
    rect(xleft = p1 - round(sd(d$BL1.H)), ybottom = 0, xright = p1 + round(sd(d$BL1.H)), ytop = 5000, col = adjustcolor(ifelse(codes$ploidy_exp[i] == "haploid", "red", "blue"), alpha.f = 0.5))
    r <- r+1
    i <- i+1
    #if(r == 129) browser()
}
dev.off()














# DS2 lines ----------------------------------------------------------

#Screening function:
get.data <- function(file, trim=T, ratio_limit = 0.1){
  log<-capture.output({
    dx<-read.FCS(paste("raw-data/FACS_november/December_dominance/", file,sep=""), emptyValue=F, alter.names=T)})
  d<-as.data.frame(dx@exprs)
  if(trim){
    # get rid of events with BL1.A or FSC.A less than or equal to 20000 or equal to the max.
    d<-d[d$BL1.H>20000 & d$BL1.H<max(d$BL1.H) & 
           d$BL1.A<max(d$BL1.A) & d$FSC.H>0 & d$FSC.H<max(d$FSC.H) & d$FSC.A<max(d$FSC.A),]
    d$flag <- sapply(1:nrow(d), function(r){if(d$FSC.A[r] > d$FSC.H[r]*(1+ratio_limit) | d$FSC.A[r] < d$FSC.H[r]*(1-ratio_limit) |
                                               d$BL1.A[r] > d$BL1.H[r]*(1+ratio_limit) | d$BL1.A[r] < d$BL1.H[r]*(1-ratio_limit)){T}else{F}})
  }
  return(d)
}

#====== Example ===========
# General idea exemplified with single sample
d<-get.data("1", 1,"A6") # diploid control
plot(FSC.H~FSC.A, d, pch=".")

dens <- density(d$BL1.H)
plot(dens)
abline(v=dens$x[which(dens$y == max(dens$y))])
plot(density(d$BL1.H), xlim = c(0,4e+05)) # density is a function in the BiocGenerics package
abline(v= peak1(d))
peak1(d)



# For DS2 lines, complete January 2020 ------------------------------------


codes <- read.csv("raw-data/FACS_november/FACS_codes_december.csv", stringsAsFactors = F)
files <- list.files("raw-data/FACS_november/December_dominance/", pattern = "Box")

codes$file <- paste("Box",codes$FACS_plate, "_Dominance_Group_", codes$Row, codes$Column, ".fcs", sep = "")
codes <- rbind(codes, c("n/a", "n/a", "n/a", "106x106", "106_106.fcs"))
codes <- rbind(codes, c("n/a", "n/a", "n/a", "111x111", "111_111.fcs"))
codes <- rbind(codes, c("n/a", "n/a", "n/a", "111x93", "111_93.fcs"))
codes <- rbind(codes, c("n/a", "n/a", "n/a", "134x134", "134_134.fcs"))
codes <- rbind(codes, c("n/a", "n/a", "n/a", "33x33", "33_33.fcs"))
codes <- rbind(codes, c("n/a", "n/a", "n/a", "35x35", "35_35.fcs"))
# There were two samples that I labelled as uncertain but based on the order I had put them right. One of them is repeated
codes$sample[codes$sample == "4(uncertain)"] <- "4"
codes <- codes[-which(codes$sample == "89(uncertain)"),]


# Attach the lagging lines, and remove them from the first set
codes2 <- read.csv("data/facs_lagginglines3.csv")
codes2$file <- paste("lagginglines3_Dominance_Group_", codes2$Well, ".fcs", sep = "")
files2 <- list.files("raw-data/FACS_november/December_dominance/", pattern = "lagginglines3")
codes2 <- codes2[codes2$file %in% files2,]
codes2$FACS_plate <- "n/a"
codes2$Row <- substr(codes2$Well,1,1)
codes2$Column <- substr(codes2$Well,2,2)
codes2$sample <- as.character(codes2$Sample)
codes2 <- codes2[,c("FACS_plate", "Row", "Column", "sample", "file")]
codes <- codes[!codes$sample %in% codes2$sample,]
codes <- rbind(codes, codes2)


codes3 <- data.frame(sample = c("33x33", "111x93", "35x35", "134x134", "Blank", "35","a-2", "111", "35x46", "Blank", "48x48", "Blank", "111x111", "dip-2"),
                     Well = c("A6", "A7", "B6", "B7", "C6", "C7", "D6", "D7", "E6",  "E7", "F6", "F7", "G6", "H6"))
codes3$file <- paste("Lagginglines4_Dominance_Group_", codes3$Well, ".fcs", sep = "")
files3 <- list.files("raw-data/FACS_november/December_dominance/", pattern = "Lagginglines4")
codes3 <- codes3[codes3$file %in% files3,] # Three blanks were not run
codes3$FACS_plate <- "n/a"
codes3$Row <- substr(codes3$Well,1,1)
codes3$Column <- substr(codes3$Well,2,2)
codes3 <- codes3[,c("FACS_plate", "Row", "Column", "sample", "file")]
codes <- codes[!codes$sample %in% codes3$sample,]
codes <- rbind(codes, codes3)


codes$ploidy_exp <- "haploid"
codes$ploidy_exp[grep("x", codes$sample)] <- "diploid"
codes$ploidy_exp[codes$sample == "dip-2"] <- "diploid"


codes$mateA <- sapply(codes$sample, function(x){strsplit(x, split = "x")[[1]][1]})
codes$mateB <- sapply(codes$sample, function(x){strsplit(x, split = "x")[[1]][2]})
codes$genotype <- ifelse(is.na(codes$mateB), "haploid", ifelse(codes$mateA == codes$mateB, "homozygote", "heterozygote"))
codes$genotype <- factor(codes$genotype, levels = c("haploid", "heterozygote", "homozygote"))

codes <- codes[order(codes$mateA, codes$genotype),]


# Stats on raw data files
codes$all_events <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r], trim = F)
  return(nrow(d))
})

codes$gate_events <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  return(nrow(d))
})

codes$flagged_events <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  return(sum(d$flag))
})

codes$sd_size <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  #d <- d[!d$flag,]
  if(nrow(d) == 0){return("na")}else{
    return(round(sd(d$FSC.H),1))}
})

codes$sd_ploidy <- sapply(c(1:nrow(codes)), function(r){
  d <- get.data(codes$file[r])
  #d <- d[!d$flag,]
  if(nrow(d) == 0){return("na")}else{
    return(round(sd(d$BL1.H),1))}
})

codes$ploidy<-sapply(c(1:nrow(codes)),function(r){
  d <- get.data(codes$file[r])
  #d <- d[!d$flag,]
  p1<-peak1(d)
  return(p1)
}
)

# Change the mateA for the control lines so that they fall on top
codes$mateA[codes$sample == "dip-2"] <- "0"
codes$mateA[codes$sample == "a-2"] <- "0"


# Order everything based on mateA
codes$mateA <- as.numeric(codes$mateA)

# There's a case of 83B x 46A which is supposed to be the 83 heterozygote
# I guess we need to relabel it so that the mateA is 83 and mateB is 46
codes$sample[codes$sample == "46x83"] <- "83x46"
codes$mateA[codes$sample == "83x46"] <- 83
codes$mateB[codes$sample == "83x46"] <- 46
codes <- codes[order(codes$mateA, codes$genotype),]

plotorder <- rep(unique(codes$mateA), each = 3)


# Good.
pdf("figures/facs_complete_january2020.pdf", onefile = TRUE)
par(mfrow = c(10,3))
par(mar = c(0,0,0,0))
par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
#par(fig=c(0.1,0.1,0.3,0.1))
r <- 1
i <- 1
while(i <= nrow(codes)){
  if(codes$mateA[i] != plotorder[r]){
    plot(0, type = 'n', axes = FALSE)
    r <- r + 1
    #if(r == 129) browser()
  }else{
    d <- get.data(codes$file[i])
    d <- d[!d$flag,]
    p1<-peak1(d)
    plot(density(d$BL1.H), axes = F, xlab = NA, ylab = NA, xlim = c(0,5e+05), ylim = c(0, 1.5e-05), main = codes$sample[i])
    abline(v = p1)
    text(x = 4e+05, y = 5e-06, labels = paste(nrow(d)," events", sep = ""))
    rect(xleft = p1 - round(sd(d$BL1.H)), ybottom = 0, xright = p1 + round(sd(d$BL1.H)), ytop = 5000, col = adjustcolor(ifelse(codes$ploidy_exp[i] == "haploid", "red", "blue"), alpha.f = 0.5))
    r <- r+1
    i <- i+1
    #if(r == 129) browser()
  }
}
dev.off()


