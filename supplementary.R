# This script is as part of the manuscript Circadian entrainment to red-light Zeitgebers and action spectrum for entrainment in the jewel wasp Nasonia vitripennis
# In this script, we plotted additional figures as part of the supplementary materials for the manuscript 

#load libraries
library(ggplot2)
library(data.table)
library(cowplot)
library(readxl)
library(dplyr)
library(officer)
library(rvg)
library(tidyverse)
library(zeitgebr)
library(sleepr)
library(ggetho)
library(behavr)
library(damr)
library(reshape)
library(reshape2)
library(ggpubr)
library(lme4)
library(nlme)
library(lsmeans)
library(drc)

# periodogram
DATA_DIR <-getwd()
list.files(DATA_DIR, pattern = ".txt|*.csv")
setwd(DATA_DIR)

metadata <- fread("metadata.csv")

metadata = metadata %>% mutate(start_datetime= "2022-04-20 10:00:00", stop_datetime = "2022-05-10 10:00:00")

metadata0.5 = subset(metadata, metadata$`temperature increase` == 0.5)
metadata1 = subset(metadata, metadata$`temperature increase` == 1)
metadata2 = subset(metadata, metadata$`temperature increase` == 2)
metadata3 = subset(metadata, metadata$`temperature increase` == 3)
metadata4 = subset(metadata, metadata$`temperature increase` == 4)
metadata5 = subset(metadata, metadata$`temperature increase` == 5)

metadata0.5 <- metadata0.5 %>% mutate(start_datetime = "2021-10-29 10:52:00", stop_datetime = "2021-11-15 12:11:00")
metadata1 <- metadata1 %>% mutate(start_datetime = "2021-11-16 12:30:00", stop_datetime = "2021-12-05 12:30:00")
metadata2 <- metadata2 %>% mutate(start_datetime = "2022-01-01 13:48:00", stop_datetime = "2022-01-20 13:48:00")
metadata3 <- metadata3 %>% mutate(start_datetime = "2022-01-21 16:47:00", stop_datetime = "2022-02-09 16:47:00")
metadata4 <- metadata4 %>% mutate(start_datetime = "2022-02-10 18:00:00", stop_datetime = "2022-03-03 18:00:00")
metadata5 <- metadata5 %>% mutate(start_datetime = "2022-03-04 21:00:00", stop_datetime = "2022-03-24 20:58:00")


metadata0.5 <- link_dam_metadata(metadata0.5,result_dir = DATA_DIR)
metadata1 <- link_dam_metadata(metadata1,result_dir = DATA_DIR)
metadata2 <- link_dam_metadata(metadata2,result_dir = DATA_DIR)
metadata3 <- link_dam_metadata(metadata3,result_dir = DATA_DIR)
metadata4 <- link_dam_metadata(metadata4,result_dir = DATA_DIR)
metadata5 <- link_dam_metadata(metadata5,result_dir = DATA_DIR)

dt0.5 <- load_dam(metadata0.5)
dt1 <- load_dam(metadata1)
dt2 <- load_dam(metadata2)
dt3 <- load_dam(metadata3)
dt4 <- load_dam(metadata4)
dt5 <- load_dam(metadata5)

metadata <- link_dam_metadata(metadata,result_dir = DATA_DIR)
dt1 <- load_dam(metadata)


dt1[, moving := activity >0]


d455 <- dt1[xmv(wavelength)==455]
d470 <- dt1[xmv(wavelength)==470]
d505 <- dt1[xmv(wavelength)==505]
d528 <- dt1[xmv(wavelength)==528]
d566 <- dt1[xmv(wavelength)==566]
d590 <- dt1[xmv(wavelength)==590]
d617 <- dt1[xmv(wavelength)==617]
d625 <- dt1[xmv(wavelength)==625]
d656 <- dt1[xmv(wavelength)==656]


d455[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d455[, .(id, uid), meta=T]
d505[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d505[, .(id, uid), meta=T]
d470[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d470[, .(id, uid), meta=T]
d528[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d528[, .(id, uid), meta=T]
d566[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d566[, .(id, uid), meta=T]
d590[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d590[, .(id, uid), meta=T]
d617[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d617[, .(id, uid), meta=T]
d625[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d625[, .(id, uid), meta=T]
d656[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d656[, .(id, uid), meta=T]



per_xsq_d455T <- periodogram(activity, 
                            d455,
                            FUN = chi_sq_periodogram)

per_xsq_d455T <- find_peaks(per_xsq_d455T)

per_xsq_d470T <- periodogram(activity, 
                            d470,
                            FUN = chi_sq_periodogram)

per_xsq_d470T <- find_peaks(per_xsq_d470T)

per_xsq_d505T <- periodogram(activity, 
                            d505,
                            FUN = chi_sq_periodogram)

per_xsq_d505T <- find_peaks(per_xsq_d505T)

per_xsq_d528T <- periodogram(activity, 
                            d528,
                            FUN = chi_sq_periodogram)

per_xsq_d528T <- find_peaks(per_xsq_d528T)

per_xsq_d566T <- periodogram(activity, 
                            d566,
                            FUN = chi_sq_periodogram)

per_xsq_d566T <- find_peaks(per_xsq_d566T)

per_xsq_d590T <- periodogram(activity, 
                            d590,
                            FUN = chi_sq_periodogram)

per_xsq_d590T <- find_peaks(per_xsq_d590T)

per_xsq_d617T <- periodogram(activity, 
                             d617,
                             FUN = chi_sq_periodogram)

per_xsq_d617T <- find_peaks(per_xsq_d617T)

per_xsq_d625T <- periodogram(activity, 
                             d625,
                             FUN = chi_sq_periodogram)

per_xsq_d625T <- find_peaks(per_xsq_d625T)
per_xsq_d656T <- periodogram(activity, 
                             d656,
                             FUN = chi_sq_periodogram)

per_xsq_d656T <- find_peaks(per_xsq_d656T)

p455 <- ggperio(per_xsq_d455T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p470 <- ggperio(per_xsq_d470T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p505 <- ggperio(per_xsq_d505T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p528 <- ggperio(per_xsq_d528T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p566 <- ggperio(per_xsq_d566T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p590 <- ggperio(per_xsq_d590T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p617 <- ggperio(per_xsq_d617T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p625 <- ggperio(per_xsq_d625T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p656 <- ggperio(per_xsq_d656T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()


# center of gravity for different T cycle

Tcycle <- c(19,24,23,22,21,20)


ldf <- list() #create a list
listtxt <- dir(pattern = "*Texp.txt") #creates the list of all the Texp. txt files in the directory (experimental groups)

for (k in 1:length(listtxt)){
  ldf[[k]] <- fread(listtxt[k])  #read all exp txt file
}

for (k in 1:length(listtxt)) {
  d <- ldf[[k]] #assign the first one to d
  d <- d[,-c(1, 5:9)] #delete not useable columns
  
  d <- d %>% dplyr::rename(Day="V2",minute="V3", DAM="V4",Light="V10",
                           "1"="V11","2"="V12","3"="V13","4"="V14","5"="V15",
                           "6"="V16","7"="V17","8"="V18","9"="V19","10"="V20",
                           "11"="V21","12"="V22","13"="V23","14"="V24","15"="V25",
                           "16"="V26","17"="V27","18"="V28","19"="V29","20"="V30",
                           "21"="V31","22"="V32","23"="V33","24"="V34","25"="V35",
                           "26"="V36","27"="V37","28"="V38","29"="V39","30"="V40",
                           "31"="V41","32"="V42")
  
  
  d <- melt(d, id.vars = c('Day', 'minute','DAM','Light'),variable.name = 'variable') #variable is FlyID
  d <- as.data.table(d)
  d <- d %>% mutate(Tcycle=Tcycle[k], group="experimental",wls="625 nm")
  
  
  #experimental group
  x <- Tcycle[k]*60 #minutes in this T cycle
  y <- length(d$minute)/32 #how many lines for each animal
  bin <- rep(seq(from=1, to=ceiling(y/10), by=1), each=10)
  bin = bin[1:y]
  
  bin2 <- rep(bin, 32)
  
  d$bin <- bin2
  
  q <- subset(d, d$variable == 1)
  q$time <- rep(seq(1,x/10,1), each = 10, ceiling(y/x), length.out = nrow(q))
  
  d$time <- rep(q$time,32)
  
  daggregate <- d[, .(mean = mean(value, na.rm = T)), by = c('bin','time','variable')]
  
  
  daggregate$Tcycle <- Tcycle[k]
  daggregate$group <- "exprimental"
  daggregate$wls <- "625 nm"
  
  day <- rep(seq(from=1, to=ceiling(y/x),by=1), each=x/10)
  z = length(daggregate$bin)/32
  day <- day[1:z]
  day <- rep(day, 32)
  
  daggregate$day <- day
  
  #calculate phase angle of entrainment
  q1 <- subset(daggregate, daggregate$variable ==1)
  angle <- rep(seq(2*pi/(x/10), 2*pi, 2*pi/(x/10)), ceiling(y/x), length.out = nrow(q1))
  #i <- length(unique(daggregate$time))
  #angle <- angle[1:i]
  
  daggregate$angle <- rep(angle, 32)
  daggregate$sin <- sin(daggregate$angle)
  daggregate$cos <- cos(daggregate$angle)
  daggregate$sinactivity <- daggregate$mean*daggregate$sin
  daggregate$cosactivity <- daggregate$mean*daggregate$cos
  
  #just plotting and check data
  ggplot(subset(daggregate, daggregate$variable == '10'), aes(x=time, y=mean))+
    geom_line()+
    facet_grid(day~ .) 
  
  daggregate2 <- daggregate[, .(meansin=mean(sinactivity), meancos=mean(cosactivity), meanact=mean(mean)), by=c("variable","day")]
  daggregate2$Tcycle <- Tcycle[k]
  daggregate2$group <- "exprimental"
  daggregate2$wls <- "625 nm"
  
  
  daggregate2$quard <- ifelse (daggregate2$meansin >0 ,
                               ifelse (daggregate2$meancos >0, 1,2),
                               ifelse(daggregate2$meancos <0, 3,4))
  daggregate2$tan <- daggregate2$meansin/daggregate2$meancos
  daggregate2$arctan <- atan(daggregate2$tan)
  daggregate2$angle <-  (daggregate2$arctan/(2*pi))*360
  
  daggregate2$correctedangle <- ifelse (daggregate2$quard == 1, daggregate2$angle,
                                        ifelse (daggregate2$quard == 2, daggregate2$angle+180,
                                                ifelse(daggregate2$quard ==3, daggregate2$angle+180,
                                                       daggregate2$angle+360 )))
  daggregate2$centralofgravity <- (daggregate2$correctedangle/360)*daggregate2$Tcycle
  daggregate2$centralofgravity_min <- (daggregate2$correctedangle/360)*(daggregate2$Tcycle*60/10)
  
  assign(paste("daggregate", k, sep="_"), daggregate)
  assign(paste("daggregate2", k, sep="_"), daggregate2)
}

d656_2 <- rbind(daggregate2_1,daggregate2_2,daggregate2_3,daggregate2_4,
                daggregate2_5,daggregate2_6)
d656 <- rbind(daggregate_1,daggregate_2,daggregate_3,daggregate_4,
              daggregate_5,daggregate_6)


d656_2a <- d656_2
d656_2a$centralofgravity <- d656_2a$centralofgravity + d656_2a$Tcycle

d656_2a <- bind_rows(d656_2,d656_2a)

d24 <- subset(d656_2a, Tcycle == 24)
d23 <- subset(d656_2a, Tcycle == 23)
d22 <- subset(d656_2a, Tcycle == 22)
d21 <- subset(d656_2a, Tcycle == 21)
d20 <- subset(d656_2a, Tcycle == 20)
d19 <- subset(d656_2a, Tcycle == 19)


p24 <- ggplot(d24, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))

p23 <- ggplot(d23, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,27), xmax = c(23,46),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,23),xmax=c(4,27),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,27),xmax=c(23,46),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,46))

p22 <- ggplot(d22, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,26), xmax = c(22,44),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,22),xmax=c(4,26),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,26),xmax=c(22,44),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,44))

p21 <- ggplot(d21, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))


p20 <- ggplot(d20, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,24), xmax = c(20,40),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,20),xmax=c(4,24),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,24),xmax=c(20,40),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40))


p19 <- ggplot(d19, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,23), xmax = c(19,38),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,19),xmax=c(4,23),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,23),xmax=c(19,38),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,38))

####################################################################
# center of gravity for different filters and wavelength

wavelength = c(455,470,505,528,540, 566, 590, 617,625,656)


ldf <- list() #create a list
listtxt <- dir(pattern = "*.txt") #creates the list of all the Texp. txt files in the directory (experimental groups)

for (k in 1:length(listtxt)){
  ldf[[k]] <- fread(listtxt[k])  #read all exp txt file
}

for (k in 1:length(listtxt)) {
  d <- ldf[[k]] #assign the first one to d
  d <- d[,-c(1, 5:9)] #delete not useable columns
  
  d <- d %>% dplyr::rename(Day="V2",minute="V3", DAM="V4",Light="V10",
                           "1"="V11","2"="V12","3"="V13","4"="V14","5"="V15",
                           "6"="V16","7"="V17","8"="V18","9"="V19","10"="V20",
                           "11"="V21","12"="V22","13"="V23","14"="V24","15"="V25",
                           "16"="V26","17"="V27","18"="V28","19"="V29","20"="V30",
                           "21"="V31","22"="V32","23"="V33","24"="V34","25"="V35",
                           "26"="V36","27"="V37","28"="V38","29"="V39","30"="V40",
                           "31"="V41","32"="V42")
  
  
  d <- melt(d, id.vars = c('Day', 'minute','DAM','Light'),variable.name = 'variable') #variable is FlyID
  d <- as.data.table(d)
  
  
  #experimental group
  x <- 21*60 #minutes in this T cycle
  y <- length(d$minute)/32 #how many lines for each animal
  bin <- rep(seq(from=1, to=ceiling(y/10), by=1), each=10)
  bin = bin[1:y]
  
  bin2 <- rep(bin, 32)
  
  d$bin <- bin2
  
  q <- subset(d, d$variable == 1)
  q$time <- rep(seq(1,x/10,1), each = 10, ceiling(y/x), length.out = nrow(q))
  
  d$time <- rep(q$time,32)
  
  daggregate <- d[, .(mean = mean(value, na.rm = T)), by = c('bin','time','variable')]
  
  
  day <- rep(seq(from=1, to=ceiling(y/x),by=1), each=x/10)
  z = length(daggregate$bin)/32
  day <- day[1:z]
  day <- rep(day, 32)
  
  daggregate$day <- day
  daggregate$wavelength = wavelength[k]
  
  #calculate phase angle of entrainment
  q1 <- subset(daggregate, daggregate$variable ==1)
  angle <- rep(seq(2*pi/(x/10), 2*pi, 2*pi/(x/10)), ceiling(y/x), length.out = nrow(q1))
  #i <- length(unique(daggregate$time))
  #angle <- angle[1:i]
  
  daggregate$angle <- rep(angle, 32)
  daggregate$sin <- sin(daggregate$angle)
  daggregate$cos <- cos(daggregate$angle)
  daggregate$sinactivity <- daggregate$mean*daggregate$sin
  daggregate$cosactivity <- daggregate$mean*daggregate$cos
  
  #just plotting and check data
  ggplot(subset(daggregate, daggregate$variable == '10'), aes(x=time, y=mean))+
    geom_line()+
    facet_grid(day~ .) 
  
  daggregate2 <- daggregate[, .(meansin=mean(sinactivity), meancos=mean(cosactivity), meanact=mean(mean)), by=c("variable","day")]
  daggregate2$Tcycle <- 21
  daggregate2$group <- "exprimental"
  daggregate2$wavelength <- wavelength[k]
  
  
  daggregate2$quard <- ifelse (daggregate2$meansin >0 ,
                               ifelse (daggregate2$meancos >0, 1,2),
                               ifelse(daggregate2$meancos <0, 3,4))
  daggregate2$tan <- daggregate2$meansin/daggregate2$meancos
  daggregate2$arctan <- atan(daggregate2$tan)
  daggregate2$angle <-  (daggregate2$arctan/(2*pi))*360
  
  daggregate2$correctedangle <- ifelse (daggregate2$quard == 1, daggregate2$angle,
                                        ifelse (daggregate2$quard == 2, daggregate2$angle+180,
                                                ifelse(daggregate2$quard ==3, daggregate2$angle+180,
                                                       daggregate2$angle+360 )))
  daggregate2$centralofgravity <- (daggregate2$correctedangle/360)*daggregate2$Tcycle
  daggregate2$centralofgravity_min <- (daggregate2$correctedangle/360)*(daggregate2$Tcycle*60/10)
  
  assign(paste("daggregate", k, sep="_"), daggregate)
  assign(paste("daggregate2", k, sep="_"), daggregate2)
}

dd <- rbind(daggregate2_1,daggregate2_2,daggregate2_3,daggregate2_4,
                daggregate2_5,daggregate2_6, daggregate2_7, daggregate2_8, daggregate2_9, daggregate2_10)
ddd<- rbind(daggregate_1,daggregate_2,daggregate_3,daggregate_4,
              daggregate_5,daggregate_6, daggregate_7, daggregate_8, daggregate_9, daggregate_10)


dd_2 <- dd
dd_2$centralofgravity <- dd_2$centralofgravity + dd_2$Tcycle

dd_2a <- bind_rows(dd,dd_2)

d455 <- subset(dd_2a, wavelength == 455)
d470 <- subset(dd_2a, wavelength == 470)
d505 <- subset(dd_2a, wavelength == 505)
d528 <- subset(dd_2a, wavelength == 528)
d566 <- subset(dd_2a, wavelength == 566)
d590 <- subset(dd_2a, wavelength == 590)
d617 <- subset(dd_2a, wavelength == 617)
d625 <- subset(dd_2a, wavelength == 625)
d656 <- subset(dd_2a, wavelength == 656)


p455 <- ggplot(d455, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p470 <- ggplot(d470, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p505 <- ggplot(d505, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p528 <- ggplot(d528, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p566 <- ggplot(d566, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p590 <- ggplot(d590, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p617 <- ggplot(d617, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p625 <- ggplot(d625, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

p656 <- ggplot(d656, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol=8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,25), xmax = c(21,42),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,21),xmax=c(4,25),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,25),xmax=c(21,42),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42))

####################################################################
# temperature periodogram

dt0.5[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
dt0.5[, .(id, uid), meta=T]
dt1[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
dt1[, .(id, uid), meta=T]
dt2[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
dt2[, .(id, uid), meta=T]
dt3[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
dt3[, .(id, uid), meta=T]
dt4[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
dt4[, .(id, uid), meta=T]
dt5[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
dt5[, .(id, uid), meta=T]


per_xsq_d0.5T <- periodogram(activity, 
                             dt0.5,
                             FUN = chi_sq_periodogram)

per_xsq_d0.5T <- find_peaks(per_xsq_d0.5T)

per_xsq_d1T <- periodogram(activity, 
                             dt1,
                             FUN = chi_sq_periodogram)

per_xsq_d1T <- find_peaks(per_xsq_d1T)

per_xsq_d2T <- periodogram(activity, 
                             dt2,
                             FUN = chi_sq_periodogram)

per_xsq_d2T <- find_peaks(per_xsq_d2T)

per_xsq_d3T <- periodogram(activity, 
                             dt3,
                             FUN = chi_sq_periodogram)

per_xsq_d3T <- find_peaks(per_xsq_d3T)

per_xsq_d4T <- periodogram(activity, 
                             dt4,
                             FUN = chi_sq_periodogram)

per_xsq_d4T <- find_peaks(per_xsq_d4T)

per_xsq_d5T <- periodogram(activity, 
                             dt5,
                             FUN = chi_sq_periodogram)

per_xsq_d5T <- find_peaks(per_xsq_d5T)


p0.5 <- ggperio(per_xsq_d0.5T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p1 <- ggperio(per_xsq_d1T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p2 <- ggperio(per_xsq_d2T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p3 <- ggperio(per_xsq_d3T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p4 <- ggperio(per_xsq_d4T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

p5 <- ggperio(per_xsq_d5T) + 
  facet_wrap(~ uid, ncol=8) +
  geom_line(aes(group = id)) +
  geom_peak(peak_rank=1,col = "black") +
  geom_line(aes(y = signif_threshold)) + 
  theme_classic()

####################################################################

# center of gravity for different filters and wavelength

temp = c(0.5,1,2,3,4,5)


ldf <- list() #create a list
listtxt <- dir(pattern = "*.txt") #creates the list of all the Texp. txt files in the directory (experimental groups)

for (k in 1:length(listtxt)){
  ldf[[k]] <- fread(listtxt[k])  #read all exp txt file
}

for (k in 1:length(listtxt)) {
  d <- ldf[[k]] #assign the first one to d
  d <- d[,-c(1, 5:9)] #delete not useable columns
  
  d <- d %>% dplyr::rename(Day="V2",minute="V3", DAM="V4",Light="V10",
                           "1"="V11","2"="V12","3"="V13","4"="V14","5"="V15",
                           "6"="V16","7"="V17","8"="V18","9"="V19","10"="V20",
                           "11"="V21","12"="V22","13"="V23","14"="V24","15"="V25",
                           "16"="V26","17"="V27","18"="V28","19"="V29","20"="V30",
                           "21"="V31","22"="V32","23"="V33","24"="V34","25"="V35",
                           "26"="V36","27"="V37","28"="V38","29"="V39","30"="V40",
                           "31"="V41","32"="V42")
  
  
  d <- melt(d, id.vars = c('Day', 'minute','DAM','Light'),variable.name = 'variable') #variable is FlyID
  d <- as.data.table(d)
  
  
  #experimental group
  x <- 24*60 #minutes in this T cycle
  y <- length(d$minute)/32 #how many lines for each animal
  bin <- rep(seq(from=1, to=ceiling(y/10), by=1), each=10)
  bin = bin[1:y]
  
  bin2 <- rep(bin, 32)
  
  d$bin <- bin2
  
  q <- subset(d, d$variable == 1)
  q$time <- rep(seq(1,x/10,1), each = 10, ceiling(y/x), length.out = nrow(q))
  
  d$time <- rep(q$time,32)
  
  daggregate <- d[, .(mean = mean(value, na.rm = T)), by = c('bin','time','variable')]
  
  
  day <- rep(seq(from=1, to=ceiling(y/x),by=1), each=x/10)
  z = length(daggregate$bin)/32
  day <- day[1:z]
  day <- rep(day, 32)
  
  daggregate$day <- day
  daggregate$temp = temp[k]
  
  #calculate phase angle of entrainment
  q1 <- subset(daggregate, daggregate$variable ==1)
  angle <- rep(seq(2*pi/(x/10), 2*pi, 2*pi/(x/10)), ceiling(y/x), length.out = nrow(q1))
  #i <- length(unique(daggregate$time))
  #angle <- angle[1:i]
  
  daggregate$angle <- rep(angle, 32)
  daggregate$sin <- sin(daggregate$angle)
  daggregate$cos <- cos(daggregate$angle)
  daggregate$sinactivity <- daggregate$mean*daggregate$sin
  daggregate$cosactivity <- daggregate$mean*daggregate$cos
  
  #just plotting and check data
  ggplot(subset(daggregate, daggregate$variable == '10'), aes(x=time, y=mean))+
    geom_line()+
    facet_grid(day~ .) 
  
  daggregate2 <- daggregate[, .(meansin=mean(sinactivity), meancos=mean(cosactivity), meanact=mean(mean)), by=c("variable","day")]
  daggregate2$Tcycle <- 21
  daggregate2$group <- "exprimental"
  daggregate2$temp <- temp[k]
  
  
  daggregate2$quard <- ifelse (daggregate2$meansin >0 ,
                               ifelse (daggregate2$meancos >0, 1,2),
                               ifelse(daggregate2$meancos <0, 3,4))
  daggregate2$tan <- daggregate2$meansin/daggregate2$meancos
  daggregate2$arctan <- atan(daggregate2$tan)
  daggregate2$angle <-  (daggregate2$arctan/(2*pi))*360
  
  daggregate2$correctedangle <- ifelse (daggregate2$quard == 1, daggregate2$angle,
                                        ifelse (daggregate2$quard == 2, daggregate2$angle+180,
                                                ifelse(daggregate2$quard ==3, daggregate2$angle+180,
                                                       daggregate2$angle+360 )))
  daggregate2$centralofgravity <- (daggregate2$correctedangle/360)*daggregate2$Tcycle
  daggregate2$centralofgravity_min <- (daggregate2$correctedangle/360)*(daggregate2$Tcycle*60/10)
  
  assign(paste("daggregate", k, sep="_"), daggregate)
  assign(paste("daggregate2", k, sep="_"), daggregate2)
}

dd <- rbind(daggregate2_1,daggregate2_2,daggregate2_3,daggregate2_4,
            daggregate2_5,daggregate2_6)
ddd<- rbind(daggregate_1,daggregate_2,daggregate_3,daggregate_4,
            daggregate_5,daggregate_6)


dd_2 <- dd
dd_2$centralofgravity <- dd_2$centralofgravity + dd_2$Tcycle

dd_2a <- bind_rows(dd,dd_2)

d0.5 <- subset(dd_2a, temp == 0.5)
d1 <- subset(dd_2a, temp == 1)
d2 <- subset(dd_2a, temp == 2)
d3 <- subset(dd_2a, temp == 3)
d4 <- subset(dd_2a, temp == 4)
d5 <- subset(dd_2a, temp == 5)


p0.5 <- ggplot(d0.5, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))

p1 <- ggplot(d1, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))

p2 <- ggplot(d2, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))

p3 <- ggplot(d3, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))

p4 <- ggplot(d4, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))

p5 <- ggplot(d5, aes(x=centralofgravity, y= day))+
  geom_point()+
  facet_wrap(~ variable, ncol = 8)+
  scale_y_reverse(breaks= c(-2, 0,2,4,6, 8,10,12,14,16,18))+
  theme_classic()+
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = c(4,28), xmax = c(24,48),
           ymin = 0, ymax = Inf) +
  annotate("rect", colour = "black", fill="white",
           xmin=c(0,24),xmax=c(4,28),ymin=0,ymax=-0.5, alpha=0.9)+
  annotate("rect", colour = "black", fill="black",
           xmin=c(4,28),xmax=c(24,48),ymin=0,ymax=-0.5, alpha=0.9)+
  scale_x_continuous(breaks = c(0,8,16,24,32,40,48))


####################################################################

#create a new powerpoint document

j = list(p455, p470, p505, p528, p566, p590, p617, p625, p656)
j = list(p0.5,p1,p2,p3,p4,p5)

for (k in 1:9) {
  
doc <- read_pptx("last experiment.pptx")
doc <- add_slide(doc, 'Title and Content', 'Office Theme')

#doc <- ph_with(doc, pp , location = ph_location(width=6,height=12))#for too big plot
# 
fig <- dml(ggobj =j[[k]])
# #add the plot
doc <- ph_with(doc, fig, location = ph_location(width=12,height=12))
# #write the document to a file

print(doc, target ="last experiment.pptx")
}

