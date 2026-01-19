# This script is as part of the manuscript Circadian entrainment to red-light Zeitgebers and action spectrum for entrainment in the jewel wasp Nasonia vitripennis
# In this script, we modeled the effect of wavelength and T-cycle on entrainment and the circadian action spectrum of entrainment.

rm(list=ls()) #clear all objects

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

###prepare actograms########
#figure 1
DATA_DIR <-"/path/to/data/directory"
list.files(DATA_DIR, pattern = ".txt|*.csv")
setwd(DATA_DIR)

metadata <- fread("metadata.csv")

metadataS <- subset(metadata, metadata$group=="experimental")
metadata1 <- metadataS[c(1:96),]
metadata2 <- metadataS[c(97:192),]

metadata1 <- metadata1 %>% mutate(start_datetime = "2021-06-16 16:00:00", stop_datetime = "2021-06-30 16:00:00")
metadata2 <- metadata2 %>% mutate(start_datetime = "2021-07-01 17:00:00", stop_datetime = "2021-07-15 17:00:00")

metadata1 <- link_dam_metadata(metadata1,result_dir = DATA_DIR)
metadata2 <- link_dam_metadata(metadata2,result_dir = DATA_DIR)
dt1 <- load_dam(metadata1)
dt2 <- load_dam(metadata2)

dt1[, moving := activity >0]
dt2[, moving := activity >0]

d24T <- dt1[xmv(photoperiod)==24]
d23T <- dt1[xmv(photoperiod)==23]
d22T <- dt1[xmv(photoperiod)==22]
d21T <- dt2[xmv(photoperiod)==21]
d20T <- dt2[xmv(photoperiod)==20]
d19T <- dt2[xmv(photoperiod)==19]

d24T[,uid:=1 : .N, meta=T] #make a new variable uid for each individual
d24T[, .(id, uid), meta=T]
d23T[,uid:=1 : .N, meta=T] 
d23T[, .(id, uid), meta=T]
d22T[,uid:=1 : .N, meta=T] 
d22T[, .(id, uid), meta=T]
d21T[,uid:=1 : .N, meta=T] 
d21T[, .(id, uid), meta=T]
d20T[,uid:=1 : .N, meta=T] 
d20T[, .(id, uid), meta=T]
d19T[,uid:=1 : .N, meta=T] 
d19T[, .(id, uid), meta=T]


#select individual to plot
T24 <- ggetho(d24T[id %in% c("2021-06-16 16:00:00|Monitor5_24Texp.txt|01",
                             "2021-06-16 16:00:00|Monitor5_24Texp.txt|21",
                             "2021-06-16 16:00:00|Monitor5_24Texp.txt|24")], aes(x=t, z=moving), multiplot = 2,
              multiplot_period = hours(24))+
  stat_bar_tile_etho(fill="black") +
  facet_wrap(~ uid, ncol=3) +
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  stat_ld_annotations(ypos="top",l_duration = hours(4),period = hours(24),ld_colours = c("white","black")) +
  stat_ld_annotations(l_duration= hours(4),period = hours(24 ),height=1, alpha=0.01, outline=NA)+
  #stat_ld_annotations(l_duration= hours(0),ypos="bottom",height=(1-8/18), alpha=0.01, outline=NA)+
  theme(strip.text.x = element_text(size=0),axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text=element_text(size=11)) +
  scale_y_discrete(breaks=c("0", "2", "4", "6","8","10","12","14","16","18") )+
  ggetho::scale_x_hours(breaks = hours(c(0,8,16,24,32,40,48)))
T24


T23 <- ggetho(d23T[id %in% c("2021-06-16 16:00:00|Monitor7_23Texp.txt|16",
                             "2021-06-16 16:00:00|Monitor7_23Texp.txt|18",
                             "2021-06-16 16:00:00|Monitor7_23Texp.txt|23")], aes(x=t, z=moving), multiplot = 2,
              multiplot_period = hours(23))+
  stat_bar_tile_etho(fill="black") +
  facet_wrap(~ uid, ncol=3) +
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  stat_ld_annotations(ypos="top",l_duration = hours(4),period = hours(23),ld_colours = c("white","black")) +
  stat_ld_annotations(l_duration= hours(4),period = hours(23),height=1, alpha=0.01, outline=NA)+
  #stat_ld_annotations(l_duration= hours(0),ypos="bottom",height=(1-8/18), alpha=0.01, outline=NA)+
  theme(strip.text.x = element_text(size=0),axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text=element_text(size=11)) +
  scale_y_discrete(breaks=c("0", "2", "4", "6","8","10","12","14") )+
  ggetho::scale_x_hours(breaks = hours(c(0,8,16,24,32,40,46)))
T23
T22 <- ggetho(d22T[id %in% c("2021-06-16 16:00:00|Monitor9_22Texp.txt|01",
                             "2021-06-16 16:00:00|Monitor9_22Texp.txt|04",
                             "2021-06-16 16:00:00|Monitor9_22Texp.txt|25")], aes(x=t, z=moving), multiplot = 2,
              multiplot_period = hours(22))+
  stat_bar_tile_etho(fill="black") +
  facet_wrap(~ uid, ncol=3) +
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  stat_ld_annotations(ypos="top",l_duration = hours(4),period = hours(22),ld_colours = c("white","black")) +
  stat_ld_annotations(l_duration= hours(4),period = hours(22),height=1, alpha=0.01, outline=NA)+
  #stat_ld_annotations(l_duration= hours(0),ypos="bottom",height=(1-8/18), alpha=0.01, outline=NA)+
  theme(strip.text.x = element_text(size=0),axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text=element_text(size=11)) +
  scale_y_discrete(breaks=c("0", "2", "4", "6","8","10","12","14") )+
  ggetho::scale_x_hours(breaks = hours(c(0,7,14,21,28,35,44)))
T22
T21 <- ggetho(d21T[id %in% c("2021-07-01 17:00:00|Monitor5_21Texp.txt|08",
                             "2021-07-01 17:00:00|Monitor5_21Texp.txt|22",
                             "2021-07-01 17:00:00|Monitor5_21Texp.txt|29")], aes(x=t, z=moving), multiplot = 2,
              multiplot_period = hours(21))+
  stat_bar_tile_etho(fill="black") +
  facet_wrap(~ uid, ncol=3) +
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  stat_ld_annotations(ypos="top",l_duration = hours(4),period = hours(21),ld_colours = c("white","black")) +
  stat_ld_annotations(l_duration= hours(4),period = hours(21),height=1, alpha=0.01, outline=NA)+
  #stat_ld_annotations(l_duration= hours(0),ypos="bottom",height=(1-8/18), alpha=0.01, outline=NA)+
  theme(strip.text.x = element_text(size=0),axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text=element_text(size=11)) +
  scale_y_discrete(breaks=c("0", "2", "4", "6","8","10","12","14","16","18","20") )+
  ggetho::scale_x_hours(breaks = hours(c(0,7,14,21,28,35,42)))
T21
T20 <- ggetho(d20T[id %in% c("2021-07-01 17:00:00|Monitor7_20Texp.txt|18",
                             "2021-07-01 17:00:00|Monitor7_20Texp.txt|32",
                             "2021-07-01 17:00:00|Monitor7_20Texp.txt|06")], aes(x=t, z=moving), multiplot = 2,
              multiplot_period = hours(20))+
  stat_bar_tile_etho(fill="black") +
  facet_wrap(~ uid, ncol=3) +
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  stat_ld_annotations(ypos="top",l_duration = hours(4),period = hours(20),ld_colours = c("white","black")) +
  stat_ld_annotations(l_duration= hours(4),period = hours(20),height=1, alpha=0.01, outline=NA)+
  theme(strip.text.x = element_text(size=0),axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text=element_text(size=11)) +
  scale_y_discrete(breaks=c("0", "2", "4", "6","8","10","12","14") )+
  ggetho::scale_x_hours(breaks = hours(c(0,8,16,24,32,40)))
T20
T19 <- ggetho(d19T[id %in% c("2021-07-01 17:00:00|Monitor9_19Texp.txt|07",
                             "2021-07-01 17:00:00|Monitor9_19Texp.txt|24",
                             "2021-07-01 17:00:00|Monitor9_19Texp.txt|28")], aes(x=t, z=moving), multiplot = 2,
              multiplot_period = hours(19))+
  stat_bar_tile_etho(fill="black") +
  facet_wrap(~ uid, ncol=3) +
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  stat_ld_annotations(ypos="top",l_duration = hours(4),period = hours(19),ld_colours = c("white","black")) +
  stat_ld_annotations(l_duration= hours(4),period = hours(19),height=1, alpha=0.01, outline=NA)+
  theme(strip.text.x = element_text(size=0),axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text=element_text(size=11)) +
  scale_y_discrete(breaks=c("0", "2", "4", "6","8","10","12","14") )+
  ggetho::scale_x_hours(breaks = hours(c(0,8,16,24,32,38)))
T19

p <- ggarrange(T24,T23,T22,T21,T20,T19,
               ncol=1,
               nrow=6,
               align = "h",
               common.legend=T,
               legend = 'right')

p

pp <- annotate_figure(p, left = text_grob("days",rot=90, size = 14, face = "bold"),
                      bottom = text_grob("time (ZT, hrs)", size = 14, face = "bold"))

pp

#############################################
#figure 2
####phase angle of entrainment calculation###

Tcycle <- c(19,20,21,22,23,24)

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
  
  
  d <- melt(d, id.vars = c('Day', 'minute','DAM','Light'),variable.name = 'FlyID')
  d <- as.data.table(d)
  d <- d %>% mutate(Tcycle=Tcycle[k], group="experimental",wls="656 nm")
  
  #experimental group
  x <- Tcycle[k]*60 #minutes in this T cycle
  y <- length(d$minute)/32 #how many lines for each animal
  Seq1 <- rep(seq(from=1, to=x, by=1), y/x) #create seq1
  z <- length(Seq1) #check how long seq1 is
  i <- y-z # how many more numbers to add
  Seq2 <- c(Seq1, seq(1, i,1)) #create seq2
  Tminute <- rep(Seq2, each=1, 32) #number for the datatable
  
  d$Tminute <- Tminute
  
  daggregate <- d[, .(mean = mean(value, na.rm = T)), by = c('Tminute','FlyID')]
  str(daggregate)
  
  bin <- rep(seq(1, length(daggregate$Tminute)/320,1), each=10, 32)
  daggregate$bin <- bin
  daggregate2 <- daggregate[, .(meanactivity = mean(mean)), by = c('bin','FlyID')]
  
  daggregate2$Tcycle <- Tcycle[k]
  daggregate2$group <- "exprimental"
  daggregate2$wls <- "656 nm"
  
  
  #calculate phase angle of entrainment
  x <- Tcycle[k]*60/10 #how many 10mins in this tcycle
  daggregate2$angle <- rep(seq(2*pi/x, 2*pi, 2*pi/x),32)
  daggregate2$sin <- sin(daggregate2$angle)
  daggregate2$cos <- cos(daggregate2$angle)
  daggregate2$sinactivity <- daggregate2$meanactivity*daggregate2$sin
  daggregate2$cosactivity <- daggregate2$meanactivity*daggregate2$cos
  
  daggregate3 <- daggregate2[, .(meansin=mean(sinactivity), meancos=mean(cosactivity)), by=c("FlyID")]
  daggregate3$Tcycle <- Tcycle[k]
  daggregate3$group <- "exprimental"
  daggregate3$wls <- "656 nm"
  
  
  daggregate3$quard <- ifelse (daggregate3$meansin >0 ,
                               ifelse (daggregate3$meancos >0, 1,2),
                               ifelse(daggregate3$meancos <0, 3,4))
  daggregate3$tan <- daggregate3$meansin/daggregate3$meancos
  daggregate3$arctan <- atan(daggregate3$tan)
  daggregate3$angle <-  (daggregate3$arctan/(2*pi))*360
  daggregate3$correctedangle <- ifelse (daggregate3$quard == 1, daggregate3$angle,
                                        ifelse (daggregate3$quard == 2, daggregate3$angle+180,
                                                ifelse(daggregate3$quard ==3, daggregate3$angle+180,
                                                       daggregate3$angle+360 )))
  daggregate3$centralofgravity <- (daggregate3$correctedangle/360)*daggregate3$Tcycle
  daggregate3$centralofgravity_min <- (daggregate3$correctedangle/360)*(daggregate3$Tcycle*60/10)
  
  assign(paste("daggregate2", k, sep="_"), daggregate2)
  assign(paste("daggregate3", k, sep="_"), daggregate3)
}

d656 <- rbind(daggregate2_1,daggregate2_2,daggregate2_3,daggregate2_4,
               daggregate2_5,daggregate2_6)
d656_2 <- rbind(daggregate3_1,daggregate3_2,daggregate3_3,daggregate3_4,
                 daggregate3_5,daggregate3_6)

#fwrite(daggregate3,"d590_21T_phaseangleofentrainment.csv")

#calculate 21 T cycle 590nm with the new data
d21_590 <- fread("Monitor10.1.txt")

#delete the first few days of non entrain
d21_590 <- d21_590[-c(0:5760),]
#fwrite(d21_590, "Monitor10.1_removefirst4days.txt")

#phase angle of entrainment on a circular scale
doveralldata <- fread("redlightentrainment_phaseangleofentrainment.csv")

doveralldata1 <- subset(doveralldata, Entrainment=="entrain")
doveralldata1$wavelength <- as.factor(doveralldata1$wavelength)

#making phase angle of entrainment double plot
doveralldata2 <- doveralldata1
doveralldata2.1 <- doveralldata1
doveralldata2$correctedangle <- doveralldata2$correctedangle+360
doveralldata2.1$correctedangle <- doveralldata2.1$correctedangle-360
doveralldata3 <- rbind(doveralldata1,doveralldata2, doveralldata2.1)

plot1 <- ggplot(doveralldata3, aes(x=Tcycle, y=correctedangle))+
  geom_point(aes(shape=wavelength, colour=wavelength))+
  # stat_summary(fun.data=mean_cl_normal) + 
  #geom_smooth(method='lm', formula= y~x, colour="black", size=0.5)+
  facet_wrap(wavelength~.)+
  theme_classic()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  labs(title="Phase angle of entrainment",x="T cycle (h)", y="phase angle") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  theme(axis.text=element_text(size=8)) +
  geom_hline(aes(yintercept=360), linetype="dashed")

plot1
#select data based on double plot to predict linear regression
list <-  ifelse(doveralldata3$wavelength==590,
                ifelse(doveralldata3$Tcycle==19 & doveralldata3$correctedangle <360 &doveralldata3$correctedangle >0, 99, 
                       ifelse(doveralldata3$Tcycle==20 & doveralldata3$correctedangle <360&doveralldata3$correctedangle >0, 99, 
                              ifelse(doveralldata3$Tcycle==21 & doveralldata3$correctedangle <360&doveralldata3$correctedangle >0, 99, 
                                     ifelse(doveralldata3$Tcycle==22 & doveralldata3$correctedangle >0 & doveralldata3$correctedangle <360, 99,
                                            ifelse(doveralldata3$Tcycle==23 & doveralldata3$correctedangle >250 & doveralldata3$correctedangle <500, 99,
                                                   ifelse(doveralldata3$Tcycle==24 & doveralldata3$correctedangle >360, 99, 0)))))),
                ifelse(doveralldata3$wavelength==625,
                       ifelse(doveralldata3$Tcycle==20 & doveralldata3$correctedangle <360&doveralldata3$correctedangle >0, 99,
                              ifelse(doveralldata3$Tcycle==21 & doveralldata3$correctedangle <360&doveralldata3$correctedangle >0, 99,
                                     ifelse(doveralldata3$Tcycle==22 & doveralldata3$correctedangle >100 & doveralldata3$correctedangle <500, 99,
                                            ifelse(doveralldata3$Tcycle==23 & doveralldata3$correctedangle >250 & doveralldata3$correctedangle <500,99,
                                                   ifelse(doveralldata3$Tcycle==24 & doveralldata3$correctedangle >250, 99, 0))))),
                       ifelse(doveralldata3$wavelength==656,
                              ifelse(doveralldata3$Tcycle==19 & doveralldata3$correctedangle <100 & doveralldata3$correctedangle >-250, 99,
                                     ifelse(doveralldata3$Tcycle==20 & doveralldata3$correctedangle <480&doveralldata3$correctedangle >0, 99,
                                            ifelse(doveralldata3$Tcycle==21 & doveralldata3$correctedangle <500&doveralldata3$correctedangle >0, 99,
                                                   ifelse(doveralldata3$Tcycle==22 & doveralldata3$correctedangle <500&doveralldata3$correctedangle >0, 99,
                                                          ifelse(doveralldata3$Tcycle==23 & doveralldata3$correctedangle >200 & doveralldata3$correctedangle <500,99,
                                                                 ifelse(doveralldata3$Tcycle==24 & doveralldata3$correctedangle >360, 99, 0)))))),0)))

doveralldata3$list <- list   
doveralldata3.1 <- subset(doveralldata3, list==99)

# fwrite(doveralldata3.1, "redlightentrainment_modelpredictiondata.csv")

#central of gravity plot
dCOG <- fread("redlightentrainment_modelpredictiondata.csv")
dCOG$wavelength <- as.factor(dCOG$wavelength)

plot1 <- ggplot(dCOG, aes(x=Tcycle, y=correctedangle))+
  geom_point(aes(shape=wavelength, colour=wavelength))+
  # stat_summary(fun.data=mean_cl_normal) + 
  #geom_smooth(method='lm', formula= y~x, colour="black", size=0.5)+
  facet_wrap(wavelength~.)+
  theme_classic()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  labs(title="Phase angle of entrainment",x="T cycle (hrs)", y="phase angle") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  theme(axis.text=element_text(size=8)) +
  geom_hline(aes(yintercept=360), linetype="dashed")+
  scale_y_continuous( breaks=c(-240,-120, 0, 120, 240, 360,480,600,720) )
plot1

dCOG1 <- dCOG[!is.na(dCOG$correctedangle)] 
dCOG1 <- dCOG1 %>% dplyr::rename(phaseangle= correctedangle)

#first model
m_0 <- lm(phaseangle ~ Tcycle*wavelength, data = dCOG1) #Tcycle+wavelength+group+Tcycle:wavelength
drop1(m_0, test='Chisq') #stepwise model deletion p=0.1112
m_1 <- update(m_0, ~.-Tcycle:wavelength)
drop1(m_1, test='Chisq') 

anova(m_0,m_1) 
AIC(m_0,m_1)

anova(m_0)
m_1$coefficients
m.lst <- lstrends(m_1, "wavelength", var="Tcycle")
m.lst
pairs(m.lst)

# anova(m2,m3)
summary(m_1)
# AIC(m0,m1)

#check posthoc test, check differences between factor levels, only when anova is significant, compare if there's difference between the slop
emtrends(m_1, c("wavelength"), var = "Tcycle")


pred.d <- with(dCOG1, data.frame(Tcycle = dCOG1$Tcycle,
                                           wavelength = dCOG1$wavelength))
pred.d$phaseangle <- predict(m_1, newdata=pred.d, type='response')

# The palette for each wavelength:
cbPalette <- c( 
  "#FDDA0D",#590# 
  "#FF0000", #625# 
  "#900C3F" ) #656#

p1 <- ggplot(dCOG1, aes(x=Tcycle, y=phaseangle))+
  geom_point(aes(shape=wavelength, colour=wavelength),size=2)+
  scale_colour_manual(values=cbPalette)+
  #geom_smooth(method='lm', formula= y~x, colour="black", size=0.5)+
  geom_line(data=pred.d, aes(x=Tcycle, y=phaseangle),size=1, linetype="dashed")+
  scale_shape_manual(values=c(17, 15,18))+
  facet_wrap(wavelength~.)+
  theme_classic()+
  theme(panel.spacing = unit(2, "lines"))+
#  theme(strip.background = element_blank(), strip.placement = "outside", strip.text.x = element_text(size=12))+
  labs(x="T cycle (h)", y="phase angle of entrainment (degree)") +
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11)) +
  theme(strip.text.x = element_blank())+
  theme(legend.position="none")

p1

###percentage of entrainment plot####
doverall <- fread("red light entrainment 3wls 6Tcycles.csv")
daggregate1 <- fread("red light entrainment percentage for graph.csv")
daggregate1$wavelength <- as.factor(daggregate1$wavelength)

doverall$count <- 1
doverallave <- doverall %>% group_by(Tcycle, wavelength, Entrainment) %>% 
                summarise(count= n())
doverallave$percentage <- doverallave$count/32
#fwrite(doverallave, "red light entrainment average 3wls 6tcycles.csv")

Entrained <- subset(daggregate1,daggregate1$Entrainment=="entrain")


m1 <- glm(Percentage ~ wavelength*Tcycle, family="quasibinomial", data=Entrained)
m2 <- glm(Percentage ~ wavelength+Tcycle, family="quasibinomial", data=Entrained)
summary(m1)
library(car)
Anova(m1)
summary(m2)
Anova(m2)

residualPlots(m1)


# The palette for each wavelength:
cbPalette <- c( 
  "#FDDA0D",#590# 
  "#FF0000", #625# 
  "#900C3F" ) #656#
p2 <-  ggplot(Entrained, aes(x=Tcycle, y=Percentage,colour=wavelength, shape=wavelength))+
  geom_point(size=3)+
  scale_shape_manual(values=c(17, 15,18))+
  scale_colour_manual(values=cbPalette)+
 # geom_line(aes(linetype=wavelength), size=1)+
  #scale_linetype_manual(values=c("dotdash", "solid"))+
  #facet_wrap(wavelength~., ncol=3)+
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se=F, aes(linetype=wavelength)) +
  theme_classic()+
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="T cycle (h)", y="percentage of entrainment") +
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11)) +
  theme(legend.text = element_text(size=16), legend.title = element_text(size=16),legend.key.size = unit(2, 'cm'))
p2



#test if dead or arrthymic animal was affected by Tcycle or wavelength
arrhythmic <- subset(doverallave, Entrainment == "arrhythmic")
m_0 <- glm(percentage ~ wavelength*Tcycle, family="quasibinomial", data=arrhythmic)
summary(m_0)
Anova(m_0)

dead <- subset(doverallave, Entrainment == "dead")
m_1 <- glm(percentage ~ wavelength*Tcycle, family="quasibinomial", data=dead)
summary(m_1)
Anova(m_1)

##figure 2c
spectrometer <- fread("red light light spectrum.csv")

spectrometer$relative590 <- spectrometer$`590nm`/mean(spectrometer$`590nm`)
spectrometer$relative625 <- spectrometer$`625nm`/mean(spectrometer$`625nm`)
spectrometer$relative656 <- spectrometer$`656nm`/mean(spectrometer$`656nm`)

spectrometer <- spectrometer[,-c(2:4)]
spectrometer2 <- melt(spectrometer, id.vars = c('wavelength'),variable.name = 'LED') #transpose the dataset for further plotting

#plot for spectrometer data
p3 <-  ggplot(data = spectrometer2, aes(x = wavelength, y = value, group=LED, fill=LED))+
  geom_line(size=0.5,aes(colour=LED))+
  scale_colour_manual(values=cbPalette)+
  geom_ribbon(aes(x=wavelength,ymax=value),ymin=0) +
  scale_fill_manual(name='', values=c( "#FDDA0D", "#FF0000", "#900C3F"))+
  labs(x="wavelength (nm)", y="relative intensity") +
  theme_classic()+xlim(500,800)+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11)) +
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) #change legend text fo
p3

p <- ggarrange(
               p2,p3,
               ncol=2,
               nrow=1,
               common.legend=T,
               legend = 'right')

p

pp <- ggarrange(p1,p,
                ncol=1,nrow=2)
pp


#############################################
##figure 3
###graphs for light entrainment
led_entrain <- fread("different LED light experiment/overall data/different LED experiment entrain percentage.csv")
led_entrain$count <-ifelse(led_entrain$Actogram=="entrained" |led_entrain$Actogram=="freerunning", 1,0) 
led_entrain2 <- subset(led_entrain, led_entrain$count==1)
led_entrain2$count <- ifelse(led_entrain2$Actogram=="entrained", 1,0) 

led_entrain <- fread("different LED experiment percentage for graphs.csv")
led_entrain$irradiance <- log10(led_entrain$`lightintensity`)
led_entrain$wavelength <- as.factor(led_entrain$wavelength)
led_entrain2 <- subset(led_entrain, led_entrain$Actogram=="Entrain")
led_entrain2 <- subset(led_entrain2, led_entrain2$wavelength!="540 nm")

# The palette with grey:
cbPalette <- c("#00008b",#455# 
               "#009dff",#470nm#
               "#7FB383", #505# 
               "#3EA055",#528#
               # "008000",#540#
               "#85BB65",#566# 
               "#FDDA0D",#590# 
               "#FF5733",#617# 
               "#FF0000", #625# 
               "#900C3F" ) #656#

plot <-  ggplot(led_entrain2, aes(x=irradiance, y=Percentage, group=wavelength))+
  geom_point(size=3, aes(x=irradiance, y=Percentage,colour=wavelength))+
  geom_line(aes(colour=wavelength), size=1)+
  theme_classic()+
  facet_wrap(~wavelength)+
  scale_colour_manual(values=cbPalette)+
  labs(x="light intensity (photons*m-2*s-1)", y="percentage of entrainment") +
  scale_x_continuous(breaks = c(12,12.5,13,13.5,14,14.5,15,15.5,16))+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=12)) +
  ylim(0,1)+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size
plot


m1 <- drm(Percentage ~ irradiance, data=led_entrain2, curveid = wavelength,
          fct = LL.2(names = c("Slope", "ED50")), 
          type="continuous")
summary(m1)
#modelFit(m1)
plot(m1)

led_entrain3 <- subset(led_entrain2, led_entrain2$wavelength!="505 nm")
m2 <- drm(Percentage ~ irradiance, data=led_entrain3, curveid = wavelength,
          fct = LL.2(names = c("Slope", "ED50")), 
          type="continuous")

summary(m2)
plot(m2)

# new data with predictions with m1
pred.d <- with(led_entrain2,data.frame(irradiance=rep(seq(10,18, 0.5),10),
                                       wavelength=rep(levels(led_entrain2$wavelength), each =17)))
pred.d <- subset(pred.d, pred.d$wavelength!="540 nm")
pm <- predict(m1, newdata=pred.d, interval="confidence") 
pred.d$p <- pm[,1]
pred.d$pmin <- pm[,2]
pred.d$pmax <- pm[,3]

# new data with predictions with m2
pred.d2 <- with(led_entrain3,data.frame(irradiance=rep(seq(10,18, 0.5),10),
                                        wavelength=rep(levels(led_entrain3$wavelength), each =17)))
pred.d2 <- subset(pred.d2, pred.d2$wavelength!="540 nm")
pred.d2 <- subset(pred.d2, pred.d2$wavelength!="505 nm")
#pred.d2 <- subset(pred.d2, pred.d2$wavelength!=505)
pm <- predict(m2, newdata=pred.d2, interval="confidence") 
pred.d2$p <- pm[,1]
pred.d2$pmin <- pm[,2]
pred.d2$pmax <- pm[,3]


p1 <- ggplot(led_entrain2,aes(x=irradiance,y=Percentage))+
  geom_point(aes(colour=wavelength))+
  facet_wrap(wavelength~.)+
  theme_classic()+
  scale_colour_manual(values=cbPalette)+
  #geom_ribbon(data=pred.d, aes(x=irradiance, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=pred.d, aes(x=irradiance, y=p,colour=wavelength)) +
  labs(x="log10 of light intensity (photons*m-2*s-1)", y="percentage of entrainment")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11)) +
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) #change legend text fo

p1

p2 <- ggplot(led_entrain3,aes(x=irradiance,y=Percentage))+
  geom_point(aes(colour=wavelength))+
  facet_wrap(wavelength~.)+
  theme_classic()+
  scale_colour_manual(values=cbPalette)+
  # geom_ribbon(data=pred.d2, aes(x=irradiance, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=pred.d2, aes(x=irradiance, y=p,colour=wavelength)) +
  labs(x="log10 of light intensity (photons*m-2*s-1)", y="percentage of entrainment")#+
#coord_trans(x="log")

p2

ed50 <- as.data.frame(ED(m2, 50, interval = "delta"))
ed50$wavelength <- c(455,470,528,566,590,617,625,656)
ed50$ED <- 50
ed50$sensitivity <- 1/(ed50$Estimate)
max = max(ed50$sensitivity)
ed50$normal <- ed50$sensitivity/max

ed30 <- as.data.frame(ED(m2, 30, interval = "delta"))
ed30$wavelength <- c(455,470,528,566,590,617,625,656)
ed30$ED <- 30
ed30$sensitivity <- 1/(ed30$Estimate)
max = max(ed30$sensitivity)
ed30$normal <- ed30$sensitivity/max

ed80 <- as.data.frame(ED(m2, 80, interval = "delta"))
ed80$wavelength <- c(455,470,528,566,590,617,625,656)
ed80$ED <- 80
ed80$sensitivity <- 1/(ed80$Estimate)
max = max(ed80$sensitivity)
ed80$normal <- ed80$sensitivity/max

ED <- rbind(ed30,ed50,ed80)
ED$ED <- as.factor(ED$ED)

cbPalette1 <- c("lightblue2", "skyblue3","steelblue4") 
p3 <- ggplot(ED, aes(x=wavelength, y=normal, shape=ED,color=ED))+
  scale_colour_manual(values=cbPalette1)+
  geom_point(size=3)+
  geom_line(size=1.05)+
  theme_classic()+
  labs(title="entrainment behaviour",x="wavelength (nm)", y="relative sensitivity")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11), plot.title = element_text(hjust = 0.5,size = 14, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))+theme(legend.position = "none") #change legend text fo

p3

#getwd()
setwd("C:/Yifan/PhD/lab book/red light entrainment")

spectrometer <- fread("all LB LED light spectrum.csv")

spectrometer$"455 nm" <- spectrometer$`455nm`/mean(spectrometer$`455nm`)
spectrometer$"470 nm" <- spectrometer$`470nm`/mean(spectrometer$`470nm`)
spectrometer$"505 nm" <- spectrometer$`505nm`/mean(spectrometer$`505nm`)
spectrometer$"528 nm" <- spectrometer$`528nm`/mean(spectrometer$`528nm`)
spectrometer$"566 nm" <- spectrometer$`566nm`/mean(spectrometer$`566nm`)
spectrometer$"590 nm" <- spectrometer$`590nm`/mean(spectrometer$`590nm`)
spectrometer$"617 nm" <- spectrometer$`617nm`/mean(spectrometer$`617nm`)
spectrometer$"625 nm" <- spectrometer$`625nm`/mean(spectrometer$`625nm`)
spectrometer$"656 nm" <- spectrometer$`656nm`/mean(spectrometer$`656nm`)

spectrometer <- spectrometer[,-c(2:11)]

spectrometer2 <- melt(spectrometer, id.vars = c('wavelength'),variable.name = 'LED') #transpose the dataset for further plotting

spectrometer2$count <- rep(seq(1,9, by = 1), each =401)

# The palette for each wavelength:
cbPalette <- c("#00008b",#455# 
               "#009dff",#470nm#
               "#7FB383", #505# 
               "#3EA055",#528#
               "#85BB65",#566# 
               "#FDDA0D",#590# 
               "#FF5733",#617# 
               "#FF0000", #625# 
               "#900C3F" ) #656#
p4 <-  ggplot(data = spectrometer2, aes(x = wavelength, y = value,  group=LED, fill=LED))+
  geom_line(size=0.5,aes(colour=LED))+
  scale_colour_manual(values=cbPalette)+
  geom_ribbon(aes(x=wavelength,ymax=value),ymin=0) +
  scale_fill_manual(name='', values=c("#00008b",#455# 
                                      "#009dff",#470nm#
                                      "#7FB383", #505# 
                                      "#3EA055",#528#
                                      "#85BB65",#566# 
                                      "#FDDA0D",#590# 
                                      "#FF5733",#617# 
                                      "#FF0000", #625# 
                                      "#900C3F" ) )+
  labs(x="wavelength (nm)", y="relative intensity") +
  theme_classic()+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=12)) +
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.position = "none")
p4

###########################
#input electrophysiology dose response curves

ele <- fread("Nasonia9wtFERG.csv")

ele$wavelength <- as.factor(ele$wavelength)

#dose response curves
all(is.na(ele$light_intensity))
all(is.na(ele$value))

ele$irradiance <- 10^ele$light_intensity

m3 <- drm(value ~ irradiance, curveid = wavelength, data = ele,
          fct = LL.2(names = c("Hill slopes","ED50")))
summary(m3)

plot(m3)

# new data with predictions
pred.d3 <- with(ele,data.frame(light_intensity=rep(seq(-5,0, 0.5),10),
                               wavelength=rep(levels(ele$wavelength), each =11)))
pred.d3$irradiance <- 10^pred.d3$light_intensity

pm <- predict(m3, newdata=pred.d3, interval="confidence") 
pred.d3$p <- pm[,1]
pred.d3$pmin <- pm[,2]
pred.d3$pmax <- pm[,3]

pred.d3 <- subset(pred.d3, pred.d3$wavelength!="540")
ele1 <- subset(ele, ele$wavelength!="540")


p5 <- ggplot(ele1,aes(x=log10(irradiance),y=value, colour=wavelength))+
  geom_point()+
  #facet_wrap(wavelength~.)+
  theme_classic()+
  scale_colour_manual(values=cbPalette)+
  # geom_ribbon(data=pred.d3, aes(x=light_intensity, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=pred.d3, aes(x=log10(irradiance), y=p)) +
  labs(x="normalized light intensity (log10(I/Imax))", y="electrophysiology response")+theme_classic()+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none")

p5


ed50 <- as.data.frame(ED(m3, 50, interval = "delta"))
ed50$wavelength <- c(455,470,505,528,540,566,590,617,625,656)
ed50$ED <- 50
ed50$sensitivity <- 1/ed50$Estimate
max = max(ed50$sensitivity)
ed50$normal <- ed50$sensitivity/max


ed30 <- as.data.frame(ED(m3, 30, interval = "delta"))
ed30$wavelength <- c(455,470,505,528,540,566,590,617,625,656)
ed30$ED <- 30
ed30$sensitivity <- 1/(ed30$Estimate)
max = max(ed30$sensitivity)
ed30$normal <- ed30$sensitivity/max

ed80 <- as.data.frame(ED(m3, 80, interval = "delta"))
ed80$wavelength <- c(455,470,505,528,540,566,590,617,625,656)
ed80$ED <- 80
ed80$sensitivity <- 1/(ed80$Estimate)
max = max(ed80$sensitivity)
ed80$normal <- ed80$sensitivity/max

ED2 <- rbind(ed30,ed50,ed80)
ED2$ED <- as.factor(ED2$ED)

p6 <- ggplot(ED2, aes(x=wavelength, y=normal, shape=ED,colour=ED))+
  geom_point(size=3)+
  geom_line(size=1.05)+
  scale_colour_manual(values=cbPalette1)+
  theme_classic()+
  labs(title="electrophysiology response", x="wavelength (nm)", y="relative sensitivity")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) +theme(legend.position = "none")#change legend text fo+

p6
###
p <- ggarrange(p4,p5,p3,p6, 
               ncol=2,nrow=2,
               align = "hv")
p

pp <- ggarrange(p1, p,
                ncol=1, nrow =2,
                common.legend = T, legend="right",
                heights = c(1,1))

pp

####supplementary figure of temperature entrainment
temp_entrain <- fread("temperature experiment entrain percentage.csv")
temp_entrain$count <-ifelse(temp_entrain$Entrain=="entrained" |temp_entrain$Entrain=="freerunning", 1,0) 
#   
temp_entrain2 <- subset(temp_entrain, temp_entrain$count==1)
temp_entrain2$count <- ifelse(temp_entrain2$Entrain=="entrained", 1,0) 
#fwrite(temp_entrain2, "temperature experiment only entrained_or_freerun.csv")

str(temp_entrain)

#plot temperature overall graph
temp_entrain <- fread("temperature experiment percentage for graphs.csv")


plot <- ggplot(temp_entrain, aes(x=Temperature, y=Percentage*100, group=Actogram))+
  geom_point(size=3)+
  geom_line(aes(linetype=Actogram), size=1)+
  theme_classic()+
  labs(x="temperature increase", y="percentage of entrainment") +
  theme(axis.title = element_text(hjust = 0.5, size = 18, face = "bold"))+
  theme(axis.text=element_text(size=16, face="bold")) +
  theme(legend.text = element_text(size=16), legend.title = element_text(size=16),legend.key.size = unit(2, 'cm'))
plot

#try dose-response curve fitting for temperature entrain exp
temp_entrain2 <- subset(temp_entrain,temp_entrain$Actogram=="Entrain")
str(temp_entrain2)
m1 <- drm(Percentage ~ Temperature, data=temp_entrain2,
          fct = LL.2(names = c("Slope", "ED50")), 
          type="continuous")
summary(m1)
plot(m1, broken=F)

demo.LL.2 <- drm(data = temp_entrain2,Percentage~Temperature,fct=LL.2()) # run model.
summary(demo.LL.2)

# predictions and confidence intervals.
demo.fits <- expand.grid(temp=exp(seq(log(1.00e-01), log(10), length=20))) 
# new data with predictions
pm <- predict(demo.LL.2, newdata=demo.fits, interval="confidence") 
demo.fits$p <- pm[,1]
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]

p <- ggplot(temp_entrain2,aes(x=Temperature,y=Percentage))+
  geom_point()+
  theme_classic()+
  geom_line(data=demo.fits, aes(x=temp, y=p)) +
  labs(x="temperature increase", y="percentage of entrainment")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11)) +
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) #change legend text fo
p


#create a new powerpoint document
doc <- read_pptx("overall data/final results.pptx")
doc <- add_slide(doc, 'Title and Content', 'Office Theme')

#doc <- ph_with(doc, T24 , location = ph_location(width=9,height=3))#for too big plot
# 
fig <- dml(ggobj =pp)
# #add the plot
doc <- ph_with(doc, fig, location = ph_location(width=10,height=10))
# #write the document to a file

print(doc, target ="overall data/final results.pptx")
