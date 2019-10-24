# http://rcompanion.org/handbook/I_09.html
# https://m-clark.github.io/docs/mixedModels/anovamixed.html
library(vegan)
library(lubridate)
library(ggplot2)
library(plyr)
library(lme4)
library(nlme)
library(codyn)
library(tidyverse)
library(mgcv)
library(gratia)


#####################################################################################################
##################################################################################################### DATA
#####################################################################################################

fsr_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2FSR13_ITS2_GenMt0.001.csv")
bk_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2BK13_ITS2_GenMt0.001.csv")
gh_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2GH13_ITS2_GenMt0.001.csv")
cc_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2CC13_ITS2_GenMt0.001.csv")

fsr_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2FSR14_ITS2_GenMt0.001.csv")
bk_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2BK14_ITS2_GenMt0.001.csv")
gh_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2GH14_ITS2_GenMt0.001.csv")
cc_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2CC14_ITS2_GenMt0.001.csv")

#####################################################################################################
##################################################################################################### DIVERSITY
#####################################################################################################

l <- list(gh_14,bk_14,cc_14,fsr_14,gh_13,bk_13,cc_13,fsr_13)
jn <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13')
rcn <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13')
evn <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13')
dvn <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13')
ol <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13')

for (i in 1:length(l)){
  site <- as.vector(rep(substring(ol[[i]],1,2),length(ol[[i]])))
  yr <- as.vector(rep(substring(ol[[i]],4,5),length(ol[[i]])))
  jn[[i]] <- as.vector(yday(as.Date(l[[i]][1:nrow(l[[i]]),1])))
  rcn[[i]] <- as.vector(specnumber(l[[i]][,2:ncol(l[[i]])])) 
  evn[[i]] <- as.vector(diversity(l[[i]][,2:ncol(l[[i]])])/log(specnumber(l[[i]][,2:ncol(l[[i]])]))) 
  dvn[[i]] <- as.vector(diversity(l[[i]][,2:ncol(l[[i]])], index = "shannon"))
  ol[[i]] <- as.data.frame(cbind(jn[[i]],rcn[[i]],evn[[i]],dvn[[i]],site,yr))
  colnames(ol[[i]]) <- c("JD","RC","EV","DV","site","yr")
}

dfDV <- ol[[1]]; for (i in 2:length(ol)){
  dfDV <- rbind(dfDV,ol[[i]])
}

dfDV$yr <- as.numeric(as.character(dfDV$yr))
dfDV$JD <- as.numeric(as.character(dfDV$JD))
dfDV$DV <- (as.numeric(as.character(dfDV$DV)))
dfDV$RC <- (as.numeric(as.character(dfDV$RC)))
dfDV$EV <- (as.numeric(as.character(dfDV$EV)))
dfDV$site = factor(dfDV$site, levels=c('gh','bk','cc','fs'))

ggplot(dfDV, aes(x=JD, y=DV, colour=site)) + geom_point(size=1, aes(shape=site)) + theme_bw() + ylab("Shannon Diversity") + xlab("Julian Date") + 
  geom_smooth(method = "loess",se=F) + facet_grid(rows = vars(yr)) + scale_colour_manual(values = c("black","darkgrey","green3","springgreen4"))
ggplot(dfDV, aes(x=site, y=DV, colour=site)) + theme_bw() + geom_boxplot() + facet_grid(rows=vars(yr)) + ylab("Shannon Diversity") + xlab("Site") + 
  geom_point(size=1) + stat_summary(fun.y=mean,color="black",geom="point",size=2,shape=24) + scale_colour_manual(values = c("black","darkgrey","green3","springgreen4"))


l <- list('cc','fs','gh','bk')
d <- c(0.545,0.073,0.983,0.954)
o <- c(0,0,0,0)
ot <- c(0,0,0,0)
for (i in 1:length(l)){
  o[[i]] <- median(subset(dfDV, dfDV$site==l[[i]] & dfDV$yr==14)$DV)
  ot[[i]] <- median(subset(dfDV, dfDV$site==l[[i]] & dfDV$yr==13)$DV)
  lmdf <- as.data.frame(cbind(d,o,ot))
  colnames(lmdf) <- c("devo","fo","th")
}
lmdf$ave <- (lmdf$th + lmdf$fo)/2

plot(lmdf$devo,lmdf$ave)
cor.test(lmdf$devo,lmdf$ave, alternative = "greater")

plot(lmdf$devo,lmdf$th)
cor.test(lmdf$devo,lmdf$th, alternative = "greater")

plot(lmdf$devo,lmdf$fo)
cor.test(lmdf$devo,lmdf$fo, alternative = "greater")

##################################################################################################### Mean
##################################################################################################### Rank
##################################################################################################### Shift

tidyup <- function(filepath, site) {
  read_csv(filepath) %>%
    gather(date, prop, -genus) %>%
    mutate(site = rep(site, n()),
           date = ymd(date),
           xday = yday(date),
           genus = str_remove(genus, "f__"),
           genus = str_remove(genus, "aceae"),
           genus = str_remove(genus, "g__"),
           year = year(date)) %>%
    unite(sitedate, site, date, remove = FALSE) %>%
    dplyr::select(site, year, date, xday, sitedate, genus, prop)
}
###-------------------------------------------------------------------------------------------
fsr_13 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2FSR13_ITS2_GenMt0.001.csv", "FSR")
bk_13 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2BK13_ITS2_GenMt0.001.csv","BK")
gh_13 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2GH13_ITS2_GenMt0.001.csv","GH")
cc_13 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2CC13_ITS2_GenMt0.001.csv","CC")

fsr_14 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2FSR14_ITS2_GenMt0.001.csv", "FSR")
bk_14 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2BK14_ITS2_GenMt0.001.csv","BK")
gh_14 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2GH14_ITS2_GenMt0.001.csv","GH")
cc_14 <- tidyup("E:/UR2019/RarefiedITS2/tabled/RareITS2CC14_ITS2_GenMt0.001.csv","CC")
###-------------------------------------------------------------------------------------------
dat13 <- full_join(gh_13, bk_13) %>%
  full_join(cc_13) %>%
  full_join(fsr_13) %>%
  spread(genus, prop) %>%
  gather(genus, prop, -site, -date, -sitedate, -year, -xday) %>%
  mutate(genus = factor(genus),
         genus = factor(genus, levels = rev(levels(genus))),
         site = factor(site, levels = c("FSR", "CC", "BK", "GH"))) %>%
  replace_na(list(prop = 0))

dat14 <- full_join(gh_14, bk_14) %>%
  full_join(cc_14) %>%
  full_join(fsr_14) %>%
  spread(genus, prop) %>%
  gather(genus, prop, -site, -date, -sitedate, -year, -xday) %>%
  mutate(genus = factor(genus),
         genus = factor(genus, levels = rev(levels(genus))),
         site = factor(site, levels = c("FSR", "CC", "BK", "GH"))) %>%
  replace_na(list(prop = 0))
###-------------------------------------------------------------------------------------------
TO_13 <- rank_shift(df = dat13,time.var = "xday",species.var = "genus",
                        abundance.var = "prop",replicate.var = "site")
TO_14 <- rank_shift(df = dat14,time.var = "xday",species.var = "genus",
                  abundance.var = "prop",replicate.var = "site")
TO_13$xday <- as.numeric(as.character(substr(TO_13$year_pair,1,3)))
TO_14$xday <- as.numeric(as.character(substr(TO_14$year_pair,1,3)))
TO_13$yr <- "13"
TO_14$yr <- "14"
dfTO <- rbind(TO_13,TO_14)
dfTO$site = factor(dfTO$site, levels=c('GH','BK','CC','FSR'))

ggplot(dfTO, aes(x=xday, y=MRS, colour=site)) + geom_point(size=1, aes(shape=site)) + theme_bw() + ylab("Mean Rank Shift") + xlab("Julian Date") +
  geom_smooth(method = "loess",se=F) + facet_grid(rows = vars(yr), scales="free") + scale_colour_manual(values = c("black","darkgrey","green3","springgreen4"))
ggplot(dfTO, aes(x=site, y=MRS, colour=site)) + theme_bw() + geom_boxplot() + facet_grid(rows=vars(yr), scales="free") + ylab("Mean Rank Shift") + xlab("Site") +
  geom_point(size=1) + stat_summary(fun.y=mean,color="black",geom="point",size=2,shape=24) + scale_colour_manual(values = c("black","darkgrey","green3","springgreen4"))


l <- list('CC','FSR','GH','BK')
d <- c(0.545,0.073,0.983,0.954)
o <- c(0,0,0,0)
ot <- c(0,0,0,0)
for (i in 1:length(l)){
  o[[i]] <- median(subset(dfTO, dfTO$site==l[[i]] & dfTO$yr==14)$MRS)
  ot[[i]] <- median(subset(dfTO, dfTO$site==l[[i]] & dfTO$yr==13)$MRS)
  lmdf <- as.data.frame(cbind(d,o,ot))
  colnames(lmdf) <- c("devo","fo","th")
}
lmdf$ave <- (lmdf$th + lmdf$fo)/2

plot(lmdf$devo,lmdf$ave)
cor.test(lmdf$devo,lmdf$ave, alternative = "greater")

plot(lmdf$devo,lmdf$th)
cor.test(lmdf$devo,lmdf$th, alternative = "greater")

plot(lmdf$devo,lmdf$fo)
cor.test(lmdf$devo,lmdf$fo, alternative = "greater")


##################################################################################################### Qualitative
##################################################################################################### 
##################################################################################################### Turnover

fsr_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2FSR13_ITS2_GenMt0.001.csv")
bk_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2BK13_ITS2_GenMt0.001.csv")
gh_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2GH13_ITS2_GenMt0.001.csv")
cc_13 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2CC13_ITS2_GenMt0.001.csv")

fsr_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2FSR14_ITS2_GenMt0.001.csv")
bk_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2BK14_ITS2_GenMt0.001.csv")
gh_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2GH14_ITS2_GenMt0.001.csv")
cc_14 <- read.csv("E:/UR2019/RarefiedITS2/RareITS2CC14_ITS2_GenMt0.001.csv")

###--------------------------------------------------------------------------

l <- list(gh_14,bk_14,cc_14,fsr_14,gh_13,bk_13,cc_13,fsr_13)
ndl <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # no date
cl <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # ceiling
al <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # apps
dl <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # disses
ol <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # output

for (i in 1:length(l)){
  site <- as.vector(rep(substring(ol[[i]],1,2),(nrow(l[[i]]) - 1)))
  yr <- as.vector(rep(substring(ol[[i]],4,5),(nrow(l[[i]]) - 1)))
  jn <- as.vector(yday(as.Date(l[[i]][2:nrow(l[[i]]),1])))
  ndl[[i]] <- l[[i]][,2:ncol(l[[i]])]
  cl[[i]] <- ceiling(ndl[[i]])
  j <- 1:(nrow(cl[[i]])-1)
  j1 <- 1 + j
  al[[i]] <- as.data.frame(cl[[i]][j1,] > cl[[i]][j,])*1
  dl[[i]] <- as.data.frame(cl[[i]][j1,] < cl[[i]][j,])*1
  al[[i]]$app <- rowSums(al[[i]]) 
  dl[[i]]$dis <- rowSums(dl[[i]]) 
  ol[[i]] <- as.data.frame(cbind(site, yr, jn, dl[[i]]$dis,al[[i]]$app))
  colnames(ol[[i]]) <- c("site","yr","day","dis","app")
}

dfQT <- ol[[1]]; for (i in 2:length(ol)){
  dfQT <- rbind(dfQT,ol[[i]])
}

dfQT$yr <- as.numeric(as.character(dfQT$yr))
dfQT$day <- as.numeric(as.character(dfQT$day))
dfQT$qTurn <- (as.numeric(as.character(dfQT$dis)) + as.numeric(as.character(dfQT$app)))
dfQT$site = factor(dfQT$site, levels=c('gh','bk','cc','fs'))

ggplot(dfQT, aes(x=day, y=qTurn, colour=site)) + geom_point(size=1, aes(shape=site)) + theme_bw() + ylab("Qualitative Turnover") + xlab("Julian Date") +
  geom_smooth(method = "loess",se=F) + facet_grid(rows = vars(yr), scales="free") + scale_colour_manual(values = c("black","darkgrey","green3","springgreen4"))
ggplot(dfQT, aes(x=site, y=qTurn, colour=site)) + theme_bw() + geom_boxplot() + facet_grid(rows=vars(yr), scales="free") + ylab("Qualitative Turnover") + xlab("Site") +
  geom_point(size=1) + stat_summary(fun.y=mean,color="black",geom="point",size=2, shape=24) + scale_colour_manual(values = c("black","darkgrey","green3","springgreen4"))


l <- list('cc','fs','gh','bk')
d <- c(0.545,0.073,0.983,0.954)
o <- c(0,0,0,0)
ot <- c(0,0,0,0)
for (i in 1:length(l)){
  o[[i]] <- median(subset(dfQT, dfQT$site==l[[i]] & dfQT$yr==14)$qTurn)
  ot[[i]] <- median(subset(dfQT, dfQT$site==l[[i]] & dfQT$yr==13)$qTurn)
  lmdf <- as.data.frame(cbind(d,o,ot))
  colnames(lmdf) <- c("devo","fo","th")
}
lmdf$ave <- (lmdf$th + lmdf$fo)/2

plot(lmdf$devo,lmdf$ave)
cor.test(lmdf$devo,lmdf$ave, alternative = "greater")

plot(lmdf$devo,lmdf$th)
cor.test(lmdf$devo,lmdf$th, alternative = "greater")

plot(lmdf$devo,lmdf$fo)
cor.test(lmdf$devo,lmdf$fo, alternative = "greater")













































##################################################################################################### 
##################################################################################################### GAMs
##################################################################################################### 

dfDV13 <- subset(dfDV, dfDV$yr==13)
dfDV14 <- subset(dfDV, dfDV$yr==14)
summary(gam(DV ~ te(JD, by=site, bs="fs") + site + as.factor(yr), data = dfDV), select=T, method="REML") # try te() allows different wiggliness between sites
summary(gam(DV ~ s(JD, by=site) + as.factor(yr), data = dfDV), select=T, method="REML") # try te() allows different wiggliness between sites
summary(gam(DV ~ s(JD, by=site, bs="re") + as.factor(yr), data = dfDV), select=T, method="REML")
g1 <- gam(DV ~ s(JD, site, bs="fs") + as.factor(yr), data=dfDV, method="REML")
summary(g1) # http://astrostatistics.psu.edu/su07/R/html/mgcv/html/summary.gam.html
draw(g1)

summary(gam(DV ~ s(JD, site, bs="fs"), data=dfDV13, method="REML"))
summary(gam(DV ~ s(JD, site, bs="fs"), data=dfDV14, method="REML"))
summary(gam(DV ~ s(JD) + site, data=dfDV13, method="REML"))
summary(gam(DV ~ s(JD) + site, data=dfDV14, method="REML"))

summary(gam(DV ~ s(JD, bs="fs") + site + as.factor(yr), data = dfDV), select=T, method="REML")


plot(gam(DV ~ s(JD, site, bs="fs") + as.factor(yr), data = dfDV), select=T, method="REML")
summary(gam(DV ~ s(JD) + site, data = dfDV14))
plot.gam(gam(DV ~ s(JD) + site + yr, data = dfDV), all.terms = T)
gam.check(gam(DV ~ s(JD) + site + yr, data = dfDV))

x <- gam(DV ~ te(JD, by=site, bs="gp") + site + as.factor(yr), data = dfDV, select=T, method="REML") # try te() allows different wiggliness between sites
plot.gam(x, pages = 1, all.terms = T)

# https://www.fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
# https://petolau.github.io/Analyzing-double-seasonal-time-series-with-GAM-in-R/
t

?plot.gam




##################################################################################################### 
##################################################################################################### dH/dT
##################################################################################################### 

fsr_13 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/FSR13_Med_GenPf0.0005Mt0.001.csv")
bk_13 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/BK13_Med_GenPf0.0005Mt0.001.csv")
gh_13 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/GH13_Med_GenPf0.0005Mt0.001.csv")
cc_13 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/CC13_Med_GenPf0.0005Mt0.001.csv")

fsr_14 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/FSR14_Med_GenPf0.0005Mt0.001.csv")
bk_14 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/BK14_Med_GenPf0.0005Mt0.001.csv")
gh_14 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/GH14_Med_GenPf0.0005Mt0.001.csv")
cc_14 <- read.csv("E:/UR2019/EcoMetrics/BC_dissimilarity/CC14_Med_GenPf0.0005Mt0.001.csv")

###--------------------------------------------------------------------------

l <- list(gh_14,bk_14,cc_14,fsr_14,gh_13,bk_13,cc_13,fsr_13)
Pi1 <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # Pi 1
Pi2 <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # Pi 2
dPi <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # delta Pi
GdHdT <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # contribution of each genus
dHdT <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # sum contribution
ol <- list('gh_14','bk_14','cc_14','fs_14','gh_13','bk_13','cc_13','fs_13') # output

for (i in 1:length(l)){
  l[[i]][l[[i]] == 0] <- NA
  site <- as.vector(rep(substring(ol[[i]],1,2),(nrow(l[[i]]) - 1)))
  yr <- as.vector(rep(substring(ol[[i]],4,5),(nrow(l[[i]]) - 1)))
  Pi1[[i]] <- t(l[[i]][1:(nrow(l[[i]])-1),2:ncol(l[[i]])])
  Pi2[[i]] <- t(l[[i]][2:(nrow(l[[i]])),2:ncol(l[[i]])])
  dPi[[i]] <- (Pi1[[i]] - Pi2[[i]])
  d1 <- yday(l[[i]][1:(nrow(l[[i]])-1),1])
  d2 <- yday(l[[i]][2:(nrow(l[[i]])),1])
  dDay <- (d2 - d1)
  GdHdT[[i]] <- (dPi[[i]]/dDay[[i]])*(1+log(Pi1[[i]]))
  dHdT <- (-1)*colSums(GdHdT[[i]], na.rm = T)
  ol[[i]] <- as.data.frame(cbind(site, yr, d1, dHdT))
  colnames(ol[[i]]) <- c("site","yr","day","dHdT")
}

dfQT <- ol[[1]]; for (i in 2:length(ol)){
  dfQT <- rbind(dfQT,ol[[i]])
}

dfQT$yr <- as.numeric(as.character(dfQT$yr))
dfQT$day <- as.numeric(as.character(dfQT$day))
dfQT$dHdT <- as.numeric(as.character(dfQT$dHdT))

ggplot(df, aes(day, dHdT, color=site)) + 
  geom_line(size= 2) + theme_bw() + facet_grid(rows = vars(yr))
ggplot(df, aes(day, dHdT, colour=site)) + 
  geom_point(size=2) + theme_bw() +
  geom_smooth(se = F) + facet_grid(rows = vars(yr))

