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
##################################################################################################### Diversity
#####################################################################################################

fsr_13 <- read.csv("RareITS2FSR13_ITS2_GenMt0.001.csv")
bk_13 <- read.csv("RareITS2BK13_ITS2_GenMt0.001.csv")
gh_13 <- read.csv("RareITS2GH13_ITS2_GenMt0.001.csv")
cc_13 <- read.csv("RareITS2CC13_ITS2_GenMt0.001.csv")

fsr_14 <- read.csv("RareITS2FSR14_ITS2_GenMt0.001.csv")
bk_14 <- read.csv("RareITS2BK14_ITS2_GenMt0.001.csv")
gh_14 <- read.csv("RareITS2GH14_ITS2_GenMt0.001.csv")
cc_14 <- read.csv("RareITS2CC14_ITS2_GenMt0.001.csv")

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
fsr_13 <- tidyup("RareITS2FSR13_ITS2_GenMt0.001.csv", "FSR")
bk_13 <- tidyup("RareITS2BK13_ITS2_GenMt0.001.csv","BK")
gh_13 <- tidyup("RareITS2GH13_ITS2_GenMt0.001.csv","GH")
cc_13 <- tidyup("RareITS2CC13_ITS2_GenMt0.001.csv","CC")

fsr_14 <- tidyup("RareITS2FSR14_ITS2_GenMt0.001.csv", "FSR")
bk_14 <- tidyup("RareITS2BK14_ITS2_GenMt0.001.csv","BK")
gh_14 <- tidyup("RareITS2GH14_ITS2_GenMt0.001.csv","GH")
cc_14 <- tidyup("RareITS2CC14_ITS2_GenMt0.001.csv","CC")
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

fsr_13 <- read.csv("RareITS2FSR13_ITS2_GenMt0.001.csv")
bk_13 <- read.csv("RareITS2BK13_ITS2_GenMt0.001.csv")
gh_13 <- read.csv("RareITS2GH13_ITS2_GenMt0.001.csv")
cc_13 <- read.csv("RareITS2CC13_ITS2_GenMt0.001.csv")

fsr_14 <- read.csv("RareITS2FSR14_ITS2_GenMt0.001.csv")
bk_14 <- read.csv("RareITS2BK14_ITS2_GenMt0.001.csv")
gh_14 <- read.csv("RareITS2GH14_ITS2_GenMt0.001.csv")
cc_14 <- read.csv("RareITS2CC14_ITS2_GenMt0.001.csv")

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


