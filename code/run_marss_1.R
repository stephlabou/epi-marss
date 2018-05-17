#test run w/ mar formatted data

#load packages
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(MAR1)
library(MARSS)

#read in formatted data
#(NAs replaced and zeros replaced)
dat <- read.csv("../data/fixed_mar_dat.csv", stringsAsFactors = FALSE)

summary(dat)
head(dat)

#reorg columns and format column names
dat_format <- dat %>% 
              select(year, month, fixed_chla_avg, fixed_chla_sum, fixed_temp_avg,
                     fixed_Cyclops_adult_count, fixed_Epischura_adult_count,
                     fixed_Epischura_copep_count, fixed_Epischura_naup_count,
                     fixed_Achnanthes_density, fixed_Aulacoseira_density,
                     fixed_Chroomonas_density, fixed_Chrysidalis_density,
                     fixed_Dinobryon_density, fixed_Nitzchia_density, fixed_picoplankton_density,
                     fixed_Romeria_density, fixed_Stephanodiscus_density,
                     fixed_Synedra_density, fixed_unid_nano_density, fixed_unid_unid_density) %>% 
              arrange(year, month)

names(dat_format) <- gsub("fixed_", "", names(dat_format))
names(dat_format) <- gsub("_density", "", names(dat_format))
names(dat_format) <- gsub("_count", "", names(dat_format))

#convert phyto from thousands of cell/liter to cells/liter
dat_format <- dat_format %>% 
              mutate(Achnanthes = Achnanthes*1000,
                     Aulacoseira = Aulacoseira*100,
                     Chroomonas = Chroomonas*1000,
                     Chrysidalis = Chrysidalis*1000,
                     Dinobryon = Dinobryon*1000,
                     Nitzchia = Nitzchia*1000,
                     picoplankton = picoplankton*1000,
                     Romeria = Romeria*1000,
                     Stephanodiscus = Stephanodiscus*1000,
                     Synedra = Synedra*1000,
                     unid_nano = unid_nano*1000,
                     unid_unid = unid_unid*1000)

#write this to csv for SH
write.csv(dat_format, "../data/Baikal_dat_mar_formatted_20171117.csv", row.names = FALSE)


#ok, now in correct format for mar
#remove time, keep total chla
dat_format <- mutate(dat_format, time_index = paste(year, month, sep = "_")) %>% select(-year, -month, -chla_avg)
head(dat_format)

#remove time index, but it's here...

#format main data
dat.mar <- dat_format %>% select(-time_index)
dat.mar <- as.matrix(dat.mar)
dat.mar.t <- t(dat.mar)
str(dat.mar.t)
head(dat.mar.t)
dat.mar.t[1:10, 1:10]

the.mean = apply(dat.mar.t,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(dat.mar.t,1,var,na.rm=TRUE))
dat.mar.z = (dat.mar.t-the.mean)*(1/the.sigma)
dat.mar.z[1:10, 1:10]
#looks like z-score applied

#covars = chla, temp, phytos



