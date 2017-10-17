# ----> redone to use monthly averages
# limited to since 1975, 0-50 meters

# [ ] what do we want to do about legitimate NA values?
#   ----> see agg_depth_notes.doc

###########################################
####  Prepare data for MARSS analysis  ####
###########################################

## We are using data since 1975 and above 50m (inclusive)

## Temperature, chla, phyto, and zoop (Epischura/adult cyclops)

#From SH, groups are:

# 1.	Cyanodictyon + Synechocystic + unid_pico
# 2.	Romeria
# 3.	unid_nano
# 4.	Chrysidalis
# 5.	Nitzchia
# 6.	Dinobryon
# 7.	Chroomonas
# 8.	Stephanodiscus
# 9.	Synedra
# 10.	Achnanthes
# 11.	Melosira/Aulacoseira (in this dataset as "Aulacoseira")
# 12.	unid_unid
# 13.	Epischura adult
# 14.	Epischura copep
# 15.	Epischura nauplii
# 16.	Cyclops adult

#(1-12 are phyto, 13-16 are zoop)

#############################################
####  Packages, working directory, etc.  ####
#############################################

## Load packages
library(dplyr)
library(lubridate)
library(reshape2)
library(MAR1)
library(MARSS)

## Directory of the long-term data Subversion repository. UPDATE THIS to run the
## code on your own machine.
dir <- "D:/Labou/Baikal/baikal/Longterm_data/"

############################################
####  Average temperature data by month ####
############################################

temp <- read.csv(paste0(dir, "temp_chl_secchi_wind/cleaned_data/temp_cleaned.csv"),
                 stringsAsFactors = FALSE)

# Want monthly average across depths 0-50m
temp_small <- temp %>% 
            mutate(year = year(as.Date(date)),
                   month = month(as.Date(date))) %>% 
            #keep since 1975, <= 50m
            filter(year >= 1975 & depth <= 50)
            

#aggregate across depths by month and year
temp_monthly <- temp_small %>% 
                group_by(year, month) %>% 
                summarize(temp_050depth_avg = mean(temp, na.rm = TRUE)) %>% 
                as.data.frame()


######################################################
####  Aggregate/average chlorophyll data by month ####
######################################################

chla <- read.csv(paste0(dir, "temp_chl_secchi_wind/cleaned_data/chla_cleaned.csv"),
                 stringsAsFactors = FALSE)

# Want monthly average across depths 0-50m
chla_small <- chla %>% 
            #keep only since 1975, 50m and above
            filter(Year >= 1975 & depth <= 50) %>% 
            rename(year = Year)

#aggregate across depths by month and year
chla_monthly <- chla_small %>% 
                group_by(year, month) %>% 
                summarize(chla_050depth_avg = mean(chla, na.rm = TRUE),
                          chla_050depth_sum = sum(chla, na.rm = TRUE)) %>% 
                as.data.frame()


##################################################
#### Wrangle zoop data and aggregate by month ####
##################################################

## Note: this is the zooplankton data with zeroes already in place when species codes weren't observed
fulldat <- read.csv(paste0(dir, "zoo/data/zoopzeroskey_alldepths.csv"),
                    stringsAsFactors = FALSE)
colnames(fulldat) <- tolower(colnames(fulldat))

## Load zooplankton key and make column names lower case
key <- read.csv(paste0(dir, "zoo/data/key.csv"), stringsAsFactors = FALSE)
colnames(key) <- tolower(colnames(key))

## ----> Keep only Epischura baicalensis & adult Cyclops 

## Codes that represent Epischura/Cyclops non-double counts according to Derek's key:
epi_key <- filter(key, genus == "Epischura" & species == "baicalensis"& doublecount %in% c("N", "M")) %>% .$kod
cyclp_key <- filter(key, genus == "Cyclops" & doublecount %in% c("N", "M")) %>% .$kod

keep_key <- c(epi_key, cyclp_key) %>% unique() #unique just in case weirdness

## Keep only the Epischura/Cyclops non-double counts
zoop <- filter(fulldat, kod %in% keep_key) %>%
  ## Rename some columns
  rename(code = kod, date = date, upper_layer = ver_gr, lower_layer = nig_gr)

#zoop %>% select(genus, species) %>% unique()
#there are NA values which are useless...so refiltering to keep only Epi and Cyclops
zoop_fix <- filter(zoop, genus %in% c("Epischura", "Cyclops"))

## ----> Convert count units  

## Count is reported as individuals * 1000 / m2. To convert to
## individuals/liter, use the following function:

m2_to_l <- function(x, interval) {
  stopifnot(is.numeric(x))                           # count should be numeric
  ## takes units in 1000 individuals/m2 and converts to individuals per liter
  individuals <- x * 1000                            # convert to individuals/m2
  count_per_liter <- individuals / (interval * 1000) # convert to indiv./liter
  return(count_per_liter)
}

## Convert units to individuals / liter - NOTE this is going to need to be fixed
## because there are  some mistakes in the upper/lower intervals  on a couple of
## dates. This  causes negative  counts (if  upper is  below lower)  or infinite
## counts (if upper and lower are the same). 
zoop_fix_corr_units <- zoop_fix %>%
  mutate(count_l = m2_to_l(count, interval = lower_layer - upper_layer)) %>%
  #select(-count) %>%
  ## Filter out data that's listed as being from 1908 -- this is clearly a mistake
  filter(as.Date(date) >= as.Date("1945-01-01"))

## For now, can ignore the layer mismatch issue, since they call outside our depth range of interest
## Other issues: instances where date or layer range is missing (real NA in orig data)
## Going to exclude these, since without date or depth, data is not useful to us


#keep only adult cyclops

zoop_ready <- zoop_fix_corr_units %>% 
              filter(ifelse(genus == "Cyclops", lifestage_cop == "adult", lifestage_cop %in% c("adult", "copep", "naup")))   


## ----> keep only since 1975, <= 50m, aggregate to monthly 

#using the sequential sampling depths
#then aggregates ACROSS the 3 sequential layers

zoop_monthly <- zoop_ready %>% 
              mutate(layer_group = paste(upper_layer, lower_layer, sep = "-")) %>% 
              #keep sequential sampling depths to 50 m
              filter(layer_group %in% c("0-10", "10-25", "25-50")) %>% 
              mutate(layer_group = ordered(layer_group, levels = c("0-10", "10-25", "25-50"))) %>% 
              #keep only since 1975
              mutate(year = year(as.Date(date)),
                     month = month(as.Date(date))) %>% 
              filter(year >= 1975) %>% 
              group_by(year, month, genus, lifestage_cop) %>% 
              summarize(count_l_sum = sum(count_l, na.rm = TRUE)) %>% 
              as.data.frame()


###################################################
#### Wrangle phyto data and aggregate by month ####
###################################################

#### RETAINED FROM ORIG ####
# phy <- read.csv(paste0(dir, "phyto/Baikal_Phyto_zeroesInc_2countRemoved_KeyUpdate_20120728.csv"), stringsAsFactors = FALSE)
# names(phy) <- tolower(names(phy))
# phy$date <- as.Date(phy$date, "%m/%d/%Y")
# 
# ## Only keep data from >=1975 due to change in algae preservation technique
# phy <- filter(phy, year >= 1975)
# 
# ## Separate "unidentified" taxa into pico, nano, and unknown groups -- first
# ## figure out which codes from Lyubov's key correspond to which group
# 
# unid_group <- function(x) {
#   ## Function to group unidentified genera into picoplankton, nanoplankton, and
#   ## unknown using size information from Lyubov's plankton key
#   if (x %in% c(526, 9526, 586, 587)) {  # pico
#     result <- "unid_pico"
#   } else if (x %in% c(527, 529, 583, 483, 550, 8550, 9550, 551, 8551, 552, 8552,
#                       553, 8553, 554, 584, 580, 581, 582, 535, 536, 537, 538,
#                       9538, 539, 9539, 540, 542, 543, 557, 570, 571, 572, 573,
#                       574, 532, 531, 563, 564, 560, 561, 562)) { # nano
#     result <- "unid_nano"
#   } else if (x %in% c(1, 119, 120, 141, 193, 206, 211, 236, 237, 250, 263, 414,
#                       443, 484, 521, 530, 533, 534, 549, 565, 566, 569, 575,
#                       579, 588, 589, 590, 593, 594, 595, 596, 597, 598, 599,
#                       8120, 8549, 8566, 9533, 9534, 9551)) { # unknown
#     result <- "unid_unid"
#   } else {
#     result <- unique(phy[phy$code == x, "genus"])
#   }
#   return(result)
# }
## Apply the funciton above to create a revised genus column that separates
## "unidentified" into "unid_pico", "unid_nano", and "unid_unid". This takes a
## long time, it needs to run overnight.
# phy_newgen <- phy %>%
#   rowwise() %>%
#   mutate(genus_revised = unid_group(code))
#
# ## Export to CSV so I don't have to rerun this every time.
# write.csv(phy_newgen,
#           "../data/phy_with_three_unidentified_groups.csv",
#           row.names = FALSE)

## ----> read in phyto data with unidentified groups

phy_newgen <- read.csv("../data/phy_with_three_unidentified_groups.csv",
                       stringsAsFactors = FALSE)

## ----> phyto groups

phyto_groups <- c("Cyanodictyon", "Synechocystis", "unid_pico",
                  "Romeria", "unid_nano", "Chrysidalis", "Nitzchia",
                  "Dinobryon", "Chroomonas", "Stephanodiscus",
                  "Synedra", "Achnanthes", "Aulacoseira", "unid_unid")


phy_monthly <- phy_newgen %>% 
                  select(-genus) %>% 
                  #keep only genera of interest
                  filter(genus_revised %in% phyto_groups) %>% 
                  #make new merged group
                  mutate(genus_groups = ifelse(genus_revised %in% c("Cyanodictyon", "Synechocystis", "unid_pico"),
                                               "picoplankton", genus_revised)) %>% 
                  #already limited to 1975; limit to <= 50 meters
                  filter(depth <= 50) %>% 
                  #aggregate within date, depth, and genus group to get total density
                  group_by(year, month, genus_groups) %>% 
                  summarize(density_genus_sum = sum(density)) %>% 
                  as.data.frame()


######################################################
#### Combine monthly data frames for MARSS format ####
######################################################

# ----> reshape phyto and zoop so genera along top

phy_monthly_wide <- phy_monthly %>% 
                        #would feel more comfortable moving forward if "density sum" or other was in name
                        #so we can keep track density vs count
                        #may need to change sep later (since _ in unid_ names...)
                        mutate(genus_groups = paste(genus_groups, "density", sep = "_")) %>% 
                        dcast(year + month ~ genus_groups, value.var = "density_genus_sum") 

zoop_monthly_wide <- zoop_monthly %>% 
                      mutate(genus_stage = paste(genus, lifestage_cop, sep = "_")) %>% 
                      select(-genus, -lifestage_cop) %>% 
                      #put count/sum in name so distinguish from phyto density cols
                      #probably need to come back and change sep so unique genus/stage/type sep
                      mutate(genus_stage = paste(genus_stage, "count", sep = "_")) %>% 
                      dcast(year + month ~ genus_stage, value.var = "count_l_sum") 

# ----> merge chla/temp/phyto/zoop

head(chla_monthly)
head(temp_monthly)
head(zoop_monthly_wide)
head(phy_monthly_wide)

dat_full <- temp_monthly %>% 
            merge(chla_monthly, by = c("year", "month"), all = TRUE) %>% 
            merge(zoop_monthly_wide, by = c("year", "month"), all = TRUE) %>% 
            merge(phy_monthly_wide, by = c("year", "month"), all = TRUE)

str(dat_full)
summary(dat_full)


# ----> mwrite to csv for later use
write.csv(dat_full, "../data/mar_dat.csv", row.names = FALSE)

