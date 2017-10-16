###########################################
####  Prepare data for MARSS analysis  ####
###########################################

## TO DO

# ----> redone to use monthly averages

# [ ] are layer groups range inclusive? (0-10 group with <10 vs <=10)

# [ ] what do we want to do about legitimate NA values?
#   ----> see agg_depth_notes.doc

# ----> these turn out to be moot, since outside 50m range
# [ ] layer mistakes, so count is negative or infinite
#     [ ] what to do about NA layers 
#     [ ] counts that are NA/NaN

#############################################
####  Packages, working directory, etc.  ####
#############################################

## We are using data since 1975 and above 50m (inclusive)

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

## Load packages
library(dplyr)
library(lubridate)
library(reshape2)
library(MAR1)
library(MARSS)

## Directory of the long-term data Subversion repository. UPDATE THIS to run the
## code on your own machine.
dir <- "D:/Labou/Baikal/baikal/Longterm_data/"

########################################################
##################  Load zoop data  ####################
########################################################

## Note: this is the zooplankton data with zeroes already in place when species codes weren't observed
fulldat <- read.csv(paste0(dir, "zoo/data/zoopzeroskey_alldepths.csv"),
                    stringsAsFactors = FALSE)
colnames(fulldat) <- tolower(colnames(fulldat))

## Load zooplankton key and make column names lower case
key <- read.csv(paste0(dir, "zoo/data/key.csv"), stringsAsFactors = FALSE)
colnames(key) <- tolower(colnames(key))

######################################################################
#########  Keep only Epischura baicalensis & adult Cyclops  ##########
######################################################################

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

###############################
####  Convert count units  ####
###############################

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

#################################
####  View lifestage groups  ####
#################################

# zoop_fix_corr_units %>%
#   select(genus, lifestage_cop) %>%
#   unique() %>%
#   arrange(genus, lifestage_cop)
# zoop_fix_corr_units %>%
#   group_by(genus, lifestage_cop) %>%
#   summarize(n = n_distinct(date))
  
#keep only adult cyclops

zoop_ready <- zoop_fix_corr_units %>% 
              filter(ifelse(genus == "Cyclops", lifestage_cop == "adult", lifestage_cop %in% c("adult", "copep", "naup")))   

# zoop_ready %>%
#   select(genus, lifestage_cop) %>%
#   unique() %>%
#   arrange(genus, lifestage_cop)
# zoop_ready %>%
#   group_by(genus, lifestage_cop) %>%
#   summarize(n = n_distinct(date))

#######################################################
####  Average temperature data by month and depth  ####
#######################################################

temp <- read.csv(paste0(dir, "temp_chl_secchi_wind/cleaned_data/temp_cleaned.csv"),
                 stringsAsFactors = FALSE)

# temp_month <- temp %>%
#   ## Create column of months
#   mutate(monthyear = floor_date(as.Date(date), "month")) %>%
#   ## Average temp by month
#   group_by(monthyear) %>%
#   summarize(temp = mean(temp, na.rm = TRUE))

######################################################
####  Average chlorophyll data by month and depth ####
######################################################

chla <- read.csv(paste0(dir, "temp_chl_secchi_wind/cleaned_data/chla_cleaned.csv"),
                 stringsAsFactors = FALSE)

# chla_month <- chla %>%
#   ## Create column of months
#   mutate(monthyear = floor_date(as.Date(date), "month")) %>%
#   ## Average chla by month
#   group_by(monthyear) %>%
#   summarize(chla = mean(chla, na.rm = TRUE))

##############################################
######  Read in and organize phyto data ######
##############################################

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

phy_newgen <- read.csv("../data/phy_with_three_unidentified_groups.csv",
                       stringsAsFactors = FALSE)

##############################################
##### Organize groups (zoop, phyto, env) #####
##############################################

## ----> phyto groups

phyto_groups <- c("Cyanodictyon", "Synechocystis", "unid_pico",
                  "Romeria", "unid_nano", "Chrysidalis", "Nitzchia",
                  "Dinobryon", "Chroomonas", "Stephanodiscus",
                  "Synedra", "Achnanthes", "Aulacoseira", "unid_unid")


phy_dat_grouped <- phy_newgen %>% 
      select(-genus) %>% 
      #keep only genera of interest
      filter(genus_revised %in% phyto_groups) %>% 
      #make new merged group
      mutate(genus_groups = ifelse(genus_revised %in% c("Cyanodictyon", "Synechocystis", "unid_pico"),
                                   "picoplankton", genus_revised)) %>% 
      #already limited to 1975; limit to <= 50 meters
      filter(depth <= 50) %>% 
      #aggregate within date, depth, and genus group
      group_by(year, month, date, depth, genus_groups) %>% 
      summarize(density_genus_sum = sum(density)) %>% 
      as.data.frame()

# filter(phy_newgen, date == "1975-01-07" & genus_revised == "Nitzchia" & depth == 20)
# filter(phy_dat_grouped, date == "1975-01-07" & genus_groups == "Nitzchia" & depth == 20)
# 
# filter(phy_newgen, date == "1999-11-04" & genus_revised == "Aulacoseira" & depth == 50)
# filter(phy_dat_grouped, date == "1999-11-04" & genus_groups == "Aulacoseira" & depth == 50)


## ----> zoop groups

zoop_grouped <- zoop_ready %>% 
        mutate(layer_group = paste(upper_layer, lower_layer, sep = "-")) %>% 
        #for now, sticking with sequential depths to 50 m
        filter(layer_group %in% c("0-10", "10-25", "25-50")) %>% 
        mutate(layer_group = ordered(layer_group, levels = c("0-10", "10-25", "25-50"))) %>% 
        #group and aggregate counts within Epi lifestages (adult, copep, naup)
        group_by(date, layer_group, upper_layer, lower_layer, genus, lifestage_cop) %>% 
        summarize(count_l_sum = sum(count_l)) %>% 
        #keep only since 1975
        mutate(year = year(as.Date(date))) %>% 
        filter(year >= 1975) %>% 
        as.data.frame()

# filter(zoop_fix_corr_units, date == "1975-01-07" & lifestage_cop == "naup" & upper_layer == 0 & lower_layer == 10)
# filter(zoop_grouped, date == "1975-01-07" & lifestage_cop == "naup" & upper_layer == 0 & lower_layer == 10)
# 
# filter(zoop_fix_corr_units, date == "2003-12-06" & lifestage_cop == "copep" & upper_layer == 10 & lower_layer == 25)
# filter(zoop_grouped, date == "2003-12-06" & lifestage_cop == "copep" & upper_layer == 10 & lower_layer == 25) 

# filter(zoop_fix_corr_units, date == "1975-03-10" & lifestage_cop == "adult" & upper_layer == 0 & lower_layer == 10)
# filter(zoop_grouped, date == "1975-03-10" & lifestage_cop == "adult" & upper_layer == 0 & lower_layer == 10)
# 
# filter(zoop_fix_corr_units, date == "1994-11-09" & lifestage_cop == "adult" & upper_layer == 25 & lower_layer == 50)
# filter(zoop_grouped, date == "1994-11-09" & lifestage_cop == "adult" & upper_layer == 25 & lower_layer == 50)


# ----> limit chla and temp to 1975 and 50m

chla_match <- chla %>% 
        #looks like chla already all since 1975...
        #but just in case...
        filter(year(as.Date(date))>=1975) %>% 
        filter(depth <= 50)

temp_match <- temp %>% 
        filter(year(as.Date(date))>=1975) %>% 
        filter(depth <= 50)

# ----> check out separate existing data frames

head(chla_match)
head(temp_match)
head(phy_dat_grouped)
head(zoop_grouped)

# ----> merge chla/temp - match date and depth

#merge and make long - NAs ok - can interpolate later

env_merge <- chla_match %>% 
        select(-month, -Year) %>% 
        merge(temp_match, by = c("date", "depth"), all = TRUE)

# ----> reshape phyto so genera along top

phy_dat_grouped_wide <- phy_dat_grouped %>% 
      select(-year, -month) %>% 
      #would feel more comfortable moving forward if "density sum" or other was in name
      #so we can keep track density vs count
      #may need to change sep later (since _ in unid_ names...)
      mutate(genus_groups = paste(genus_groups, "density", sep = "_")) %>% 
      dcast(date + depth ~ genus_groups, value.var = "density_genus_sum") 

# ----> merge env with phyto and make approx layer groups

phyto_env_wide_merge <- merge(env_merge, phy_dat_grouped_wide,
                              by = c("date", "depth"))

phyto_env_ready <- phyto_env_wide_merge %>% 
                    mutate(approx_layer = ifelse(depth <= 10, "0-10", NA),
                           approx_layer = ifelse(depth >10 & depth <= 25, "10-25", approx_layer),
                           approx_layer = ifelse(depth > 25 & depth <= 50, "25-50", approx_layer)) %>% 
                    rename(phyto_env_depth = depth)

phyto_env_ready <- phyto_env_ready[,c(1, 2, 17, 3:16)]


# ----> reshape zoop so genera along top

zoop_ready <- zoop_grouped %>% 
                    mutate(genus_stage = paste(genus, lifestage_cop, sep = "_")) %>% 
                    select(-genus, -lifestage_cop, - year) %>% 
                    #put count/sum in name so distinguish from phyto density cols
                    #probably need to come back and change sep so unique genus/stage/type sep
                    mutate(genus_stage = paste(genus_stage, "count", sep = "_")) %>% 
                    dcast(date + layer_group + upper_layer + lower_layer ~ genus_stage, value.var = "count_l_sum") %>% 
                    rename(zoop_upper_layer = upper_layer, zoop_lower_layer = lower_layer)

##############################################
###### Merge for MARSS dataset draft 1 #######
##############################################

head(phyto_env_ready)
head(zoop_ready)

# ----> full merge

mardat <- merge(phyto_env_ready, zoop_ready,
                by.x = c("date", "approx_layer"),
                by.y = c("date", "layer_group"),
                all = TRUE) %>% 
          rename(layer = approx_layer) %>% 
          #reorder columns
          select(date, layer, zoop_upper_layer, zoop_lower_layer, phyto_env_depth,
                 chla, temp, Achnanthes_density, Aulacoseira_density, Chroomonas_density,
                 Chrysidalis_density, Dinobryon_density, Nitzchia_density,
                 picoplankton_density, Romeria_density, Stephanodiscus_density,
                 Synedra_density, unid_nano_density, unid_unid_density,
                 Cyclops_adult_count, Epischura_adult_count,
                 Epischura_copep_count, Epischura_naup_count) %>% 
          #make date a date
          mutate(date = as.Date(date)) %>% 
          #make layer grouped ordered factor
          mutate(layer = ordered(layer, levels = c("0-10", "10-25", "25-50")))
                 
head(mardat)
str(mardat)
summary(mardat)

write.csv(mardat, "../data/mar_dat.csv", row.names = FALSE)


# example MAR data
# data(L4.AllDates)
# head(L4.AllDates)
