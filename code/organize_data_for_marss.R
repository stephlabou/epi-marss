## TO DO

# [ ] are layer groups range inclusive? (0-10 group with <10 vs <=10)


###########################################
####  Prepare data for MARSS analysis  ####
###########################################

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
####################  Load data  #######################
########################################################

## Load the main zooplankton data file and make column names lower case 
## Note: this is the data with zeroes already in place when species codes weren't observed
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

# a <- keep_key %>% sort()
# b <- zoop$code %>% unique() %>% sort()
# #no 32, 88...
# 
# missing <- a[!(a %in% b)]
# b[!(b %in% a)]
# 
# filter(fulldat, kod %in% missing) #ok, so it's just not there
# filter(key, kod %in% missing) %>% select(notes) %>% unique() #indeed, i see that NOW
# 
# filter(fulldat, kod %in% keep_key) %>% select(genus) %>% unique()
# filter(fulldat, kod %in% keep_key) %>% filter(is.na(genus)) %>% head()
# filter(fulldat, kod == 100) %>% head()
# filter(key, kod == 100) %>% head()
# 
# #arggg there's lots of NA rows in data...
# 
# filter(fulldat, kod %in% epi_key) %>% select(genus) %>% unique()
# filter(fulldat, kod %in% cyclp_key) %>% select(genus) %>% unique() #ah HA the na is coming from the cyclops key
# filter(fulldat, kod %in% cyclp_key) %>% filter(is.na(genus)) %>% head()
# weird_key <- filter(fulldat, kod %in% cyclp_key) %>% filter(is.na(genus)) %>% select(kod) %>% unique()
# filter(key, kod %in% weird_key) #??????
# filter(fulldat, kod %in% weird_key)

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
epi_corr_units <- epi %>%
  mutate(count_l = m2_to_l(count, interval = lower_layer - upper_layer)) %>%
  select(-count) %>%
  ## Filter out data that's listed as being from 1908 -- this is clearly a
  ## mistake
  filter(as.Date(date) >= as.Date("1945-01-01"))

## count_l has negatives (when upper layer > lower layer), and Inf (one one instance, where upper layer = lower layer)

#################################
####  View lifestage groups  ####
#################################

## List of life stage groups for Epischura. Eventually we'll need this
## information to group the Epischura data by stage.
unique(epi_corr_units$lifestage_gen)    # juv vs adult
unique(epi_corr_units$lifestage_cop)    # naup, copep, adult
unique(epi_corr_units[, c("code", "lifestage_cop")]) %>% arrange(lifestage_cop, code)





#############################################
####  Average temperature data by month  ####
#############################################

temp <- read.csv(paste0(dir, "temp_chl_secchi_wind/cleaned_data/temp_cleaned.csv"),
                 stringsAsFactors = FALSE)

temp_month <- temp %>%
  ## Create column of months
  mutate(monthyear = floor_date(as.Date(date), "month")) %>%
  ## Average temp by month
  group_by(monthyear) %>%
  summarize(temp = mean(temp, na.rm = TRUE))

#############################################
####  Average chlorophyll data by month  ####
#############################################

chla <- read.csv(paste0(dir, "temp_chl_secchi_wind/cleaned_data/chla_cleaned.csv"),
                 stringsAsFactors = FALSE)

chla_month <- chla %>%
  ## Create column of months
  mutate(monthyear = floor_date(as.Date(date), "month")) %>%
  ## Average chla by month
  group_by(monthyear) %>%
  summarize(chla = mean(chla, na.rm = TRUE))

####################################
####  View max depth histogram  ####
####################################

## How deep do observations usually go? What depth should we aggregate over for
## this analysis? Here's a histogram of max depth to think about.

epi_corr_units %>%
  group_by(date) %>%
  summarize(maxdepth = max(lower_layer, na.rm = TRUE)) %>%
  ggplot(aes(x = maxdepth)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Maximum sampling depth histogram: Epischura") +
  ggsave("../figs/epi_depth_hist.png")
#highest count for max depth is 500
#i.e., most epi counts are from max depth 500

## It's probably different for temp and chla

temp %>%
  group_by(date) %>%
  summarize(maxdepth = max(depth)) %>%
  ggplot(aes(x = maxdepth)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Maximum sampling depth histogram: temperature") +
  ggsave("../figs/temp_depth_hist.png")
#highest count for maxdepth is 50
#i.e., temp counts max depth usually only 50

chla %>%
  group_by(date) %>%
  summarize(maxdepth = max(depth)) %>%
  ggplot(aes(x = maxdepth)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Maximum sampling depth histogram: chlorophyll") +
  ggsave("../figs/chla_depth_hist.png")
#max depth usually 150
#i.e., chla max depth usually only 150
#so, epi counts are from deepest, then chla, then temp
#if want all, probably stay within 50m, since that's as deep as temp usually goes


################################################
####  View epi abundance and chla by depth  ####
################################################

## View average Epischura abundance by depth interval (note that overlapping
## depth intervals are considered separately with this approach -- the values
## for depth intervals 0-10, 10-25, and 25-50 do NOT go into 0-50; the only
## observations considered in 0-50 are the ones that came from samples of the
## whole 0-50 layer).
epi_depth_abund <- epi_corr_units %>%
  group_by(upper_layer, lower_layer) %>%
  summarize(mean_count = mean(count_l, na.rm = TRUE),
            n = length(count_l)) %>%
  ungroup() %>%
  arrange(lower_layer) %>%
  mutate(layer = paste(upper_layer, lower_layer, sep = "-")) %>% 
  as.data.frame()

#want to tag ones that are sequential vs overlapping
#e.g., 0-150 is overlapping whereas 0-10, 10-25, etc. is sequential
#and the ones that are negative or Inf
#grrr I know I did this before, stupid hard drive wipe...

sequential <- c("0-10", "10-25", "25-50", "50-100", "100-150", "150-250", "250-500")

epi_depth_abund <- epi_depth_abund %>% 
  mutate(tag = ifelse(layer %in% sequential, "sequential",
                      ifelse(mean_count < 0 | is.infinite(mean_count), "weird", "overlapping")))
#epi_depth_abund %>% arrange(tag, layer)

epi_depth_abund$layer <- ordered(epi_depth_abund$layer, levels = epi_depth_abund$layer)

## Plot these average abundances
ggplot(epi_depth_abund, aes(x = layer, y = mean_count)) +
  geom_bar(stat = "identity", aes(fill = tag)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Rather awkward representation of Epischura abundance by depth")

## Well there's something weird here: a row where upper_layer is 250,
## lower_layer is 150, and abundance is negative. View this data:
epi_corr_units %>%
  filter(upper_layer == 250 & lower_layer == 150)

## This issue comes up on on 2002-06-13. Look at the original data for that date:
unique(epi[epi$date == as.Date("2002-06-13"), c("date", "upper_layer", "lower_layer")])
## Okay, it looks like the upper and lower layer got switched. There is data
## from 100 to 150 and 250 to 500, so there should be 150 to 250. Ughhhhhh. And
## no one caught this before. Sigh.

## Another problem is a row where the upper and lower layers are both 150.

## I'll have to fix the above but for now I will just remove rows where
## upper_layer is greater than or equal to lower_layer
epi_depth_abund %>%
  filter(upper_layer < lower_layer) %>%
  ggplot(aes(x = layer, y = mean_count)) +
  geom_bar(stat = "identity", aes(fill = tag)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Rather awkward representation of Epischura abundance by depth") +
  ylab("Mean Epischura count (individuals/liter)") +
  xlab("Depth interval (m)") +
  ggsave("../figs/epi_avg_count_by_depth.png")

## Another thing we could do is just group Epischura into above 50 m and below
## 50 m and compare averages
epi_corr_units %>%
  filter(upper_layer > 50 | lower_layer < 50) %>%
  mutate(over_under_50 = ifelse(upper_layer > 50, "over 50", "under 50")) %>%
  ## Remove infinite value in count_l - should do this more permanently later
  mutate(count_l = ifelse(is.infinite(count_l), NA, count_l)) %>%
  group_by(over_under_50) %>%
  summarize(mean_count = mean(count_l, na.rm = TRUE))

## View mean chlorophyll by depth
chla_by_depth <- chla %>%
  group_by(depth) %>%
  summarize(mean_chla = mean(chla, na.rm = TRUE))

## Plot only <= 100 m depth because with the greater depths it becomes impossible to see anything
ggplot(chla_by_depth[chla_by_depth$depth <= 250, ], aes(x = depth, y = mean_chla)) +
  geom_bar(stat = "identity") +
  ggtitle("Average chlorophyll at depths <= 250") +
  ggsave("../figs/chla_depth_avg.png")

## Table of the few deepest chlorophyll means
filter(chla_by_depth, depth > 250) %>%
  data.frame()

########################################
####  View commonly sampled depths  ####
########################################

## Commonly sampled depths for chlorophyll
chla %>%
  group_by(depth) %>%
  tally() %>%
  print(n = 60)

## Commonly sampled depths for temperature
temp %>%
  group_by(depth) %>%
  tally() %>%
  print()

## Commonly sampled depths for Epischura - first create single column for depth
## layer
epi_layers <- epi_corr_units %>%
  arrange(upper_layer) %>%
  mutate(depth = paste(upper_layer, lower_layer, sep = "-"))

## Then make the column an ordered factor so the histogram bars are displayed in
## order by upper_layer
epi_layers$depth <- ordered(epi_layers$depth, levels = unique(epi_layers$depth))

## Plot
epi_layers %>%
  mutate(tag = ifelse(depth %in% sequential, "sequential",
                      ifelse(count_l < 0 | is.infinite(count_l), "weird", "overlapping"))) %>% 
  group_by(depth) %>%
  ggplot(aes(x = depth, fill = tag)) +
  geom_bar() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Sampling frequency of depth layers for Epischura") +
  ggsave("../figs/epi_sampling_depth_freq.png")

#############################################
####  Find out most abundant phyto taxa  ####
#############################################

## Steps:

## - Separate "unidentified" taxa into pico, nano, and unknown

## - For each date and depth, find which genera constitute >= 10% of abundance

## - Create a vector of these genera

## - Subset overall data set to only these genera

## - Subset data to Jan/Feb/Mar and Jul/Aug/Sep, average abundances of each
## genera, show top 10

phy <- read.csv(paste0(dir, "phyto/Baikal_Phyto_zeroesInc_2countRemoved_KeyUpdate_20120728.csv"), stringsAsFactors = FALSE)
names(phy) <- tolower(names(phy))
phy$date <- as.Date(phy$date, "%m/%d/%Y")

## Only keep data from >=1975 due to change in algae preservation technique
phy <- filter(phy, year >= 1975)

## Separate "unidentified" taxa into pico, nano, and unknown groups -- first
## figure out which codes from Lyubov's key correspond to which group

unid_group <- function(x) {
  ## Function to group unidentified genera into picoplankton, nanoplankton, and
  ## unknown using size information from Lyubov's plankton key
  if (x %in% c(526, 9526, 586, 587)) {  # pico
    result <- "unid_pico"
  } else if (x %in% c(527, 529, 583, 483, 550, 8550, 9550, 551, 8551, 552, 8552,
                      553, 8553, 554, 584, 580, 581, 582, 535, 536, 537, 538,
                      9538, 539, 9539, 540, 542, 543, 557, 570, 571, 572, 573,
                      574, 532, 531, 563, 564, 560, 561, 562)) { # nano
    result <- "unid_nano"
  } else if (x %in% c(1, 119, 120, 141, 193, 206, 211, 236, 237, 250, 263, 414,
                      443, 484, 521, 530, 533, 534, 549, 565, 566, 569, 575,
                      579, 588, 589, 590, 593, 594, 595, 596, 597, 598, 599,
                      8120, 8549, 8566, 9533, 9534, 9551)) { # unknown
    result <- "unid_unid"
  } else {
    result <- unique(phy[phy$code == x, "genus"])
  }
  return(result)
}

## Apply the funciton above to create a revised genus column that separates
## "unidentified" into "unid_pico", "unid_nano", and "unid_unid". This takes a
## long time, it needs to run overnight.
# phy_newgen <- phy %>%
#   rowwise() %>%
#   mutate(genus_revised = unid_group(code))

## ghastly results of system.time() on the above:
##     user   system  elapsed 
## 16078.62 10263.24 26544.29 

# ## Export to CSV so I don't have to rerun this every time.
# write.csv(phy_newgen,
#           "../data/phy_with_three_unidentified_groups.csv",
#           row.names = FALSE)

phy_newgen <- read.csv("../data/phy_with_three_unidentified_groups.csv",
                       stringsAsFactors = FALSE)

## For each date and depth, find which genera constitute >= 10% of abundance
phy_sml <- phy_newgen %>%
  ## Remove old genus column so I don't accidentally use it
  select(-genus) %>%
  ## Total phytoplankton density for each date and depth
  group_by(date, depth) %>%
  mutate(dens_sum = sum(density)) %>%
  ## Sum density for each genus at each date and depth
  group_by(date, depth, genus_revised) %>%
  summarize(density = sum(density),
            ## Need to keep the dens_sum column for later. The below is kind of
            ## a hacky way to do that, but it should be fine because there
            ## should be only one value of dens_sum for each date and depth. If
            ## this errors, something has gone wrong before this point.
            dens_sum = unique(dens_sum)) %>%
  ## Keep genera that make up >=10% of the total phytoplankton density
  # *for any depth and date*
  filter(density >= dens_sum / 10)

## List of all genera that at some point constitute >=10% of abundance
genera <- unique(phy_sml$genus_revised)

## Top 10 winter and summer genera above 150 m
#SL: adjusted to be 50 m - match lowest temp depth
phy_final <- phy_newgen %>%
  filter(genus_revised %in% genera &
           month %in% c(1, 2, 3, 7, 8, 9) &
           #depth <= 150) %>%
           depth <= 50) %>% 
  ## First sum within genera for each date and depth
  group_by(date, depth, genus_revised) %>%
  summarize(density = sum(density)) %>%
  ## Then average across dates and depths
  mutate(season = ifelse(month(date) %in% c(1, 2, 3), "winter", "summer")) %>%
  group_by(season, genus_revised) %>%
  summarize(density = mean(density)) %>%
  top_n(n = 10, wt = density) %>%
  arrange(desc(density))

## Convert revised genus to factor
phy_final$genus_revised <- factor(phy_final$genus_revised, levels = unique(phy_final$genus_revised))

## Plot
ggplot(phy_final, aes(x = genus_revised, y = density)) +
  facet_grid(. ~ season, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Most abundant phytoplankton genera (Jul/Aug/Sep vs. Jan/Feb/Mar) - 50 m") +
  ggsave(paste0("../figs/phyto_genera_winter_summer_", Sys.Date(), ".png"),
         width = 10, height = 7)

## Different version with winter defined as Feb/Mar/April/May instead of
## Jan/Feb/Mar. This may be more appropriate for Baikal's seasons.
phy_feb <- phy_newgen %>%
  filter(genus_revised %in% genera &
           month %in% c(2, 3, 4, 5, 7, 8, 9) &
           depth <= 50) %>% 
  #depth <= 150) %>%
  ## First sum within genera for each date and depth
  group_by(date, depth, genus_revised) %>%
  summarize(density = sum(density)) %>%
  ## Then average across dates and depths
  mutate(season = ifelse(month(date) %in% c(2, 3, 4, 5), "winter", "summer")) %>%
  group_by(season, genus_revised) %>%
  summarize(density = mean(density)) %>%
  top_n(n = 10, wt = density) %>%
  arrange(desc(density))

phy_feb$genus_revised <- factor(phy_feb$genus_revised, levels = unique(phy_feb$genus_revised))

## Plot
ggplot(phy_feb, aes(x = genus_revised, y = density)) +
  facet_grid(. ~ season, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Most abundant phytoplankton genera (Jul/Aug/Sep vs. Feb/Mar/Apr/May) - 50 m") +
  ggsave(paste0("../figs/phyto_genera_winter_summer_FMAM_", Sys.Date(), ".png"),
         width = 10, height = 7)

## Genera in phy_feb not present in phy_final == Stephanodiscus
as.character(phy_feb[phy_feb$season == "winter", ]$genus_revised)[!as.character(phy_feb[phy_feb$season == "winter", ]$genus_revised) %in% as.character(phy_final[phy_final$season == "winter", ]$genus_revised)]

## Genera in phy_final not present in phy_feb == Aulacoseira
as.character(phy_final[phy_final$season == "winter", ]$genus_revised)[!as.character(phy_final[phy_final$season == "winter", ]$genus_revised) %in% as.character(phy_feb[phy_feb$season == "winter", ]$genus_revised)]


## ----> Same, but for species

## Common winter and summer species
phy_sp_common <- phy_newgen %>%
  ## Keep only common genera, months of interest, and depths of interest
  filter(genus_revised %in% genera &
           month %in% c(1, 2, 3, 7, 8, 9) &
           depth <= 50) %>% 
  ## Create genus_species column
  mutate(genus_species = paste(genus_revised, species, sep = ".")) %>% 
  ## First sum within genus_species for each date and depth
  ## (Even though shouldn't be necessary, there are a few instances where
  ## code is different but date, depth, and genus_species have duplicates;
  ## also especially necessary because of all the unid_unid)
  group_by(date, depth, genus_species) %>%
  summarize(density = sum(density)) %>%
  ## Then average across dates and depths
  mutate(season = ifelse(month(date) %in% c(1, 2, 3), "winter", "summer")) %>%
  group_by(season, genus_species) %>%
  summarize(density = mean(density)) %>%
  top_n(n = 10, wt = density) %>%
  arrange(season, desc(density)) %>% 
  as.data.frame()

## Convert revised genus_species to factor
phy_sp_common$genus_species <- factor(phy_sp_common$genus_species, levels = unique(phy_sp_common$genus_species))

## Plot
ggplot(phy_sp_common, aes(x = genus_species, y = density)) +
  facet_grid(. ~ season, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Most abundant phytoplankton species (Jul/Aug/Sep vs. Jan/Feb/Mar) - 50 m") +
  ggsave(paste0("../figs/phyto_species_winter_summer_", Sys.Date(), ".png"),
         width = 10, height = 7)

## Different version with winter defined as Feb/Mar/April/May instead of
## Jan/Feb/Mar. This may be more appropriate for Baikal's seasons.

## Common winter and summer species - alt winter/summer seasons
phy_sp_common_feb <- phy_newgen %>%
  ## Keep only common genera, months of interest, and depths of interest
  filter(genus_revised %in% genera &
           month %in% c(2, 3, 4, 5, 7, 8, 9) &
           depth <= 50) %>% 
  ## Create genus_species column
  mutate(genus_species = paste(genus_revised, species, sep = ".")) %>% 
  ## First sum within genus_species for each date and depth
  ## (Even though shouldn't be necessary, there are a few instances where
  ## code is different but date, depth, and genus_species have duplicates;
  ## also especially necessary because of all the unid_unid)
  group_by(date, depth, genus_species) %>%
  summarize(density = sum(density)) %>%
  ## Then average across dates and depths
  mutate(season = ifelse(month(date) %in% c(2, 3, 4, 5), "winter", "summer")) %>%
  group_by(season, genus_species) %>%
  summarize(density = mean(density)) %>%
  top_n(n = 10, wt = density) %>%
  arrange(season, desc(density)) %>% 
  as.data.frame()

## Convert revised genus_species to factor
phy_sp_common_feb$genus_species <- factor(phy_sp_common_feb$genus_species, levels = unique(phy_sp_common_feb$genus_species))

## Plot
ggplot(phy_sp_common_feb, aes(x = genus_species, y = density)) +
  facet_grid(. ~ season, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Most abundant phytoplankton species (Jul/Aug/Sep vs. Feb/Mar/Apr/May) - 50 m") +
  ggsave(paste0("../figs/phyto_species_winter_summer_FMAM_", Sys.Date(), ".png"),
         width = 10, height = 7)


##############################################################
####  Find out more information on "Unidentified" phytos  ####
##############################################################

## Figure out which groups are represented by "Unidentified"
unid <- phy %>%
  filter(genus == "Unidentified") %>%
  select(code, group, genus, species) %>%
  unique()

unique(unid$group)
## [1] "Cyano"        "Green"        "Chryso"       "Bacteria"     "Flagellate"  
## [6] "PicoAlgae"    "NanoAlgae"    "Cysts"        "Unidentified" "Alga"        

## I checked the phyto key and there's no more information about these
## unidentified taxa there. Let's see which are most abundant

unid_abund <- phy %>%
  filter(genus == "Unidentified") %>%
  ## Less than 150m depth
  filter(depth <= 150) %>%
  ## Sum within group/date/depth
  group_by(date, depth, group) %>%
  summarize(density = sum(density)) %>%
  ## Average by group
  group_by(group) %>%
  summarize(density = mean(density)) %>% 
  arrange(desc(density))

#################################################
####  Limit depth across data - exploratory  ####
#################################################

hist(temp$depth)
hist(chla$depth)
hist(phy_sml$depth) 
#phyto keeps data from >=1975 due to change in algae preservation technique
#limit other as well to match

head(epi_corr_units)

temp50 <- filter(temp, depth <= 50 & year(as.Date(date))>= 1975)
chla50 <- filter(chla, depth <= 50 & year(as.Date(date))>= 1975)
#should really redo so keep any genera >=10% *for depths <= 50m*
phy_sml50 <- filter(phy_sml, depth <= 50) %>% as.data.frame() #already limited to >=1975

epi_corr_units50 <- epi_corr_units %>% 
  filter(upper_layer <= 50 & lower_layer <= 50 & year(as.Date(date))>= 1975) %>% 
  mutate(layers = paste(upper_layer, lower_layer, sep = "-"))

unique(epi_corr_units50$layers) %>% sort()
# [1] "0-10"  "0-25"  "0-50"  "10-25" "25-50"
#so, have sequential: 0-10, 10-25, 25-50; and overlap 0-25, 0-50

#check out counts by depth layers
epi_corr_units50 %>% 
  group_by(layers) %>% 
  summarize(min_count = min(count_l), 
            max_count = max(count_l),
            mean_count = mean(count_l),
            median_count = median(count_l))
#yeah, there's something up with 0-25, I don't trust it

ggplot(epi_corr_units50, aes(layers, count_l)) +
  geom_boxplot() +
  ylim(0, 5)

ggplot(epi_corr_units50, aes(layers, count_l)) +
  geom_violin() +
  geom_jitter()

#SH: I agree re depth and year. 1975 it is, and 50m. 




#########################################
####  Grouping phyto and epi stages  ####
#########################################

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

#1-12 are phyto, 13-16 are zoop

unique(phy_newgen$genus_revised) %>% sort()


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

head(phy_dat_grouped)

## ----> zoop groups

head(epi_corr_units) #only Epischura - will need to get cyclops data from elsewhere

epi_grouped <- epi_corr_units %>% 
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

# filter(epi_corr_units, date == "1975-01-07" & lifestage_cop == "naup" & upper_layer == 0 & lower_layer == 10)
# filter(epi_groups, date == "1975-01-07" & lifestage_cop == "naup" & upper_layer == 0 & lower_layer == 10)
# 
# filter(epi_corr_units, date == "2003-12-06" & lifestage_cop == "copep" & upper_layer == 10 & lower_layer == 25)
# filter(epi_groups, date == "2003-12-06" & lifestage_cop == "copep" & upper_layer == 10 & lower_layer == 25) 

head(cyclop_corr_units)

cyclop_grouped <- cyclop_corr_units %>% 
  mutate(layer_group = paste(upper_layer, lower_layer, sep = "-")) %>% 
  #for now, sticking with sequential depths to 50 m
  filter(layer_group %in% c("0-10", "10-25", "25-50")) %>% 
  mutate(layer_group = ordered(layer_group, levels = c("0-10", "10-25", "25-50"))) %>% 
  #keep only cyclops adults
  filter(lifestage_cop == "adult") %>% 
  #group and aggregate counts within cyclops adults
  group_by(date, layer_group, upper_layer, lower_layer, genus, lifestage_cop) %>% 
  summarize(count_l_sum = sum(count_l)) %>% 
  #keep only since 1975
  mutate(year = year(as.Date(date))) %>% 
  filter(year >= 1975) %>% 
  as.data.frame()

# filter(cyclop_corr_units, date == "1975-03-10" & lifestage_cop == "adult" & upper_layer == 0 & lower_layer == 10)
# filter(cyclop_grouped, date == "1975-03-10" & lifestage_cop == "adult" & upper_layer == 0 & lower_layer == 10)
# 
# filter(cyclop_corr_units, date == "1994-11-09" & lifestage_cop == "adult" & upper_layer == 25 & lower_layer == 50)
# filter(cyclop_grouped, date == "1994-11-09" & lifestage_cop == "adult" & upper_layer == 25 & lower_layer == 50)


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
head(epi_grouped)
head(cyclop_grouped)

# ----> merge - match date and depth

#merge and make long - NAs ok - can interpolate later

#temp and chla (env vars)
env_merge <- chla_match %>% 
  select(-month, -Year) %>% 
  merge(temp_match, by = c("date", "depth"), all = TRUE)

#merge with phyto
phyto_env_merge <- phy_dat_grouped %>% 
  select(-year, -month) %>% 
  merge(env_merge, by = c("date", "depth"), all = TRUE) %>% 
  select(date, depth, chla, temp, genus_groups, density_genus_sum)

# #epi and cyclops
# epi_test <- epi_grouped %>% select(date, layer_group) %>% unique() %>% arrange(date, layer_group)
# cyclops_test <- cyclop_grouped %>% select(date, layer_group) %>% unique() %>% arrange(date, layer_group)
# identical(epi_test, cyclops_test) #TRUE
# #should be ok to stack (makes sense, since from same source data frame anyways)

pred_stack <- rbind(epi_grouped, cyclop_grouped) %>% select(-year)

# ----> check these out

head(phyto_env_merge)
head(pred_stack)

#need to approximate layer group for env/phyto
#well, well, well...I have no idea whether the layer groups are range inclusive
#is 0-10 up to and including 10? is 0-10 <=10 and then 10-25 is >10?

#for now, going to go with <=x and >x
#since need 50 included as <=50...

phyto_env_merge_layers <- phyto_env_merge %>% 
  mutate(approx_layer = ifelse(depth <= 10, "0-10", NA),
         approx_layer = ifelse(depth >10 & depth <= 25, "10-25", approx_layer),
         approx_layer = ifelse(depth > 25 & depth <= 50, "25-50", approx_layer))

# phyto_env_merge_layers %>% select(depth, approx_layer) %>% unique() %>% arrange(depth)

head(phyto_env_merge_layers)
head(pred_stack)

pred_stack_ready <- pred_stack %>% 
  select(-upper_layer, -lower_layer)

phyto_env_ready <- phyto_env_merge_layers %>% 
  select(-depth)

full_merge <- merge(phyto_env_ready, pred_stack_ready,
                    by.x = c("date", "approx_layer"),
                    by.y = c("date", "layer_group")) %>% 
  rename(phyto_genera = genus_groups,
         zoop_genera = genus,
         zoop_lifestage = lifestage_cop)

#hmm...leave phyto and zoop genera in separate columns or combine?
#count vs density 
#would retain both and then have lots of NAs...
#thinking remain separate since "phyto" vs "zoop"
#sort of like "chla" vs "temp" rather than "variable" vs "value"
#when value has variable (ha) units

#need to reshape so genera along top
#will need to be very careful since phyto has density sum and zoop has count sum...

#example MAR data
data(L4.AllDates)
head(L4.AllDates)

head(full_merge)

# ----> reshape so genera along top

full_merge_wide1 <- dcast(full_merge, 
                          date + approx_layer + chla + temp + zoop_genera +
                            zoop_lifestage + count_l_sum ~ phyto_genera, value.var = "density_genus_sum")

#ah yes, should do this above so doesn't barf...

head(phyto_env_merge_layers)

phyto_env_wide <- dcast(phyto_env_merge_layers,
                        date + depth + chla + temp + approx_layer ~ genus_groups, value.var = "density_genus_sum")

#gah, even further back to avoid the unnecessary NAs...
head(phy_dat_grouped)

phy_dat_grouped_wide <- phy_dat_grouped %>% 
  select(-year, -month) %>% 
  dcast(date + depth ~ genus_groups, value.var = "density_genus_sum")

head(pred_stack)

pred_stack_wide <- pred_stack %>% 
  mutate(genus_stage = paste(genus, lifestage_cop, sep = "_")) %>% 
  select(-genus, -lifestage_cop) %>% 
  dcast(date + layer_group + upper_layer + lower_layer ~ genus_stage, value.var = "count_l_sum")

head(phy_dat_grouped_wide)
head(pred_stack_wide)
head(env_merge)

#merge env with phyto

phyto_env_wide_merge <- merge(phy_dat_grouped_wide, env_merge,
                              by = c("date", "depth"))

