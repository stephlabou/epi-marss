###########################################
####  Prepare data for MARSS analysis  ####
###########################################

## This script loads data and examines various aspects to see which data we
## should use for analysis. Eventually it will output data formatted for use
## with MARSS.

library("dplyr")
library("ggplot2")
library("lubridate")

## Need to have all life stages of Epischura, minus those that are known to be
## double counts

## Directory of the long-term data Subversion repository. UPDATE THIS to run the
## code on your own machine.
dir <- "/Users/Kara/projects/baikal_svn/Longterm_data/"

## Load the main zooplankton data file and make column names lower case Note:
## this is the data with zeroes already in place when species codes weren't
## observed
fulldat <- read.csv(paste0(dir, "zoo/data/zoopzeroskey_alldepths.csv"),
                    stringsAsFactors = FALSE)
colnames(fulldat) <- tolower(colnames(fulldat))

## Load zooplankton key and make column names lower case
key <- read.csv(paste0(dir, "zoo/data/key.csv"), stringsAsFactors = FALSE)
colnames(key) <- tolower(colnames(key))

########################################################
####  Keep only Epischura and remove double counts  ####
########################################################

## Codes that represent Epischura non-double counts according to Derek's key:
nodbl <- filter(key, genus == "Epischura" & species == "baicalensis"
                & doublecount %in% c("N", "M")) %>% .$kod

## Keep only the Epischura non-double counts
epi <- filter(fulldat, kod %in% nodbl) %>%
  ## Rename some columns
  rename(code = kod, date = date, upper_layer = ver_gr, lower_layer = nig_gr)

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

#################################
####  View lifestage groups  ####
#################################

## List of life stage groups for Epischura. Eventually we'll need this
## information to group the Epischura data by stage.
unique(epi_corr_units$lifestage_gen)    # juv vs adult
unique(epi_corr_units$lifestage_cop)    # naup, copep, adult
unique(epi_corr_units[, c("code", "lifestage_cop")])

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

## It's probably different for temp and chla

temp %>%
  group_by(date) %>%
  summarize(maxdepth = max(depth)) %>%
  ggplot(aes(x = maxdepth)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Maximum sampling depth histogram: temperature") +
  ggsave("../figs/temp_depth_hist.png")

chla %>%
  group_by(date) %>%
  summarize(maxdepth = max(depth)) %>%
  ggplot(aes(x = maxdepth)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Maximum sampling depth histogram: chlorophyll") +
  ggsave("../figs/chla_depth_hist.png")

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
  summarize(mean_count = mean(count_l, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(lower_layer) %>%
  mutate(layer = paste(upper_layer, lower_layer, sep = "-"))

epi_depth_abund$layer <- ordered(epi_depth_abund$layer, levels = epi_depth_abund$layer)

## Plot these average abundances
ggplot(epi_depth_abund, aes(x = layer, y = mean_count)) +
  geom_bar(stat = "identity") +
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
  geom_bar(stat = "identity") +
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
  summarize(mean_count = mean(count_l, na.rm = TRUE, ))

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
  group_by(depth) %>%
  ggplot(aes(x = depth)) +
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
phy_newgen <- phy %>%
  rowwise() %>%
  mutate(genus_revised = unid_group(code))

## ghastly results of system.time() on the above:
##     user   system  elapsed 
## 16078.62 10263.24 26544.29 

## Export to CSV so I don't have to rerun this every time.
write.csv(phy_newgen,
          "../data/phy_with_three_unidentified_groups.csv",
          row.names = FALSE)

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
  filter(density >= dens_sum / 10)

## List of all genera that at some point constitute >=10% of abundance
genera <- unique(phy_sml$genus_revised)

## Top 10 winter and summer genera above 150 m
phy_final <- phy_newgen %>%
  filter(genus_revised %in% genera &
         month %in% c(1, 2, 3, 7, 8, 9) &
         depth <= 150) %>%
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
  ggtitle("Most abundant phytoplankton genera (Jul/Aug/Sep vs. Jan/Feb/Mar)") +
  ggsave(paste0("../figs/phyto_genera_winter_summer_", Sys.Date(), ".png"),
         width = 10, height = 7)

## Different version with winter defined as Feb/Mar/April/May instead of
## Jan/Feb/Mar. This may be more appropriate for Baikal's seasons.
phy_feb <- phy_newgen %>%
  filter(genus_revised %in% genera &
         month %in% c(2, 3, 4, 5, 7, 8, 9) &
         depth <= 150) %>%
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
  ggtitle("Most abundant phytoplankton genera (Jul/Aug/Sep vs. Feb/Mar/Apr/May)") +
  ggsave(paste0("../figs/phyto_genera_winter_summer_FMAM_", Sys.Date(), ".png"),
         width = 10, height = 7)

## Genera in phy_feb not present in phy_final
as.character(phy_feb[phy_feb$season == "winter", ]$genus)[!as.character(phy_feb[phy_feb$season == "winter", ]$genus) %in% as.character(phy_final[phy_final$season == "winter", ]$genus)]

## Genera in phy_final not present in phy_feb
as.character(phy_final[phy_final$season == "winter", ]$genus)[!as.character(phy_final[phy_final$season == "winter", ]$genus) %in% as.character(phy_feb[phy_feb$season == "winter", ]$genus)]

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
  filter(depth <= 150) %>%
  ## Sum within group/date/depth
  group_by(date, depth, group) %>%
  summarize(density = sum(density)) %>%
  ## Average by group
  group_by(group) %>%
  summarize(density = mean(density))
