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

## Directory of the long-term data repo on my machine. Obviously this will need
## to be changed if anyone else wants to run the code, but for right now it's
## just me...
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
  stopifnot(is.numeric(x))              # count should be numeric
  ## takes units in 1000 individuals/m2 and converts to individuals per liter
  individuals <- x * 1000 # convert to individuals/m2
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

## Eventually we'll need this information to group the Epischura data by stage
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

epi_depth_abund <- epi_corr_units %>%
  group_by(upper_layer, lower_layer) %>%
  summarize(mean_count = mean(count_l, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(lower_layer) %>%
  mutate(layer = paste(upper_layer, lower_layer, sep = "-"))

epi_depth_abund$layer <- ordered(epi_depth_abund$layer, levels = epi_depth_abund$layer)

ggplot(epi_depth_abund, aes(x = layer, y = mean_count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Rather awkward representation of Epischura abundance by depth")

## Well there's something weird here: a row where upper_layer is 250,
## lower_layer is 150, and abundance is negative
epi_corr_units %>%
  filter(upper_layer == 250 & lower_layer == 150)

## This happens on 2002-06-13
unique(epi[epi$date == as.Date("2002-06-13"), c("date", "upper_layer", "lower_layer")])
## Okay, it looks like the upper and lower layer got switched. There is data
## from 100 to 150 and 250 to 500, so there should be 150 to 250. Ughhhhhh. And
## no one caught this before. Sigh.

## Another problem is a row where the upper and lower layers are both 150.
## Remove this for now from the plot; fix later

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

## Chlorophyll
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
  geom_histogram(stat = "bin") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Sampling frequency of depth layers for Epischura") +
  ggsave("../figs/epi_sampling_depth_freq.png")

#############################################
####  Find out most abundant phyto taxa  ####
#############################################

## Steps:

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

## For each date and depth, find which genera constitute >= 10% of abundance
phy_sml <- phy %>%
  group_by(date, depth) %>%
  ## Sum density for each date and depth
  mutate(dens_sum = sum(density)) %>%
  ## Sum density for each genus at each date and depth
  group_by(date, depth, genus) %>%
  summarize(density = sum(density),
            ## Need to keep the dens_sum column for later. The below is kind of
            ## a hacky way to do that, but it should be fine because there
            ## should be only one value of dens_sum for each date and depth. If
            ## this errors, something has gone wrong before this point.
            dens_sum = unique(dens_sum)) %>%
  ## Keep genera that make up >=10% of the density
  filter(density >= dens_sum / 10)

## List of all genera that at some point constitute >=10% of abundance
genera <- unique(phy_sml$genus)

## Top 10 winter and summer genera above 150 m
phy_final <- phy %>%
  filter(genus %in% genera &
         month %in% c(1, 2, 3, 7, 8, 9) &
         depth <= 150) %>%
  ## First sum within genera for each date and depth
  group_by(date, depth, genus) %>%
  summarize(density = sum(density)) %>%
  ## Then average across dates and depths
  mutate(season = ifelse(month(date) %in% c(1, 2, 3), "winter", "summer")) %>%
  group_by(season, genus) %>%
  summarize(density = mean(density)) %>%
  top_n(n = 10, wt = density) %>%
  arrange(desc(density))

phy_final$genus <- factor(phy_final$genus, levels = unique(phy_final$genus))

## Plot
ggplot(phy_final, aes(x = genus, y = density)) +
  facet_grid(. ~ season, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ggtitle("Most abundant phytoplankton genera (Jul/Aug/Sep vs. Jan/Feb/Mar)") +
  ggsave(paste0("../figs/phyto_genera_winter_summer_", Sys.Date(), ".png"),
         width = 10, height = 7)

