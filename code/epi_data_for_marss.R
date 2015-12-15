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

