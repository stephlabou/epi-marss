## Once data is organized into monthly 0-50m averages/totals
## need to limit/subset data to best possible

# ----> how to deal with NA values?
#from SH:
# In the past I have interpolated when it is just one month missing, and not two in a row. 
# Where missing more than 2 in a row - in the past I used averages for the whole time period for that month. 
# If it's missing more than 2 in a row for any variable, the MAR is capable of just skipping 
# those months that have missing values. But we want to know when it is happening.

# So: interpolate when missing one NA, use average when 2 or more

# [X] tag rows in columns where missing 2 or more in a row
# [X] find averages for whole time period for that month
# [X] add column with these avgs inserted in places where missing 2 or more
# (then can make decisions about what to do for these in MARSS - have options)

# ----> handling real zeros
# For real zeroes, MAR also will blow up on zeroes if there are too many in a row... 
# a few is not a big deal, as I recall. 
# So the rule of thumb is to calculate half the lowest observable (or observed non-zero) value, 
# plus or minus some random number, to replace zeroes.  For now, I would say leave real zeroes in 
# and we can always transform them at the time of analysis.

# SL: doing zero replacements now to be ready for analysis

# [X] for each variable, find 1/2 lowest non-zero value + random number
# [X] new column with these replacing zero values

#======================================================#
## =========== load packages, read in data ============#
#======================================================#

#load packages
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

#read in data, output from "organize_data_marss_monthly.R"
dat <- read.csv("../data/mar_dat.csv")

#======================================================#
## ========= subset to most appropriate years =========#
#======================================================#

summary(dat)

#When merged chla/temp/phyto/zoop in marss monthly script,
#kept all - but this leads to NA values

#chla: 1979-2008
#temp: 1975-2008
#zoop: 1975-2003
#phyto: 1975-1999

#will limit to 1979-1999 (inclusive) - 20 years when all vars of interest
#will still have missing, but this will limit instances

#first, want to be sure every month and year is represented
yrs <- seq(1975, 2008, 1)
mons <- seq(1, 12, 1)
full_yr_mon <- expand.grid(yrs, mons) %>% rename(year = Var1, month = Var2)

full_dates <- merge(dat, full_yr_mon, by = c("year", "month"), all = TRUE)

#now filter to be 1979-1999 and force chronological
dat_yrs <- full_dates %>% 
  filter(year >= 1979 & year <= 1999) %>% 
  arrange(year, month)

#check result
summary(dat_yrs)
#a more achievable number of NAs, mostly for chla

#======================================================#
##============= tag runs of missing data =============##
#======================================================#

#I have proved to myself that this works

source("Functions/rle_updated.R")
source("Functions/tag_runs_incl_na.R")

# hmm, any instances where chla NA for avg and not for sum?
# shouldn't be
# filter(dat_yrs, is.na(chla_050depth_avg) & !is.na(chla_050depth_sum)) 
# 0 rows
# ok, so can tag chla based on chla avg columns

#all the missing zoops are from same date (most of 1991)
#yes, these are the 22 NA rows for zoop
#filter(dat_yrs, is.na(Synedra_density))

#only one temp column

#only one zoop NA - no need to tag and replace (?)
#filter(dat_yrs, is.na(Epischura_adult_count))

#tag runs for chla, temp, phyto, zoop

dat_runs <- dat_yrs %>% 
  #yes, this function format must be like this
  #small price to pay
  mutate(chla_runs = wrangle_na_runs_col(dat_yrs$chla_050depth_avg),
         temp_runs = wrangle_na_runs_col(dat_yrs$temp_050depth_avg),
         #pick any phyto - all have NAs on same date
         phyto_runs = wrangle_na_runs_col(dat_yrs$Synedra_density),
         zoop_runs = wrangle_na_runs_col(dat_yrs$Epischura_adult_count)) %>% 
  select(year, month, chla_runs, temp_runs, phyto_runs, zoop_runs)

head(dat_runs)
summary(dat_runs)
#chla missing a while bunch, all of 1984 and some surrounding months
#(confirmed by looking back at orig data, missing all of 1984)

#============================================================#
##====== find monthly averages over whole time period ======##
#============================================================#

# dat_avgs <- dat_yrs %>%
#             group_by(month) %>%
#             summarize(temp_avg = mean(temp_050depth_avg, na.rm = TRUE),
#                       chla_avg_avg = mean(chla_050depth_avg, na.rm = TRUE),
#                       chla_sum_avg = mean(chla_050depth_sum, na.rm = TRUE)) %>%
#             as.data.frame()
# 
# #still need phyto avgs...
# test <- ddply(dat_yrs, "month", function(x) {
#   y <- subset(x, select = -c(month, year))
#   apply(y, 2, mean, na.rm = TRUE)
#   })
# 
# #dat_avgs: month 2, temp_avg = 0.7220759, month 11: chla_avg_avg = 0.5336223
# #test:     month 2, temp =     0.7220759 ; month 11: chla_avg =    0.5336223
# #yep, works
# 
# all.equal(dat_avgs$chla_sum_avg, test$chla_050depth_sum)
# identical(dat_avgs$chla_sum_avg, test$chla_050depth_sum)

#groups by month, finds mean for each column, excludes NA values
#(I have proved to myself that this works)
dat_avgs <- ddply(dat_yrs, "month", function(x) {
  y <- subset(x, select = -c(month, year))
  apply(y, 2, mean, na.rm = TRUE)
})

#rename so makes sense when merge
colnames(dat_avgs) <- paste("month_mean", colnames(dat_avgs), sep = "_")

head(dat_avgs)


#============================================================#
##========== updated df to replace missing values ==========##
#============================================================#

# Where missing more than 2 in a row - in the past I used averages for the whole time period for that month. 
# If it's missing more than 2 in a row for any variable, the MAR is capable of just skipping 
# those months that have missing values. But we want to know when it is happening.

head(dat_runs) #runs of missing
head(dat_avgs) #month avgs
head(dat_yrs) #main data

#merge
dat_merge <- merge(dat_yrs, dat_runs, by = c("year", "month"), all = TRUE) %>% arrange(year, month)
head(dat_merge)

#would this be easier to tag if maybe long form?
# 
# dat_long <- melt(dat_merge, id.vars = c("year", "month", "chla_runs", "temp_runs", "phyto_runs")) %>% 
#             arrange(year, month, variable) %>% 
#             rename(orig_var = variable, orig_value = value)
# 
# head(dat_long)
# head(dat_avgs)
# 
# dat_avgs_long <- melt(dat_avgs, id.vars = "month_mean_month") %>% 
#                   rename(month = month_mean_month,
#                          avg_month_variable = variable,
#                          avg_month_value = value)
# 
# head(dat_avgs_long)
# head(dat_long)

#split by groups?
dat_temp <- select(dat_merge, year, month, temp_050depth_avg, temp_runs)

dat_chla <- select(dat_merge, year, month, chla_050depth_avg, chla_050depth_sum, chla_runs)

dat_phyto <- dat_merge[,c(1:2, 10:(ncol(dat_merge)-4), (ncol(dat_merge)-1))]

dat_zoop <- select(dat_merge, year, month, Cyclops_adult_count, Epischura_adult_count,
                   Epischura_copep_count, Epischura_naup_count, zoop_runs)
#no runs for zoop, since only missing one ever (only one NA total)

#omg there has to be a better way
#save yourself, find a more efficient solution
avg_temp <- select(dat_avgs, month_mean_month, month_mean_temp_050depth_avg)
avg_chla <- select(dat_avgs, month_mean_month, 
                   month_mean_chla_050depth_avg, month_mean_chla_050depth_sum)
avg_phyto <- dat_avgs[, c(1, 9:ncol(dat_avgs))]
avg_zoop <- select(dat_avgs, month_mean_month,
                   month_mean_Cyclops_adult_count, month_mean_Epischura_adult_count,
                   month_mean_Epischura_copep_count, month_mean_Epischura_naup_count)

#combine these

# ----> temperature
head(dat_temp)
head(avg_temp)

full_temp <- merge(dat_temp, avg_temp, by.x = "month", by.y = "month_mean_month", all.x = TRUE) %>% arrange(year, month)

fixed_temp <- full_temp %>% 
  mutate(fixed_temp_avg = ifelse(is.na(temp_050depth_avg) & temp_runs >= 2, month_mean_temp_050depth_avg, temp_050depth_avg))

# head(fixed_temp)
# filter(fixed_temp, year == 1995)
# yep, works

#format for later stacking
fixed_temp_formatted <- fixed_temp %>% 
  mutate(variable = "temp_avg") %>% 
  select(year, month, temp_runs, variable, temp_050depth_avg, fixed_temp_avg) %>% 
  rename(orig_value = temp_050depth_avg, 
         fixed_value = fixed_temp_avg,
         num_runs = temp_runs)

# filter(fixed_temp_formatted, year == 1995)

# ----> chla
head(dat_chla)
head(avg_chla)

#make dat long
dat_chla_long <- dat_chla %>% 
  melt(id.vars = c("year", "month", "chla_runs")) %>% 
  rename(orig_value = value)

#make avgs long
avg_chla_long <-  melt(avg_chla, id.vars = "month_mean_month") %>% 
  mutate(variable = gsub("month_mean_", "", variable)) %>% 
  arrange(month_mean_month) %>% 
  rename(month = month_mean_month, 
         month_avg_value = value)

#merge
full_chla <- merge(dat_chla_long, avg_chla_long, by = c("month", "variable"), all = TRUE) %>% 
  arrange(year, month, variable) %>% 
  select(year, month, chla_runs, variable, orig_value, month_avg_value)

#fix values
fixed_chla <- full_chla %>% 
  mutate(fixed_value = ifelse(is.na(orig_value) & chla_runs >= 2, month_avg_value, orig_value))

#format for later stacking
fixed_chla_formatted <- fixed_chla %>% 
  select(-month_avg_value) %>% 
  rename(num_runs = chla_runs)

# ----> phyto
head(dat_phyto)
head(avg_phyto)

#make phyto data long
dat_phyto_long <- melt(dat_phyto, id.vars = c("year", "month", "phyto_runs")) %>%
  rename(orig_value = value)
head(dat_phyto_long)

#make monthly averages long
avg_phyto_long <- melt(avg_phyto, id.vars = "month_mean_month") %>% 
  mutate(variable = gsub("month_mean_", "", variable)) %>% 
  arrange(month_mean_month) %>% 
  rename(month = month_mean_month, 
         month_avg_value = value)
#head(avg_phyto_long)

#merge long forms together
full_phyto <- merge(dat_phyto_long, avg_phyto_long, by = c("month", "variable"), all = TRUE) %>% 
  arrange(year, month, variable) %>% 
  select(year, month, phyto_runs, variable, orig_value, month_avg_value)

#head(full_phyto, 20)

# filter(full_phyto, is.na(orig_value)) %>% head(20)
# filter(full_phyto, variable == "Romeria_density")

#phyto runs technically from Synedra
#so ONLY MAKES SENSE TO INTERPRET PHYTO RUNS WHEN VALUES ARE NA
#OTHERWISE THIS IS VERY MISLEADING SO BE SUPER CAREFUL ABOUT THIS
# filter(full_phyto, phyto_runs == 12 & variable == "Stephanodiscus_density")
# filter(full_phyto, phyto_runs == 12 & variable == "Synedra_density")

fixed_phyto <- full_phyto %>% 
  mutate(fixed_value = ifelse(is.na(orig_value) & phyto_runs >= 2, month_avg_value, orig_value))

# filter(fixed_phyto, is.na(orig_value) & phyto_runs > 2) %>% head(20)
# filter(fixed_phyto, year == 1991 & variable == "Romeria_density")
# filter(fixed_phyto, !is.na(orig_value)) %>% tail(20)
# filter(fixed_phyto, is.na(orig_value) & phyto_runs <= 2) %>% head(20)
# filter(fixed_phyto, year == 1984 & variable == "Synedra_density")

head(fixed_phyto)

#format for stacking
fixed_phyto_formatted <- fixed_phyto %>% 
  select(-month_avg_value) %>% 
  rename(num_runs = phyto_runs)

# ----> format zoop data
head(dat_zoop)
head(avg_zoop)

#make long
dat_zoop_long <- melt(dat_zoop, id.vars = c("year", "month", "zoop_runs")) %>% 
  rename(orig_value = value)

#make long
avg_zoop_long <- melt(avg_zoop, id.vars = "month_mean_month") %>% 
  mutate(variable = gsub("month_mean_", "", variable)) %>% 
  arrange(month_mean_month) %>% 
  rename(month = month_mean_month, 
         month_avg_value = value)

#merge
full_zoop <- merge(dat_zoop_long, avg_zoop_long, by = c("month", "variable"), all = TRUE) %>% 
  arrange(year, month, variable) %>% 
  select(year, month, zoop_runs, variable, orig_value, month_avg_value)

#fix NAs (not happening w/ zoops - more to get in proper format)
fixed_zoop <- full_zoop %>% 
  mutate(fixed_value = ifelse(is.na(orig_value) & zoop_runs >= 2, month_avg_value, orig_value))

#format for stacking
fixed_zoop_formatted <- fixed_zoop %>% 
  select(-month_avg_value) %>% 
  rename(num_runs = zoop_runs)

#filter(fixed_zoop_formatted, orig_value != fixed_value)
#fixed all same as orig, as expected, since no runs of NAs in zoop

# ----> stack these together
head(fixed_temp_formatted)
head(fixed_chla_formatted)
head(fixed_phyto_formatted)
head(fixed_zoop_formatted)

fixed_na_dat <- rbind(fixed_temp_formatted, fixed_chla_formatted, 
                      fixed_phyto_formatted, fixed_zoop_formatted) %>% 
  arrange(year, month, variable) %>% 
  rename(fixed_na_value = fixed_value)

# head(all_fixed_dat, 20)
# filter(all_fixed_dat, is.na(orig_value) & num_runs == 1) %>% head()
# yep, looking good

# can make wide again to run in MARSS when needed


#======================================================#
##=============== replace zero values ================##
#======================================================#

# So the rule of thumb is to calculate half the lowest observable (or observed non-zero) value, 
# plus or minus some random number, to replace zeroes.  For now, I would say leave real zeroes in 
# and we can always transform them at the time of analysis.

head(fixed_na_dat)

#find lowest non-zero for each varialbe

min_origvals <- fixed_na_dat %>% 
  #filter out all zeros
  filter(orig_value > 0) %>% 
  #group by var, find min value (across all months and years)
  group_by(variable) %>% 
  summarize(min_origval = min(orig_value, na.rm = TRUE),
            half_min_origval = min_origval / 2) %>% #,
  #except for each time replace zero...should put random number in for THAT
  #so...bump that + runif(1) step to mutate ifelse statement
  #half_min_plus_rand = half_min_val + runif(1)) %>% 
  as.data.frame() %>% 
  select(-min_origval)

# filter(fixed_na_dat, variable == "picoplankton_density") %>% summary()
# filter(fixed_na_dat, variable == "picoplankton_density") %>% select(orig_value) %>% unique() %>% min(na.rm = TRUE)
# min picoplankton is still pretty high and that's a legit value, FYI

#so for some, there won't be any zeros, so these min vals won't get used
#that's fine, just keep in mind

dat_long_minval <- merge(fixed_na_dat, min_origvals, by = "variable") %>% arrange(year, month)

# head(dat_long_minval, 20)
# filter(dat_long_minval, variable == "Epischura_adult_count")
# #gah, need to put zoop back in...didn't have NAs to replace, but need to go back and insert back in
# filter(dat_long_minval, variable == "unid_nano_density")

#each group should get its own random value
#use runif to get half min orig val + random value between 0 and 1
#set seed so same random set each time
set.seed(1)
#for each variable, get a random number between 0 and 1 
#to be random add for gorup half min vals to replace zeros
random_df <- data.frame(unique(dat_long_minval$variable), runif(n_distinct(dat_long_minval$variable)))
names(random_df) <- c("variable", "random_val")

#merge and replace zeros
dat_replace_zeros <- merge(dat_long_minval, random_df, by = "variable", all = TRUE) %>% 
  #for 0 values (after fixing NAs), replace with half min orig val + random val
  mutate(fixed_zero_value = ifelse(fixed_na_value == 0, (half_min_origval + random_val), fixed_na_value)) 

#in this method, within a variable, using same random number
#therefore, within each variable, all zero values with now be the SAME non-zero value
#(group min val / 2 + group random number)

# ----> IMPORTANT INFORMATION ABOUT NA AND ZERO VALUES

#so ok - NA is replaced appropriately with a zero, which is then replaced appropriately with a non-zero
#somewhat...squishy, but all legit
#just be sure to be clear about workflow process - NAs replaced first, then zeros
#some of those are NOT TRUE ZEROS but it's how the workflow is handling them

# filter(dat_replace_zeros, orig_value == 0) %>% arrange(year, month) %>% head(10)
# filter(dat_replace_zeros, orig_value == 0 & variable == "Achnanthes_density") %>% head()
# filter(dat_replace_zeros, is.na(orig_value)) %>% head(20)
# filter(dat_replace_zeros, year == 1991 & variable == "Achnanthes_density")
# #at some point NAs became zeros?
# filter(fixed_na_dat, year == 1991 & variable == "Achnanthes_density")
# filter(avg_phyto, month_mean_month == 6)
# #ok, well the average for Achnanthes in month 6 is indeed zero...
# filter(dat_phyto, month == 6) %>% select(Achnanthes_density)


#======================================================#
##=============== back to MARSS format ===============##
#======================================================#

#have replaced NA values and zero values
head(dat_replace_zeros)

dat_full_wide <- dat_replace_zeros %>% 
  select(year, month, variable, fixed_zero_value) %>% 
  dcast(year + month ~ variable, value.var = "fixed_zero_value")

names(dat_full_wide) <- paste("fixed", names(dat_full_wide), sep = "_")

dat_full_wide <- rename(dat_full_wide, year = fixed_year, month = fixed_month)

head(dat_full_wide)
summary(dat_full_wide)

#rename some
dat_full_wide <- dat_full_wide %>% 
  rename(fixed_chla_avg = fixed_chla_050depth_avg,
         fixed_chla_sum = fixed_chla_050depth_sum)


# #what??
# filter(dat_full_wide, is.na(fixed_Achnanthes_density))
# filter(dat_replace_zeros, year == 1981 & month == 1 & variable == "Achnanthes_density")
# #how are there still NA values?
# #oh wait yeah - there are NAs when it was only one NA
# #ok, that's fine

write.csv(dat_full_wide, "../data/fixed_mar_dat.csv", row.names = FALSE)
