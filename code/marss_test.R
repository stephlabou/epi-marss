###########################################
############ MARSS analysis  ##############
###########################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(MAR1)
library(MARSS)

dat <- read.csv("../data/mar_dat.csv")

head(dat)
summary(dat)


# ----> how to deal with NA values?
#from SH:
# In the past I have interpolated when it is just one month missing, and not two in a row. 
# Where missing more than 2 in a row - in the past I used averages for the whole time period for that month. 
# If it's missing more than 2 in a row for any variable, the MAR is capable of just skipping 
# those months that have missing values. But we want to know when it is happening.

# ----> handling real zeros
# For real zeroes, MAR also will blow up on zeroes if there are too many in a row... 
# a few is not a big deal, as I recall. 
# So the rule of thumb is to calculate half the lowest observable (or observed non-zero) value, 
# plus or minus some random number, to replace zeroes.  For now, I would say leave real zeroes in 
# and we can always transform them at the time of analysis.









