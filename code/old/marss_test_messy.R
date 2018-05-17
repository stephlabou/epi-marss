###########################################
############ MARSS analysis  ##############
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

#load packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(MAR1)
library(MARSS)

#read in data
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

summary(dat)

#most of phytos have 115 NAs
phyto_na <- dat %>% 
            filter(is.na(Achnanthes_density))

summary(phyto_na) #115 NAs
str(phyto_na) #115

head(phyto_na)

phyto_na %>% select(year, month) 
#so...no phyto data for most months in 2000s
#hmmm...that's quite a few missing months in a row
#so really, we would want to stop in mid-2008
#and 1991 is just a wash, same with 2000 and 2001
#really, the whole 2000s are spotty
#go back and check what's up with those months/years in the main data

#well well well...zoop data only goes through Dec 2003
#so 12/2003 will be our upper limit for dates
#whoa who whoa phyto only goes through 12/1999
#so we won't even GET to the 2000s, stop at Dec 1999

#ok then, limit to 12/1999 (or year < 2000) and see what that looks like

dat_pre2000 <- dat %>% filter(year < 2000)

summary(dat_pre2000)

#ok, still a few missing Phyto and a LOT of missing chla, what the heck

misschla <- filter(dat_pre2000, is.na(chla_050depth_avg))
misschla %>% select(year, month) %>% unique()
#ugh, gross, lots are missing
#well darn it
#ok, so chla doesn't start until 1979
#totally jumps over 1984, 1987 is mostly missing
#at least the 1990s are pretty solid

#ok, so now we want to start in 1979 and go until end of 1999
#so...20 years

#I'm sure there's a way to plot this out buuuuut 
#for now just eyeballing and getting a sense of where missing values are

# chla_monthly %>% group_by(year) %>% summarize(n_mon = length(month)) %>% as.data.frame()
# #1985-1989 aren't great - other years are fine
# 
# temp_monthly %>% group_by(year) %>% summarize(n_mon = length(month)) %>% as.data.frame()
# #temp looks really good
# 
# zoop_monthly %>% group_by(year) %>% summarize(n_mon = n_distinct(month)) %>% as.data.frame()
# #zoop looks great
# 
# phy_monthly %>% group_by(year) %>% summarize(n_mon = n_distinct(month)) %>% as.data.frame()
# #phyto looks great

#so chla is the weak link
#decisions, decisions...

# SH:
# Where missing more than 2 in a row - in the past I used averages for the whole time period for that month. 
# If it's missing more than 2 in a row for any variable, the MAR is capable of just skipping 
# those months that have missing values.

#ok, so try this approach

#make sure every year represented
yrs <- seq(1979, 2008, 1)
mons <- seq(1, 12, 1)

full_yr_mon <- expand.grid(yrs, mons) %>% rename(year = Var1, month = Var2)

dat_time <- merge(dat, full_yr_mon, by = c("year", "month"), all = TRUE)

#this results in lots of NAs, so first need to limit to years of interest 
#(from previous explore)

dat_79_99 <- dat_time %>% 
                 filter(year >= 1979 & year <= 1999)

summary(dat_79_99)

chla_monthly_avg <- dat_79_99 %>% 
                    select(year, month, chla_050depth_avg, chla_050depth_sum) %>% 
                    #find monthly average across this 20 year time span
                    group_by(month) %>% 
                    summarize(chla_050depth_avg_month = mean(chla_050depth_avg, na.rm = TRUE),
                              chla_050depth_sum_month = mean(chla_050depth_sum, na.rm = TRUE)) %>% 
                    as.data.frame()

#note: doesn't check whether missing 2 in a row...
dat_79_99_fillna <- dat_79_99 %>% 
                        #merge monthly chla values (over entire 1979-1999)
                        merge(chla_monthly_avg, by = "month", all.x = TRUE) %>% 
                        arrange(year, month) %>% 
                        #replace chla NA values
                        #note: does for *all* NAs, not just 2 skips and more
                        mutate(chla_avg = ifelse(is.na(chla_050depth_avg), chla_050depth_avg_month, chla_050depth_avg),
                               chla_sum = ifelse(is.na(chla_050depth_sum), chla_050depth_sum_month, chla_050depth_sum)) %>% 
                        #add note column
                        mutate(notes = ifelse(is.na(chla_050depth_avg) | is.na(chla_050depth_sum), 
                                              "chla temporal monthly", NA)) %>% 
                        #remove old columns
                        select(-chla_050depth_avg, -chla_050depth_sum, 
                               -chla_050depth_avg_month, -chla_050depth_sum_month)

#can always go back and adjust based on if 2 or more skips then replace
   
summary(dat_79_99_fillna)    

#need to deal next with missing temp, missing phytos

#ok, what if I wanted to count number of missing?

chla_skips <- dat_79_99 %>% 
              summarize(test = rle(chla_050depth_avg))

rle(dat_79_99$chla_050depth_avg) %>% head()
rle(is.na(dat_79_99$chla_050depth_avg)) %>% head()

?lag

#lead = next
#lag = previous
#could maybe use this?

test <- c(1.0, 5.3, NA, NA, 2.7, NA, 4.1, NA, NA, 1.8)
test2 <- c(2, 3, 4, 5, 6, 7, 8, 9, 5, 3)

df <- data.frame(test, test2)
df <- rename(df, col1 = test, col2 = test2)

# NA == NA
# match(NA, NA)
# match(1, 2)

df %>% 
  mutate(cat = ifelse(((identical(lag(col1), col1)) == TRUE), "yes", "no"))

#hmmm, keeps resulting in no...

df %>% 
  mutate(cat = ifelse(identical(lag(col1), col1) == TRUE, "yes", "no"))


a <- rle(is.na(df$col1))
a
a$lengths
a$values

#hmmm get T/F and convert them sum?

df %>% mutate(idea = ifelse(is.na(col1), 1, 0))

df %>% mutate(idea = ifelse(is.na(col1), 1, 0)) %>% group_by(col1) %>% mutate(runs = sum(idea))
#no, that's not right...

df %>% mutate(idea = ifelse(is.na(col1), 1, 0)) %>% group_by(col1, idea) %>% mutate(runs = sum(idea))

#grrr....

rle(df$col1)$lengths
rle(df$col1)$values

df2 <- df %>% mutate(col3 = ifelse(is.na(col1), -99999, col1))
rle(df2$col3)$lengths
rle(df2$col3)$values
#ok, but the output isn't same length, ugh

sequence(rle(df2$col3)$lengths)
df2$col3
#not quite...
sequence(rle(df2$col3)$values)


#From: https://stat.ethz.ch/pipermail/r-help/2011-June/281730.html
#updated function that deals with Na values
rle_updated <-function (x)
{
  if (!is.vector(x) && !is.list(x))
    stop("'x' must be an atomic vector")
  n <- length(x)
  if (n == 0L)
    return(structure(list(lengths = integer(), values = x),
                     class = "rle"))
  
  #### BEGIN NEW SECTION PART 1 ####
  naRepFlag<-F
  if(any(is.na(x))){
    naRepFlag<-T
    IS_LOGIC<-ifelse(typeof(x)=="logical",T,F)
    
    if(typeof(x)=="logical"){
      x<-as.integer(x)
      naMaskVal<-2
    }else if(typeof(x)=="character"){
      naMaskVal<-paste(sample(c(letters,LETTERS,0:9),32,replace=T),collapse="")
    }else{
      naMaskVal<-max(0,abs(x[!is.infinite(x)]),na.rm=T)+1
    }
    
    x[which(is.na(x))]<-naMaskVal
  }
  #### END NEW SECTION PART 1 ####
  
  y <- x[-1L] != x[-n]
  i <- c(which(y), n)
  
  #### BEGIN NEW SECTION PART 2 ####
  if(naRepFlag)
    x[which(x==naMaskVal)]<-NA
  
  if(IS_LOGIC)
    x<-as.logical(x)
  #### END NEW SECTION PART 2 ####
  
  structure(list(lengths = diff(c(0L, i)), values = x[i]),
            class = "rle")
}


df

lens <- rle_updated(df$col1)$lengths
vals <- rle_updated(df$col1)$values

df3 <data.frame(rle_updated(df$col1)$lengths)

rle_updated(df$col1) 

do.call(rbind.data.frame, rle_updated(df$col1) )

df3 <- df

df3$lengths <- rle_updated(df$col1)$lengths

#ok, but need to expand if >1...
lens
vals

size <- rle_updated(df$col1)$lengths %>% length()

df_test <- data.frame(nrow=size, ncol=2)

df_test

test1 <- data.frame(lens = integer(8),
                    vals = numeric(8))

test1
  
test1$lens <- rle_updated(df$col1)$lengths
test1$vals <- rle_updated(df$col1)$values

test1

#need to repeat rows with lengths >1 for number of times value of length is

test1[rep(seq_len(nrow(test1)), each=2),]


test1[rep(seq_len(nrow(filter(test1, lengths))), each=2),]

#hmm another test from stack overflow...
df12 <- transform(df,id=as.numeric(factor(col1)))
df12
#nope, NAs remains as NA...

?rep

??expandRows

# try lapply([i], rep(val, val times))
#basically try lapply

test1

lapply(test1$lens, rep_len(lens, length.out = lens))

for (i in 1:nrow(test1)) {
  
  dat <- test1[i,]
  
  new <- rep_len(dat$lens, length.out = dat$lens)
  
  if (i == 1) {
    new_lens <- new
  }
  
  if (i > 1) {
    new_lens <- c(new_lens, new)
  }
}

#yyyyaaaaaayyyyy it looks like this works
#need to test more but yeah sure

test1
df
new_lens

df_out <- cbind(df, new_lens)
df_out

#yes - test again with slightly larger df with differing NA runs
#to test out before run on whole

#######################################################################################

a <- c(3.56, NA, 23.4, 98.0, NA, NA, NA, 3.4, 3.4, 98, NA, NA)
df <- data.frame(a)

#"original" data frame
df_orig <- df

#updated rle fxn that counts runs of NA
rle_updated <-function (x)
{
  if (!is.vector(x) && !is.list(x))
    stop("'x' must be an atomic vector")
  n <- length(x)
  if (n == 0L)
    return(structure(list(lengths = integer(), values = x),
                     class = "rle"))
  
  #### BEGIN NEW SECTION PART 1 ####
  naRepFlag<-F
  if(any(is.na(x))){
    naRepFlag<-T
    IS_LOGIC<-ifelse(typeof(x)=="logical",T,F)
    
    if(typeof(x)=="logical"){
      x<-as.integer(x)
      naMaskVal<-2
    }else if(typeof(x)=="character"){
      naMaskVal<-paste(sample(c(letters,LETTERS,0:9),32,replace=T),collapse="")
    }else{
      naMaskVal<-max(0,abs(x[!is.infinite(x)]),na.rm=T)+1
    }
    
    x[which(is.na(x))]<-naMaskVal
  }
  #### END NEW SECTION PART 1 ####
  
  y <- x[-1L] != x[-n]
  i <- c(which(y), n)
  
  #### BEGIN NEW SECTION PART 2 ####
  if(naRepFlag)
    x[which(x==naMaskVal)]<-NA
  
  if(IS_LOGIC)
    x<-as.logical(x)
  #### END NEW SECTION PART 2 ####
  
  structure(list(lengths = diff(c(0L, i)), values = x[i]),
            class = "rle")
}

#use rle_updated on orig df
lens <- rle_updated(df_orig$a)$lengths
vals <- rle_updated(df_orig$a)$values

#make df of rel_updated output
df_rle <- data.frame(lens = integer(length(lens)),
                     vals = numeric(length(vals)))
df_rle$lens <- lens
df_rle$vals <- vals

df_rle

#use for loop to expand when runs multiple
for (i in 1:nrow(df_rle)) {
  
  dat <- df_rle[i,]
  
  new <- rep_len(dat$lens, length.out = dat$lens)
  
  if (i == 1) {
    new_lens <- new
  }
  
  if (i > 1) {
    new_lens <- c(new_lens, new)
  }
}

new_lens

#bind back to orig df
df_updated <- cbind(df_orig, new_lens)

df_updated

#yep, that definitely looks like it's working properly

#now..can I turn this WHOLE THING into an easy function...
source("Functions/rle_updated.R")
source("Functions/tag_runs_incl_na.R")

df

wrangle_na_runs(df, "a")

wrangle_na_runs(df, df$a)

b <- c(3.5, 2.3, NA, NA, 7.3, 8.1, 2.3, NA, 1.2, NA, NA, NA)

test <- cbind(df, b)

test

lapply(test$a, function (x) wrangle_na_runs_apply(x))

test
       
wrangle_na_runs_apply(test, test$a)    

test <- test %>% 
          mutate(calc = wrangle_na_runs_col(test$a))

wrangle_na_runs_col <- function(x) {
  
  if (!is.vector(x))
    stop("'x' should be a dataframe column in form df$x")
  
  #col1 <- paste("'", col1, "'", sep = "")
  
  #vec <- df_orig[, a]
  
  vec <- test$a
  
  #use rle_updated on orig df
  lens <- rle_updated(vec)$lengths
  vals <- rle_updated(vec)$values
  
  #make df of rel_updated output
  df_rle <- data.frame(lens = integer(length(lens)),
                       vals = numeric(length(vals)))
  df_rle$lens <- lens
  df_rle$vals <- vals
  
  #use for loop to expand when runs multiple
  for (i in 1:nrow(df_rle)) {
    
    dat <- df_rle[i,]
    
    new <- rep_len(dat$lens, length.out = dat$lens)
    
    if (i == 1) {
      new_lens <- new
    }
    
    if (i > 1) {
      new_lens <- c(new_lens, new)
    }
  }
  
}
       
out <- wrangle_na_runs_col(test$a)    
       
test %>% mutate(out = wrangle_na_runs_col(test$a))    

###########################

chla <- c(1.2, 3.6, NA, 9.2, NA, NA, 3.4, 5.8, NA, NA, NA, 4.3, 4.3, 1.9)
tmp <- c(NA, NA, 3.5, 6.2, 3.4, NA, 9.8, 9.8, NA, 3.4, NA, NA, NA, NA)

dat_env <- data.frame(chla, tmp)

dat_env

source("Functions/tag_runs_incl_na.R")

dat_env_updated <- dat_env %>% 
                    mutate(chla_runs = wrangle_na_runs_col(dat_env$chla),
                           tmp_runs = wrangle_na_runs_col(dat_env$tmp))

dat_env_updated <- dat_env_updated %>% 
                    select(chla, chla_runs, tmp, tmp_runs)

dat_env_updated

#does this work with grouping?
#wait, does it even matter?
#as long as arranged chronologically, should be fine




       
       