# ----> function to wrangle runs of numbers, incl. NA runs

# x <- SET COLUMN OF INTEREST (df$x)

wrangle_na_runs_col <- function(x) {
  
  if (!is.vector(x))
    stop("'x' should be a dataframe column in form df$x")
  
  vec <- x
  
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
  #this is NOT best practice
  #but I am tired of fighting this
  runs_vec <<- new_lens
  
}


# # ----> function to wrangle runs of numbers, incl. NA runs
# 
# # df_orig <- SET DF OF INTEREST
# # x <- SET COLUMN OF INTEREST (df$x)
# 
# wrangle_na_runs <- function(df_orig, x) {
#   
#   if (!is.vector(x))
#     stop("'x' should be a dataframe column in form df$x")
#   
#   #col1 <- paste("'", col1, "'", sep = "")
#   
#   #vec <- df_orig[, a]
#   
#   vec <- x
#   
#   #use rle_updated on orig df
#   lens <- rle_updated(vec)$lengths
#   vals <- rle_updated(vec)$values
#   
#   #make df of rel_updated output
#   df_rle <- data.frame(lens = integer(length(lens)),
#                        vals = numeric(length(vals)))
#   df_rle$lens <- lens
#   df_rle$vals <- vals
#   
#   #use for loop to expand when runs multiple
#   for (i in 1:nrow(df_rle)) {
#     
#     dat <- df_rle[i,]
#     
#     new <- rep_len(dat$lens, length.out = dat$lens)
#     
#     if (i == 1) {
#       new_lens <- new
#     }
#     
#     if (i > 1) {
#       new_lens <- c(new_lens, new)
#     }
#   }
#   
#   new_lens
#   
#   #bind back to orig df
#   df_updated <- cbind(df_orig, new_lens)
#   
#   #remove unnecessary stuff from environment
#   rm(lens, vals, df_rle, dat, new, new_lens)
#   
#   attach(df_updated)
#   
# }
# 
# # ----> function to wrangle runs of numbers, incl. NA runs
# 
# # df_orig <- SET DF OF INTEREST
# # x <- SET COLUMN OF INTEREST (df$x)
# 
# wrangle_na_runs_apply <- function(df_orig, x) {
#   
#   if (!is.vector(x))
#     stop("'x' should be a dataframe column in form df$x")
#   
#   #col1 <- paste("'", col1, "'", sep = "")
#   
#   #vec <- df_orig[, a]
#   
#   vec <- x
#   
#   #use rle_updated on orig df
#   lens <- rle_updated(vec)$lengths
#   vals <- rle_updated(vec)$values
#   
#   #make df of rel_updated output
#   df_rle <- data.frame(lens = integer(length(lens)),
#                        vals = numeric(length(vals)))
#   df_rle$lens <- lens
#   df_rle$vals <- vals
#   
#   #use for loop to expand when runs multiple
#   for (i in 1:nrow(df_rle)) {
#     
#     dat <- df_rle[i,]
#     
#     new <- rep_len(dat$lens, length.out = dat$lens)
#     
#     if (i == 1) {
#       new_lens <- new
#     }
#     
#     if (i > 1) {
#       new_lens <- c(new_lens, new)
#     }
#   }
#   
#   new_lens
#   
#   #bind back to orig df
#   df_orig <- cbind(df_orig, new_lens)
#   
#   #remove unnecessary stuff from environment
#   #rm(lens, vals, df_rle, dat, new, new_lens)
#   
# }
