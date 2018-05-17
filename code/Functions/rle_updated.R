# ----> updated rle fxn that counts runs of NA

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