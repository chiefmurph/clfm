split.ata.periods <- function(charvec) {
  # 4/23/10 Hmisc beging taken out of commission
  # I think it was only used for sedit. Replaced by gsub
  #require(Hmisc)
  # sedit requires Hmisc package
  
  # Determine the ages in a character vector of development periods 
  #   that look like "12-24"   "24-36"   "36-48" "48-Ult"  or
  #   "12 to 24"   "24 to 36"   "36 to 48"  "48 to Ult."
  # The characters '-' and 'to' are in the list splitlist.
  # Return a numeric matrix with the first row being the beginning ages
  #   and the second row being the ending ages.
  # If "Ult" or "Ult." or "Ultimate" appear in charvec, convert to Inf
  
  # Note: Ult. must precede Ult in following character vector, otherwise
  #       it will become Inf. which is non-numeric
  # 4/23/10 change from vector for sedit to string for gsub
  #ult.synonyms <- c("INF","ULTIMATE","ULT.","ULT","INFINITY")
  ult.synonyms <- "INF|ULTIMATE|ULT.|ULT|INFINITY"
  splitlist <- c("-","TO")
  
  L<-length(charvec)
  # convert to upper case
  period<-toupper(charvec) 
  #If any "-"'s, convert to "TO"'s
  # period <- sedit(period,"-","TO") 4/23/10 Hmisc being taken out of commission
  period <- gsub("-","TO",period) 
  # Now split the names, trim away blanks.
  agesc<-trim(unlist(strsplit(period,"TO")))
  # Two cases: original names only had one name, or had more than one name
  if (L==1) { # only one name
    beginage <- suppressWarnings(as.numeric(agesc[1]))
    if (length(agesc)>2) stop("Too many separators in tail name")
    if (length(agesc)>1) { # There was a separator
      # age after separator must be synonomous with Inf
      if (!(agesc[2] %in% ult.synonyms))
        stop("Tail age not to Ultimate")
    }
    endage <- Inf
    return(matrix(c(beginage,endage),ncol=2))
  }
  # Three cases:
  #   split length=L -> no separators
  #   split length=2L -> each name had a separator
  #   split length=2*(L-1)+1 -> the L-1 names other than tail had a 
  #             separator, tail name must be synonymous with "ultimate"
  agescL<-length(agesc)
  if (agescL==L) { # names = just a vector of ages.
    # Two cases:
    #   Last name is blank or "Tail"
    #   Last name is an age
    if (agesc[L]=="" | agesc[L]=="TAIL") {        
      ages <- suppressWarnings(as.numeric(agesc[-L]))
      if (NA %in% ages) stop("Non-numeric names for ages")
      # guess age of tail
      dif <- diff(ages)
      if (!all(diff(dif)==0)) 
        stop("Non-constant age differences, cannot guess tail age")
      beginage <- c(ages,ages[length(ages)]+dif[1])
      endage <- c(beginage[-1],Inf)
    }
    else {
      agesc <- troUltToInf(agesc)# 10/13/09 "ULT" can be a valid age
      ages <- suppressWarnings(as.numeric(agesc))
      if (NA %in% ages) stop("Non-numeric names for ages:",agesc)
      beginage <- ages
      endage <- c(ages[-1],Inf)
    }
  }
  else
    if (agescL==2*L) { # separators in all names
      # Convert all "ultimates" to infinity
      # 4/23/10: Hmisc going out of commission
      # agesc <- sedit(agesc,ult.synonyms,"Inf")
      agesc <- gsub(ult.synonyms, "Inf", agesc, ignore.case=TRUE)
      ages <- suppressWarnings(as.numeric(agesc))
      if (NA %in% ages) stop("Non-numeric names for ages")
      dim(ages) <- c(2,L)
      beginage<-ages[1,]
      endage<-ages[2,]
    }
  else
    if (agescL==(2*(L-1)+1)) {
      # last name must be "TAIL" or a synonym for "Ultimate"
      if (toupper(agesc[agescL]) != "TAIL" & 
          !(agesc[agescL] %in% ult.synonyms))
        stop("Cannot interpret tail name '",charvec[L],"'")
      agesc <- agesc[-agescL]
      ages <- suppressWarnings(as.numeric(agesc))
      if (NA %in% ages) stop("Non-numeric names for ages")
      dim(ages) <- c(2,L-1)
      beginage<-c(ages[1,],ages[2,ncol(ages)])
      endage <- c(ages[2,],Inf)
    }
  else stop("Cannot interpret ata names as ages")
  
  # A couple of final checks
  # Begin ages must be increasing.
  if (!is.agevector(beginage)) stop("Ages not strictly increasing")
  # Do begin ages_k = end ages_k-1?
  if (L>1) {
    for (k in 2:L) if (beginage[k] != endage[k-1]) 
      stop("Begin/end age mismatch")
  }
  matrix(c(beginage,endage),ncol=2)
}

trim<-function(x) sub("[[:space:]]+$","",sub("^[[:space:]]+","",x))

troUltToInf <- function(x) {
  # Wherever "Ult" or "Ultimate" occurs in character vector x,
  # replace it with "Inf".
  x[which(toupper(x) %in% c("ULT","ULTIMATE"))] <- "Inf"
  return(x)
}

is.agevector <- function(x,what=c("test","convert")) {
  # x must be convertible to numeric vector, positive, strictly increasing
  # If "test", return T/F.
  # If "convert", return the numeric vector.
  # To test result and "convert", use 
  #   if (is.logical(a<-is.agevector(v,"convert"))) #failed test
  what <- match.arg(what)
  x<-troUltToInf(x)
  x<-suppressWarnings(as.numeric(x))
  if (any(is.na(x))) return(FALSE)
  if (any(x<0)) return(FALSE)
  if (length(x)>1) 
    if (any(x[-1]<=x[-length(x)])) return(FALSE)
  if (what=="test") return(TRUE)
  return(x)
}

VWata.tda <- function(tda) {
  d<-tda[!is.na(tda$ata),] # ata=NA only for unmatched begin, end values
  tapply(d$value.next,d$age,sum)/tapply(d$value,d$age,sum)
}
SAata.tda <- function(tda) tapply(tda$ata,tda$age,mean,na.rm=TRUE)
RAata.tda <- function(tda) {
  d<-tda[!is.na(tda$ata),] # ata=NA only for unmatched begin, end values
  tapply(d$value.next*d$value,d$age,sum)/tapply(d$value^2,d$age,sum)
}

approxeq <- function (x, y, tol=.Machine$double.eps^0.5) abs(x-y)<tol

approx.equal <- function (x, y, tol=.Machine$double.eps^0.5) abs(x-y)<tol


