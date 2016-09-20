###################################################
######### Utilities
###################################################
mvexpo <- function(M,v) {
  # exponentiate all elements in column of matrix M by corresponding
  # elements of vector v. 
  # V can be a scalar
  # If length(V) > ncols(M), R gives a warning, can get unexpected results.
  # If length(v)< ncols(M), elements of v are recycled.
  # If v is a matrix, must have same shape as M or R gags 'non-conformable'
  M<-as.matrix(M)
  ifelse(is.na(M),NA,t(t(M)^v))
  }

mvmult <- function(M,v) t(t(M)*v)
   # returns matrix where columns of matrix M are multiplied
   # by corresponding elements of vector v
flip.matrix <- function(M) M[nrow(M):1,]
vector.all.equal <- function(x)
  isTRUE(all.equal(x, rep(x[1], length(x)),check.attributes=FALSE))
triCumulate <- function(M) t(apply(M,1,cumsum))
triDecumulate <- function(M) {
  dm <- cbind(M[,1],t(apply(M,1,diff)))
  colnames(dm) <- colnames(M)
  return(dm)
  }

quotevector<-function(v) { # quotes both sides, one long string
  if (length(v) > 1)
    paste("'",
          paste(head(c(v),-1),"' '",collapse="",sep=""),
          tail(c(v),1),
          "'",sep="")
  else paste("'",v,"'",sep="")
  }

prefix.quote<-function(v) paste("'",v,sep="")

# 4/23/10 Hmisc going out of commission
all.is.numeric<-function (x, what = c("test", "vector"), extras = c(".", "NA")){
    what <- match.arg(what)
    old <- options(warn = -1)
    on.exit(options(old))
    x <- sub("[[:space:]]+$", "", x)
    x <- sub("^[[:space:]]+", "", x)
    xs <- x[x %nin% c("", extras)]
    isnum <- !any(is.na(as.numeric(xs)))
    if (what == "test") 
        isnum
    else if (isnum) 
        as.numeric(x)
    else x
    }
`%nin%`<-function (x, table) match(x, table, nomatch = 0) == 0

vector.order <- function(v)
  if (is.logical(v<-all.is.numeric(v,"vector"))) order(v) else order(v)

is.Ultimate <- function(x)
  # x is a character vector of ages 
  x=="Inf"
  
is.Finite <- function(x)
  # x is a character vector of ages 
  x!="Inf"

is.Infinite <- function(x)
  # x is a character vector of ages 
  x=="Inf"
  
troUltToInf <- function(x) {
  # Wherever "Ult" or "Ultimate" occurs in character vector x,
  # replace it with "Inf".
  x[which(toupper(x) %in% c("ULT","ULTIMATE"))] <- "Inf"
  return(x)
  }

nexttdfage <- function(tdf,age) {
  ages <- attr(tdf,"ages")
  return(ages[match(age,ages)+1])
  }
  
prevtdfage <- function(tdf,age) {
  ages <- attr(tdf,"ages")
  return(c(NA,ages)[match(age,ages)])
  }
  
ata.to.mat <- function(tda,what="ata") {
  ages <- attr(tda,"ages")
  tda<-tda[!is.na(tda$value.next),]
  M <- tdf.as.tri(tda,what=what)
  names(M)<-paste(ages[-length(ages)],ages[-1],sep="-")
  return(M)
  }

ata.to.matOld <- function(tda,what="ata") {
  ages <- attr(tda,"ages")
  tda<-tda[!is.na(tda$beginvalue),]
  M <- tdf.as.tri(tda,what=what)
  names(M)<-paste(ages[-length(ages)],ages[-1],sep="-")
  return(M)
  }

reshapeATADFtoMat <- function(tda,what="ata",colid="dev.period") {
  tda<-tda[!is.na(tda$value.prev),]
  M <- tdf.as.tri(tda,what=what,colid=colid)
  if (colid=="dev.period") 
    names(M) <- aggregate(tda$dev.period,by=list(tda$age),unique)[,2]
  return(M)
  }

append.list <- function(list1,list2) {
  # concatenating lists preserves names only. must copy over other attributes
  list3 <- c(list1,list2)
  attributes(list3) <- c(attributes(list3),attr.x.names(list1))
  return(list3)
  }
  
attributes.excluding <- function(obj,excl=NULL)
  attributes(obj)[-match(excl,names(attributes(obj)))]
  
append.df <- function(df1,df2) {
  # Append columns of df2 to right side of df1.
  # Similar to cbind but preserves attributes of df1.
  df <- df1
  for (nam in names(df2)) df[[nam]] <- df2[[nam]]
  return(df)
  }
  
cbind.df <- function(df1,df2) {
  # Append columns of df2 to right side of df1.
  # Similar to cbind but preserves attributes of df1.
  df1[(length(df1)+1):(length(df1)+length(df2))]<-df2
  df1
  }
  
tdfSumByAge <- function(tdf,what="value",sumname="Sum") {
  # return a 1xn matrix of sums by age of data frame tdf column named what
  what <- match.arg(what,names(tdf))
  if (missing(sumname)) sumname <- paste("Sum(",what,") by age",sep="")
#  m<-t(aggregate(tdf[what],by=list(as.numeric(tdf$age)),sum)[what])
#  m<-t(aggregate(tdf[what],by=list(as.numeric(tdf$age)),sum))
  m<-aggregate(tdf[what],by=list(age=tdf$age),sum)
  m<-m[order(m$age),]
  age<-m$age
  m<-m[[what]]
  dim(m)<-c(1,length(m))
  rownames(m)<-sumname
  colnames(m)<-age
  return(m)
  }

tdfAvgByAge <- function(tdf,what="value",sumname="Avg",...) {
  # return a vector of averages by age of data frame tdf column named what
  what <- match.arg(what,names(tdf))
  if (missing(sumname)) sumname <- paste("Avg(",what,") by age",sep="")
#  m<-t(aggregate(tdf[what],by=list(as.numeric(tdf$age)),sum)[what])
#  m<-t(aggregate(tdf[what],by=list(as.numeric(tdf$age)),sum))
  m<-aggregate(tdf[what],by=list(age=tdf$age),mean,...)
  m<-m[order(m$age),]
  age<-m$age
  m<-m[[what]]
#  dim(m)<-c(1,length(m))
#  rownames(m)<-sumname
#  colnames(m)<-age
  names(m)<-age
  attr(m,"Name")<-sumname
  return(m)
  }

trim<-function(x) sub("[[:space:]]+$","",sub("^[[:space:]]+","",x))
#gsub("([^-0-9])"," ",c("a5","ab5","FSHUIE7HYB EFS90--5cFESBUFEHSU6    TO 7 GSBKHRSFE"))
#gsub("-+","-","-----5")
#gsub("-+","-",gsub("([^-0-9])"," ",c("a5","ab5","FSHUIE7HYB EFS90--5cFESBUFEHSU6    TO 7 GSBKHRSFE")))
#x<-"7. --5  -- 8  TO 7 GSInfRSFE"
#y<-gsub("-+","-",x)
#strsplit(y,"([^-Inf0-9])")

# following is commented out, can be found in agevector.r
#regexpRun <- function(needle="([-+]?[0-9]*\\.?[0-9]+|Infinity|Inf\\.|Inf|Ultimate|Ult\\.|Ult)", 
#                      haystack, ignore.case=TRUE, ...) {
#	myRegexp<-gregexpr(needle, haystack, ignore.case=ignore.case, ...)
#	result<-character(length(myRegexp[[1]]))
#	for (i in 1:length(myRegexp[[1]])) {
#		pos<-myRegexp[[1]][i]
#		len<-attr(myRegexp[[1]],"match.length")[i]
#		tot<-pos+len-1
#		result[i]<-substr(haystack,pos,tot)
#	}
#	result
#}
#ultToInf<-function(x, ignore.case=TRUE, ...) {
#	result<-character(length(x))
#	for (i in 1:length(x)) {
#		result[i]<-gsub("Infinity|Inf\\.|Ultimate|Ult\\.","Inf",x[i] , ignore.case=ignore.case, ...)
#	}
#	result
#}



c.factor = function(x,y) { # concatenate two factors
  newlevels = union(levels(x),levels(y)) 
  m = match(levels(y), newlevels) 
  ans = c(unclass(x),m[unclass(y)]) 
  levels(ans) = newlevels 
  class(ans) = "factor" 
  ans 
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

#
##############################################################
# fata
##############################################################
#
CBE <- function(C,zero.future=FALSE) {
  #   Given a "triangle" C stored as a matrix, calculate
  #      losses at beginning and end of all development periods.

  #
  #   Return list: $begin and $end will be matrices of losses 
  #      as of the beginning and end of the development period,
  #      respectively.
  #
  #   Simply drop final column of C for $begin, first col for $end.
  #      Then, wherever the values as of the beginning or end
  #      of the period are missing, record $begin and $end as
  #      "missing" also because we need pairs of observations
  #      to compute development statistics.
  #
    Cbegin <- C[,-ncol(C)]
    Cend   <- C[,-1]
    missing.pairs <- is.na(Cbegin) | is.na(Cend)
    if (zero.future) missing.pairs <- missing.pairs | Cend==0
    Cbegin[missing.pairs] <- NA
    Cend[missing.pairs]   <- NA
    list(begin=Cbegin,end=Cend)
   }

Fata <- function(C) {
  # Calc triangle of actual link ratios matrix triangle C
  dev <- CBE(C) # Get losses @ begin, end of all development periods
  if (is.null(dev)) return(dev) # Something wrong with C
  Fmat <- dev$end/dev$begin
  colnames(Fmat) <- paste(colnames(dev$begin),colnames(dev$end),sep="-")
  # remove any rows that are all NA
  Fmat <- Fmat[apply(Fmat,1,function(x) !all(is.na(x))),]
  return(Fmat)
  }

LinkRatios <- function(X,nalastrow.rm=FALSE) {
    if (!is.matrix(X))
        if (is.null(
            X<-tryCatch(as.matrix(X,dimnames=dimnames(X)),
                error=function(...) {cat("Not a matrix\n");NULL})
                    )) return(NULL)
#    Y<-array(NA,dim(X)-c(0,1)) # NA's for cells with no future experience
    Y<-array(NA,dim=dim(X)-c(0,1),dimnames=list(dimnames(X)[[1]],head(dimnames(X)[[2]],-1)))
    # Find right edge of each row
    w<-lapply(1:nrow(X),function(i) tail(which(X[i,]!=0),1))
    # if row is all zeros, that entry in w list will have length 0
    for (i in 1:nrow(X)) {
        if (length(w[[i]])>0)
        if (w[[i]]>1) for (j in 1:(w[[i]]-1)) {
            z <-tryCatch(X[i,j+1]/X[i,j],error=function(...) {cat("Problem calculating link ratios from loss matrix\n");return(NULL)})
            if (is.null(z)) return(NULL)
            Y[i,j]<-z
            }
        }
    # In "usual" case, the most recent AY will be a row of all NS's, so
    #   remove it if so specified (default).
    if (nalastrow.rm) {
        w <- apply(Y,1,function(x) !all(is.na(x)))
        if (any(w)) Y<-Y[1:tail(which(w),1),] 
        else Y<-Y[FALSE,]
        }
    Y
    }

numata <- function(C) 
  # returns a vector which counts the number of non-NA link ratios 
  # available in each column of triangle C
  colSums(as.matrix(!is.na(Fata(C))))

Wata <- function(C,alpha=1) {
  # The link ratio weights per Bardis, Majidi, Murphy paper
  # Default of alpha=1 corresponds to volume weighted averages
  dev <- CBE(C) # Get losses @ begin, end of all development periods
  if (is.null(dev)) return(dev) # Something wrong with C
  # To get the weights, raise beginning period losses to power 2-alpha,
  # then divide by column sums, ignoring missing values.
  return(matrixcolpct(mvexpo(dev$begin,2-alpha)))
  }

LRfcn <- function(C,alpha=1,zero.future=FALSE) {
  # The Link Ratio function per Bardis, Majidi, Murphy FFM paper.
  # C can be a triangle, or just a pair of columns
  #   (if pair of columns, need "as.matrix" below so colSums works)
  # Default of alpha=1 corresponds to volume weighted averages
  # Returns vector of link ratios for that triangle and alpha.
  #
  # If alpha is a vector, return a matrix of link ratios where 
  # the link ratios are stored in columns, each element of alpha
  # corresponding to a different column of the matrix.
   
  # For some reason, NA^0=1 in R, but we want those to be NA's. 
  # Thus, the ifelse functions in denominators below.

  # Get losses @ begin, end of all development periods
  dev <- CBE(C,zero.future=zero.future)
  if (is.null(dev)) return(dev) # Something wrong with C; see CBE

  # Case where alpha is a scalar
  if (length(alpha)==1) {
    if (is.matrix(dev$begin) | is.data.frame(dev$begin)) {
      avgata <- colSums(mvexpo(dev$begin,(1-alpha))*dev$end,na.rm=TRUE) /
                colSums(mvexpo(dev$begin,(2-alpha))        ,na.rm=TRUE)
      }
    else # C was just a pair of columns, so begin and end are simply vectors
      avgata <- sum(mvexpo(dev$begin,(1-alpha))*dev$end,na.rm=TRUE) /
                sum(mvexpo(dev$begin,(2-alpha))        ,na.rm=TRUE)
    names(avgata) <- colnames(dev$beg)
    }
  # alpha not a scalar. Not programmed if array, so return NULL in that case
  # Case where alpha is a vector
  else
  if (is.matrix(alpha)) return (NULL)
  else { # alpha is a vector
    # If C has more than one development period, 
    # idea is that length(alpha)=ncols(C.begin), so the alpha corresponding
    # to a column of C comes from the appropriate element of alpha
    if (is.matrix(dev$begin) | is.data.frame(dev$begin))
      avgata <- colSums(mvexpo(dev$begin,(1-alpha))*dev$end,na.rm=TRUE) /
                colSums(mvexpo(dev$begin,(2-alpha))        ,na.rm=TRUE)  
      # C is a pair of columns, return LRfcn evaluated at each element of alpha
    else avgata <- sapply(alpha,FUN=function(x)
                  sum(mvexpo(dev$begin,(1-x))*dev$end,na.rm=TRUE)/
                  sum(mvexpo(dev$begin,(2-x) )       ,na.rm=TRUE))
    names(avgata) <- colnames(dev$beg)
    }
  # Case where alpha is a matrix. Not accomodated at this time.
  # alpha is an array of >= 3 dimensions. Not accomodated at this time.
  if (any(is.na(avgata))) {
    trowarning("LRfcn: NA's in ata set to 1.000: ",paste(sprintf("%.3f",avgata),collapes=""))
    avgata[is.na(avgata)] <- 1
    }
  avgata
  }

VWata <- function(C,zero.future=FALSE) LRfcn(C,1,zero.future=zero.future)
  # returns vector of volume weighted average link ratios

SAata <- function(C,zero.future=FALSE) LRfcn(C,2,zero.future=zero.future)
  # returns vector of simple average/straight average link ratios

RAata <- function(C,zero.future=FALSE) LRfcn(C,0,zero.future=zero.future)
  # returns vector of linear regression ("regression average") link ratios
  
VWata.tda <- function(tda) {
    d<-tda[!is.na(tda$ata),] # ata=NA only for unmatched begin, end values
    tapply(d$value.next,d$age,sum)/tapply(d$value,d$age,sum)
    }
SAata.tda <- function(tda) {
  d<-tda[!is.na(tda$ata),] # ata=NA only for unmatched begin, end values
  tapply(d$ata,d$age,mean,na.rm=TRUE)
}
RAata.tda <- function(tda) {
    d<-tda[!is.na(tda$ata),] # ata=NA only for unmatched begin, end values
    tapply(d$value.next*d$value,d$age,sum)/tapply(d$value^2,d$age,sum)
    }
    
  
info<-function(msg,show=F)
  if (show) cat(msg,"\n")
trowarning<-function(...,show.msg=TRUE) {
    msg<-paste("In '", as.character(sys.call(sys.parent())[[1L]]),"': ",...,"\n",sep="")
#prn(missing(show.msg))
#prn(getShowWarningMsg.ffm())
    if (missing(show.msg)) {
        if (getShowWarningMsg.ffm()) warning(msg, call.=FALSE)
        }
    else if (show.msg) warning(msg)
    }
troerrormsg<-function(...,show.msg=TRUE) {
#prn(missing(show.msg))
#prn(getShowWarningMsg.ffm())
    if (missing(show.msg)) {
        if (!exists("FFMenv")) cat(paste(...,"\n",sep=""))
        else
        if (getShowWarningMsg.ffm()) cat(paste(...,"\n",sep=""))
        }
    else if (show.msg) cat(paste(...,"\n",sep=""))
    }

# Two functions from R-devel mailing list to use 
#   in place of x %in% list when x, list are not integers.
# absolute tolerance
approxin <- function(x,list,tol=.0001) any(abs(list-x)<tol)
# relative tolerance; only exact 0 will match 0
rapproxin <- function(x,list,tol=.0001) 
    (x==0 && 0 %in% list) || any(abs((list-x)/x)<=tol,na.rm=TRUE)
# Approximate match -- vectorized -- inspired by approxin above.
match.approx  <- function(x,list,tol=.0001)
    sapply(apply(abs(outer(list,x,"-"))<tol,2,which),"[",1)   

# In response to my email to R-Devel:
#   http://tolstoy.newcastle.edu.au/R/e6/devel/09/03/1036.html
#   Here are 2 other implentations of that match.approx function
#   which use much less memory (and are faster) when the length
#   of 'x' and 'list' are long (>100, say).  The first uses
#   approx(method="const") to figure out which entries in the
#   list are just below and above each entry in x and the second
#   uses sorting tricks to do the same thing.  Then you only have
#   to figure out if the closest of those 2 entries is close enough.

#The original one above fails when tol>min(diff(sort(list))).

match.approx2 <- function(x,list,tol=.0001) {
    o1 <- rep.int(c(FALSE,TRUE), c(length(x),length(list)))[order(c(x,list))]
    o2 <- rep.int(c(FALSE,TRUE), c(length(x),length(list)))[order(c(x,list))]
 
    below <- approx(list, list, xout=x, method="constant", f=0)$y
    above <- approx(list, list, xout=x, method="constant", f=1)$y
    stopifnot(all(below<=x, na.rm=TRUE), all(above>=x, na.rm=TRUE))
    closestInList <- ifelse(x-below < above-x, below, above)
    closestInList[x<min(list)] <- min(list)
    closestInList[x>max(list)] <- max(list)
    closestInList[abs(x-closestInList)>tol] <- NA
    match(closestInList, list)
    }
match.approx3 <- function(x, list, tol=.0001){
    stopifnot(length(list)>0, !any(is.na(x)), !any(is.na(list)))
    oox <- order(order(x)) # essentially rank(x)
    i <- rep(c(FALSE,TRUE), c(length(x),length(list)))[order(c(x,list))]
    i <- cumsum(i)[!i] + 1L
    i[i > length(list)] <- NA
    i <- order(list)[i]
    leastUpperBound <- i[oox]
    i <- rep(c(TRUE,FALSE), c(length(list),length(x)))[order(c(list,x))]
    i <- cumsum(i)[!i]
    i[i < 1L] <- NA
    i <- order(list)[i]
    greatestLowerBound <- i[oox]
    closestInList <- ifelse(is.na(greatestLowerBound),
            leastUpperBound, # above max(list)
            ifelse(is.na(leastUpperBound),
                greatestLowerBound, # below min(list)
                ifelse(x-list[greatestLowerBound]<list[leastUpperBound]-x,
                        greatestLowerBound,
                        leastUpperBound)))
    if (tol<Inf) closestInList[abs(x - list[closestInList])>tol] <- NA
    closestInList
    }
    
is.wholenumber <- function(x, tolerance = .Machine$double.eps^0.5)
    return(abs(x - round(x)) < tolerance)
    # Further to above [R-Devel] 4/23/09:
    # Every value of double type which is larger than 2^53 is an even integer. 
    # The number has at least 54 digits mathematically, but only the first 53
    # of them are stored, so the last one is zero. Since rounding is done 
    # towards an even value, we get x+1 == x.
    # Double type is safe for integer values in the interval [-2^53, 2^53].
    # Inside this interval, basic arithmetic operations (+,-,*,/) on integers,
    # whose result is an integer mathematically, yield integer result also in
    # the machine.
is.zero <- function(x, tol=.Machine$double.eps^0.5) abs(x) < tol
approx.zero <- function(x, tol=.Machine$double.eps^0.5) abs(x) < tol
approxeq <- function (x, y, tol=.Machine$double.eps^0.5) abs(x-y)<tol
approx.equal <- function (x, y, tol=.Machine$double.eps^0.5) abs(x-y)<tol
vector.approx.equal <- function(x, tol=.Machine$double.eps^0.5)
    all(approx.zero(x-x[1], tol=tol))


troUnityTail <- function(Tail,tolerance= .Machine$double.eps^0.5)
    is.wholenumber(Tail-1,tolerance=tolerance)
    
is.unity <- function(x,tolerance= .Machine$double.eps^0.5)
    is.wholenumber(x-1,tolerance=tolerance)
    
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

ataRatio<-function(x) {
    # x is a vector of values corresponding to a row in a triangle
    # Result will be a vector of age-to-age ratios
    len <- length(x)
    if (len<2) return(x[0])
    x[-1]/x[-len]
    }

mcd<-function(tri,future.zeros=FALSE) {
    if (future.zeros)
        apply(tri,1,function(x) 
            ifelse(length(w<-which(x!=0)) > 0,
                    x[tail(w,1)],
                    x[1]))
    else
        apply(tri,1,function(x) 
            ifelse(length(w<-which(!is.na(x))) > 0,
                    x[tail(w,1)],
                    x[1]))
    }

right.edge.col<-function(tri,future.zeros=FALSE) {
    if (future.zeros)
        apply(tri,1,function(x)
            ifelse(length(w<-which(x!=0)) > 0,
                    tail(w,1),
                    0))
    else
        apply(tri,1,function(x)
            ifelse(length(w<-which(!is.na(x))) > 0,
                    tail(w,1),
                    0))
    }

convert.young.NAs.to.zero <- function(x, future.zeros=FALSE) {
    x[is.na(x) & col(x)<right.edge.col(x,future.zeros=future.zeros)] <- 0
    x
    }
