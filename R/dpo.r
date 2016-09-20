as.dpo <- function (ata,ages=NULL){
    # ata is a vector of age-to-age factors with Tail as last element
    # ages is a vector of ata ages
    # ata can also have two attributes:
    #   names: the beginning ages of the development periods, like 12, 24, ... 
    #          can also be 12-24, 24-36, ... .
    #   Name: pattern name, e.g., "Selected", that can be printed in reports
    #
    # If development periods not provided in names(ata),
    # default to seq_along(ata), indicating age of losses at the beginning of 
    # the development period. 
  
    
  # Value returned: a development pattern object

    # A development pattern object is a data frame with components:
    #   devperiod: look like "age-endage"
    #   ata: age-to-age factors, aka RTR (report-to-report) factors 
    
  
  
  
  #         or link ratios; Tail in last entry
    
  #   ldf: "to-ultimate" factors (cumulative product of ata's)
    #   age: ages from which to develop losses;
    #   endage: ages to which to develop losses; Inf in last entry
    #   agendx: row number in dpo of that age
    #   nextagendx: row number of next age
    #         agendx and nextagendx are useful for stepping 
    #         forward and backward through a dpo
    #
    # attributes:
    #   Name: useful attribute for reports
    #   Ka:   # of link ratios excluding tail; same as nrow(dpo)-1
    #   Tail: same as dpo$ata[nrow(dpo)]
    #
    # other:
    #   rownames=ages; important so can look up ata's by age
    #
  
    if (is.null(ata)) stop("NULL ata")
    if (length(ata)==0L) stop("zero length ata")
  
    K  <- length(ata)
    Ka <- K-1

    # If ages is NULL, then see if development period ages 
    #   are stored in ata's names. 
    #   If if ata's 'names' is empty, use sequential integers.
    if (missing(ages))
        if (is.null(names(ata))) nm <- as.character(1:K)
        else nm <- names(ata)
    else nm<-as.character(ages)

    # ata "names" nm can be a vector of numbers, or a vector
    # of development periods like "12-24","24-36" etc.
    # split.ata.periods returns a 2-col matrix (beginage,endage)
    nx <- split.ata.periods(nm)
    beginage<-nx[,1]
    endage<-nx[,2]
    devperiod <- paste(beginage,endage,sep="-")
    agendx     <- 1:K
    nextagendx <- if (K>1) c(2:K,K) else K
    ldf <- rev(cumprod(rev(ata)))
print(devperiod)
print (ata)
print(ldf)
print(beginage)
print(endage)
print(agendx)
print(nextagendx)

    dpo<- data.frame(devperiod=devperiod,
                     ata=ata,
                     ldf=ldf,
                     age=beginage,
                     endage=endage,
                     agendx=agendx,
                     nextagendx=nextagendx,
                     stringsAsFactors=FALSE)
    attr(dpo,"Name") <- attr(ata,"Name")
    attr(dpo,"Ka")   <- Ka
    attr(dpo,"Tail") <- ata[K]
    rownames(dpo) <- beginage # so can lookup ata's, etc. by age
    class(dpo) <- c("dpo",class(dpo))
    return(dpo)
    #  attr(dpo,"acpwidth") <- acpwidth
    #  attr(dpo,"timeunits") <- timeunits
    }

dpoVWata <- function(tri,Tail=1.000) {
    if (!is.null(colnames(tri)))
        structure(as.dpo(c(VWata(tri),Tail),colnames(tri)),
                  "Name"="All-Yr Vol Wtd")
    else
        structure(as.dpo(c(VWata(tri),Tail),1:ncol(tri)),
                  "Name"="All-Yr Vol Wtd")
    }

dpoSAata <- function(tri,Tail=1.000) {
    if (!is.null(colnames(tri)))
        structure(as.dpo(c(SAata(tri),Tail),colnames(tri)),
                  "Name"="All-Yr Simple Avg")
    else
        structure(as.dpo(c(SAata(tri),Tail),1:ncol(tri)),
                  "Name"="All-Yr Simple Avg")
    }

dpoRAata <- function(tri,Tail=1.000) {
    if (!is.null(colnames(tri)))
        structure(as.dpo(c(RAata(tri),Tail),colnames(tri)),
                  "Name"="All-Yr Regression Avg")
    else
        structure(as.dpo(c(RAata(tri),Tail),1:ncol(tri)),
                  "Name"="All-Yr Regression Avg")
    }

is.dpo <- function(obj) inherits(obj,"dpo")

lookup.dpo <- function(dpo,age,what="ata") dpo[as.character(age),what]
  
summary.dpo <- function(dpo,...) {
  # show ata's and ldf's and reverse order of the rows at the end
  cbind(ATA=dpo$ata,Age=dpo$age,LDF=dpo$ldf)[nrow(dpo):1,]
  }
  
age.ndx.dpo <- function(x,age) {
  if (missing(age)) return (1:nrow(x))
  return (match(age,x$age))
  }

prevdpoage <- function(dpo,age) dpo$age[match(age,dpo$endage)]

dpoModelIndex <- function(dpo,age) {
  if (missing(age)) return(NULL)
  sapply(age, function(x) which(dpo$endage>x)[1])
  }

dpoSelectMackSEs <- function(dpo) {
  dpo$sef <- dpo$sfMack
  dpo$sef2 <- dpo$sfMack^2
  dpo$sigma <- dpo$srMack
  dpo$sigma2 <- dpo$srMack^2
  return(dpo)
  }

dpoInterpSE <- function(dpo,age,what=c("sef2","sigma2","oa","sef","sigma")) {
    if (is.logical(agen<-is.agevector(age,"convert"))) 
        stop("age", head(age)," ... not usable as age vector") 
    what <- match.arg(what)
#   approxfun(dpo$age,dpo[[what]],method="constant",rule=2)(agen)
    # Changed 3/11/09 to linear interpolation
#   approxfun(dpo$age,dpo[[what]],method="linear",rule=2)(agen)
    # Changed 3/12/09 to approx (S+ has no approxfun)
    approx(dpo$age,dpo[[what]],xout=agen,method="linear",rule=2)$y
    }
  
