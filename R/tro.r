#
# tro.r
#
# Daniel Murphy
#
# copywrite 2009
#
new.tro <- function() list()

tdf.mcd <- function(tdf) {
    # triangle's most current diagonal in df format
    # selection based on fact that rownames(tdf)="acp.age"
    # 4-10-09: separator changed to '-' because rownames that look like
    #   "2000.12" and "2000.120" can be interpreted by R as numeric columns
    #   when read back from a database, which would be a problem because
    #   those two have the same numeric value.
    v<-attr(tdf,"maxage")
#  tdf[paste(names(v),v,sep="."),]
    tdf[paste(names(v),v,sep="-"),]
    }

troDFToMatrix <- function(tdf)
  # Here, cast does a better job of naming columns, putting in correct order
  # Note: acp's convert as factors
  cast(tdf,acp~age,value="value")

tdf.as.trimat <- function(x,what="value",rowvar="acp",colvar="age") {
    require(Hmisc)
    rownum <- match(rowvar,names(x))
    if (is.na(rownum)) {
        cat("Row id",rowvar,"nonexistent.\n")
        return(NULL)
        }
    colnum <- match(colvar,names(x))
    if (is.na(colnum)) {
        cat("Column id",colvar,"nonexistent.\n")
        return(NULL)
        }
    what <- match(what,names(x))
    if (is.na(what)) { 
        cat("Display value",what,"nonexistent.\n")
        return(NULL)
        }
    reShape(x[[what]],id=x[[rownum]],colvar=x[[colnum]])
    }

tdf.as.tri <- function(x,what="value",rowid="acp",colid="age") {
  require(Hmisc)
  rownum <- match(rowid,names(x))
  if (is.na(rownum)) {
    cat("Row id",rowid,"nonexistent.\n")
    return(NULL)
    }
  colnum <- match(colid,names(x))
  if (is.na(colnum)) {
    cat("Column id",colid,"nonexistent.\n")
    return(NULL)
    }
  what <- match(what,names(x))
  if (is.na(what)) { 
    cat("Display value",what,"nonexistent.\n")
    return(NULL)
    }
  M <- reShape(x[[what]],id=x[[rownum]],colvar=x[[colnum]])
  nams<-colnames(M)
  M <- as.data.frame(M[vector.order(rownames(M)),vector.order(colnames(M))])
  names(M)<-nams
  M
  }

matrix.as.tri <- function(M) {
  if (is.null(colnames(M))) ages <- 1:ncol(M)
  else
  if (is.logical(ages<-is.agevector(colnames(M),"convert"))) {
    cat("Column names must be increasing ages\n")
    return(NULL)
    }
  M <- as.data.frame(M)
  names(M) <- attr(M,"ages") <- ages # character
  return(M)
  }

tri.as.tdf <- function(tri,timeunits) {
    # Convert tri from df or matrix to "long" triangle data frame object
    
    # Modified 4/24/09:
    #   Now it is expected that the triangle's rownames are acp's
    #       and colnames are age's
    #   triGetFromExcel was modified to store attr(tri, "timeunits").
    #       Copy to tdf.
    # Modified 4/27/09:
    #   Remove any other processing. Calculations and attributes that 
    #       depend on other tro components will become fields and attributes
    #       of tda. Limit tdf to be simply the "long" version of
    #       the triangle

    # Check that column names are convertible to a vector of ages in 
    # increasing order.
#    rownams <- rownames(tri)
#    if (is.null(rownams)) rownams<-troDefaultFirst.acp()-1+1:nrow(tri)
#    colnams <- colnames(tri)
#    if (is.null(colnams)) {
#        colnams<-1:ncol(tri)
#        attr(tri,"timeunits")<-"years"
#        }
    if (is.null(rownames(tri))) rownams<-troDefaultFirst.acp()-1+1:nrow(tri)
    else rownams <- rownames(tri)
    if (is.null(colnames(tri))) colnams<-1:ncol(tri)
    else colnams <- colnames(tri)
    if (is.logical((age<-is.agevector(colnams,"convert")))) {
        cat("Column names ",colnams," not usable as ages\n")
        return(NULL)
        }
  
    tdf<-reshape(
        # in case tri is a matrix, not data.frame
        as.data.frame(tri), 
        # triangles are inherently in wide format
        direction="long", 
        # all columns (ages) of a triangle are time-varying
        varying=colnams, 
        # Assign the name 'value' to the data.frame column that holds the
        #   amounts observed at various points in time.
        #   (Note: if v.names were to be omitted, R would try
        #   to guess a name to use from the column names of tri, which
        #   would only work if they looked like X.12, X.24, ..., in which case
        #   R would name the column 'X'. Since tri's columns are simply numbers,
        #   that "guessing" algorithm won't work, so we'll just assign the name.)
        v.names="value",
        # idvar = variable that identifies the unique individual (acc period)
        #   whose values are observed at different points in time;
        #   name data.frame column 'acp' and get its values from rownames of tri
        idvar="acp",ids=rownams, 
        # timevar = variable that identifies the different times at which
        #   the accident period's values were observed;
        # name the data.frame column 'age' and get its values from 'age' vector
        #   assigned in previous command
        timevar="age",times=age
        )
    # Replace '.' in rownames with '-'.
    rownames(tdf)<-paste(tdf$acp,"-",tdf$age,sep="")
    # Remove NA's and put columns in desired order
    tdf<-subset(tdf,subset=!is.na(value),select=c(acp,age,value))
    # Let's see if we lost any all-NA columns
    ages <- unique(tdf$age) # numeric here
    K <- length(ages)
    if (K<ncol(tri)) 
        warning("Warning: Some age(s) had no observations. Potentially inaccurate interpolations.\n")

    # If timeunits specified, use it. Otherwise, use tri's "timeunits"
    #   attribute if it exists. If it does not, use a 'default decision
    #   algorithm' which says "If ages differ by 1, 'years', otherwise 
    #   'months' -- unless only one age, then 'years'."
    if (missing(timeunits)) attr(tdf,"timeunits") <- 
        ifelse(!is.null(attr(tri,"timeunits")), attr(tri,"timeunits"),
                ifelse(length(ages)==1, "years",
                       ifelse((ages[2]-ages[1])==1, "years", "months")))
    else attr(tdf,"timeunits") <- timeunits

    # Other properties of the triangle object ...
    attr(tdf,"ages") <- ages
    attr(tdf,"K") <- K    # number of columns=# of ages w/ at least 1 obs
    acps <-unique(tdf$acp)
    attr(tdf,"acps") <- acps[order(acps)]
    maxage  <-tapply(tdf$age,tdf$acp,max) # vector of max ages by acp
    maxage <- maxage[order(names(maxage))]
    attr(tdf,"maxage") <- maxage

    return(tdf)

    acpwidth <- attr(tri,"acpwidth")
    acpendmonth <- attr(tri,"acpendmonth")

    # Store the average accident date of each accident period
    acpyr<-sapply(strsplit(tdf$acp,"-"),function(x) x[1])
    # Assume the accident period's end date is the last day of the month
    # where the year is the row name and the month either December or 
    # is the evaluation 'age' if < 12
    # ASSUMPTION: The year = row name indicates the year in which the 
    # accident period ends.
    tdf$acpenddt <- 
        as.mondate(lastdateofmonth(as.Date(paste(acpyr,acpendmonth,"1",sep="-"))))
    tdf$acpbegindt <- tdf$acpenddt-acpwidth # width months/years earlier
    tdf$dola <- tdf$acpbegindt + ifelse(tdf$age<acpwidth,tdf$age,acpwidth)/2
    tdf$evaldt <- tdf$acpbegindt+tdf$age
    tdf$acpwidth <- tdf$acpenddt-tdf$acpbegindt
    attr(tdf,"acpwidth")  <- acpwidth

    tdf
    
    }

tri.as.tdfOld <- function(tri) {
    # Convert from df or matrix to "long" triangle data frame object

    # Check that column names are convertible to a vector of ages in 
    # increasing order.
    if (is.logical((age<-is.agevector(colnames(tri),"convert")))) {
        cat("Column names not usable as ages\n")
        return(NULL)
        }
  
    tdf<-reshape(
        # in case tri is a matrix, not data.frame
        as.data.frame(tri), 
        # triangles are inherently in wide format
        direction="long", 
        # all columns (ages) of a triangle are time-varying
        varying=colnames(tri), 
        # Assign the name 'value' to the data.frame column that holds the
        #   amounts observed at various points in time.
        #   (Note: if v.names were to be omitted, R would try
        #   to guess a name to use from the column names of tri, which
        #   would only work if they looked like X.12, X.24, ..., in which case
        #   R would name the column 'X'. Since tri's columns are simply numbers,
        #   that "guessing" algorithm won't work, so we'll just assign the name.)
        v.names="value",
        # idvar = variable that identifies the unique individual (acc period)
        #   whose values are observed at different points in time;
        #   name data.frame column 'acp' and get its values from rownames of tri
        idvar="acp",ids=rownames(tri), 
        # timevar = variable that identifies the different times at which
        #   the accident period's values were observed;
        # name the data.frame column 'age' and get its values from 'age' vector
        #   assigned in previous command
        timevar="age",times=age
        )
    # Replace '.' in rownames with '-'.
    rownames(tdf)<-paste(tdf$acp,"-",tdf$age,sep="")
    # Remove NA's and put columns in desired order
    tdf<-subset(tdf,subset=!is.na(value),select=c(acp,age,value))
    # Let's see if we lost any all-NA columns
    ages <- unique(tdf$age) # numeric here
    K <- length(ages)
    if (K<ncol(tri)) {
        cat("Warning: Some age(s) had no observations.\n")
        cat("Potentially inaccurate interpolations.\n")
        }

    # If tri was read from excel without names, then it was automatically
    #   given age's = 1,2,... and acp's = 1,2,... .
    # In that case, we'll assume that timeunits="years" and acpwidth=1.
    # Otherwise, we'll assume that timeunits="months" and acpwidth=12.
    # Obviously, this needs to be made more flexible.
    if (is.null(attr(tri,"has.names"))) {
        timeunits <- "months"
        acpwidth <- 12 
        }
    else
    if (attr(tri,"has.names")) {
        timeunits <- "months"
        acpwidth <- 12 
        }
    else {
        timeunits <- "years"
        acpwidth <- 1
        }

    # Store the average accident date of each accident period
    acpyr<-sapply(strsplit(tdf$acp,"-"),function(x) x[1])
    # Assume the accident period's end date is the last day of the month
    # where the year is the row name and the month either December or 
    # is the evaluation 'age' if < 12
    # ASSUMPTION: The year = row name indicates the year in which the 
    # accident period ends.
    # ASSUMPTION: Accident period ends on the last day of the month.
    acpendmonth <- 12 # could be 6 -- June -- for public entities, for example
    tdf$acpenddt <- 
        as.mondate(lastdateofmonth(as.Date(paste(acpyr,acpendmonth,"1",sep="-"))))
    tdf$acpbegindt <- tdf$acpenddt-acpwidth # width months/years earlier
    tdf$dola <- tdf$acpbegindt + ifelse(tdf$age<acpwidth,tdf$age,acpwidth)/2
    tdf$evaldt <- tdf$acpbegindt+tdf$age
    tdf$acpwidth <- tdf$acpenddt-tdf$acpbegindt
    # Other properties of the triangle object ...
    acps    <-unique(tdf$acp)
    acps<-acps[order(acps)]
    maxage  <-tapply(tdf$age,tdf$acp,max) # vector of max ages by acp
    maxage <- maxage[order(names(maxage))]

    attr(tdf,"acps")      <- acps
    attr(tdf,"ages")      <- ages
    attr(tdf,"K")         <- K    # number of columns=# of ages w/ at least 1 obs
    attr(tdf,"maxage")    <- maxage
#   attr(tdf,"mcd.row")   <- mcd.row
    attr(tdf,"acpwidth")  <- acpwidth
    attr(tdf,"timeunits") <- timeunits

    tdf
    
    }

tdaSetEvaldt <- function(tda,acpTable=NULL) {
    # This function is necessary in the scenario where loss data is
    #   entered in acp/age format -- the typical loss triangle.
    #   The evaluation date, therefore, is a calculation. 
    # It is customary to calculate the age of an accident period ("acp") 
    #   as the length of time from the beginning of the accident period
    #   to the evaluation date of the loss amount.
    # The beginning date of the acp is stored in tro's acpTable.
    #   That date represents the first date of exposure for that acp.
    #   A difference in time as measured by 'age' must be measured from
    #   the beginning of that day.
    #   Since tro's "mondate's" store a date as of the instant of time when
    #   that day ends, the instant when exposure begins for an accident period
    #   is the instant when the previous day ends.
    # Therefore, the age of an accident period is the length of time from
    #   the end of the day just before the period begins to the 
    #   evaluation date.
    # The age of the losses in tda is the numeric value stored in 
    #   that data.frame's 'age' column, as measured in either "months"
    #   or "years" depending on its "timeunits" attribute.
    # If the acpTable doesn't exist, we'll assume that tda$acp represents 
    #   a valid year. If it does exist, then its acp values must match
    #   those in tda.
    if (is.null(acpTable$acpbegindt))  # requires that tda$acp be a valid year
        LBdt<-as.mondate(paste(as.numeric(tda$acp)-1,"12-31",sep="-"))
    else LBdt <- '-'(as.mondate(acpTable[tda$acp,"acpbegindt"]),1,"days") # "the day before"
    '+'(LBdt,tda$age,attr(tda,"timeunits"))
#    if (attr(tda,"timeunits")=="months") LBdt+tda$age
#    else LBdt+tda$age/12 # "years" = "months"/12
    }

tdaSetDola <- function(tda,acpTable=NULL){
    # Calculate the age of the average date of loss (DOL)
    # for each row in the tda data.frame. 
    # The average DOL is the midpoint of the time period defined by
    #   the half-open/half-closed interval 
    #   (lowerbound date,upperbound date] 
    #   where those two dates are stored in the acpTable.
    # The DOL age (dola) is the length of time between the 
    #   tda evaluation date and the average date of loss.   
    if (is.null(acpTable$acpbegindt)) # requires that tda$acp be a valid year
        LBdt<-as.mondate(paste(as.numeric(tda$acp)-1,"12-31",sep="-"))
    else LBdt <- '-'(as.mondate(acpTable[tda$acp,"acpbegindt"]),1,"days") # "the day before"
    if (is.null(acpTable$acpenddt)) # requires that tda$acp be a valid year
        UBdt <- min.mondate.vectors(tda$evaldt,
                    as.mondate(paste(tda$acp,"12-31",sep="-")))
    else
        UBdt <- min.mondate.vectors(tda$evaldt,
                        as.mondate(acpTable[tda$acp,"acpenddt"]))
    DOL <- mean.mondate.vectors(LBdt,UBdt)
    tda$evaldt - DOL
    }

as.tda <- function(tro, acp = "acp", age = "age", value = "value") {

    # Given triangle in df format, tro$tdf, return an "expanded version"
    #   with calculated fields.
    
    # Start tda with a properly sorted ("ordered") version of tdf.
    #   This helps with the 'tapply' functions below.
    #   tda inherits tdf's attributes.
    tda<-tro$tdf[order(tro$tdf[[acp]],tro$tdf[[age]]),]
    attr(tda,"timeunits") <- attr(tro$tdf,"timeunits")
    
    # Calculated fields 
#    tda$evaldt<-tdaSetEvaldt(tda,tro$acpTable)
#    tda$dola <- tdaSetDola(tda,tro$acpTable)

    # Now "link up" an observation on an acp with its next observation.
    tda$age.next<-unlist(tapply(tda[[age]], tda[[acp]], 
                                function(x) c(x[-1],NA)))
    tda$value.next<-unlist(tapply(tda[[value]], tda[[acp]], 
                                function(x) c(x[-1],NA)))
#    tda$evaldt.next<-as.mondate(unlist(tapply(tda$evaldt, tda$acp, 
#                                function(x) c(x[-1],NA))))
#    tda$dola.next<-unlist(tapply(tda$dola, tda$acp, 
#                                function(x) c(x[-1],NA)))

    # Now calculate some "change" fields.
    # The incremental value of loss in a triangle is the change
    #   from the previous value, or the value itself if it is the
    #   first observation.
    # The age-to-age factor ('ata') is the ratio of the next value 
    #   and the value.
    tda$value.inc<-unlist(tapply(tda[[value]], tda[[acp]], 
                                function(x) c(x[1],diff(x))))
    tda$ata <- tda$value.next/tda[[value]]

    # A helpful label.
    tda$dev.period<-unlist(tapply(tda[[age]], tda[[acp]], 
          function(x) paste(x,c(x[-1],NA),sep="-")))
    
    # List of ages, and number of ages
    # For each acp, find all the ages, ignore the last one, then find
    #   the unique ones.
    begin.ages <- unique(unlist(tapply(tda[[age]], tda[[acp]], function(x) head(x,-1))))
    attr(tda,"begin.ages") <- begin.ages
    attr(tda,"Ka") <- length(begin.ages)
    
    maxage  <-tapply(tro$tdf[[age]],tro$tdf[[acp]],max) # vector of max ages by acp
    maxage <- maxage[order(names(maxage))]
    attr(tda,"maxage") <- maxage
    
    # rename to hardcoded names used elsewhere
    w <- which(names(tda) %in% c(acp, age, value))
    names(tda)[w] <- c("acp", "age", "value")

    # Done
    tda
    
    }

as.tdaOld <- function(tro) { # tdf,dpo=NULL) {
    # Given triangle in df format, return
    # same df with "age-to-age" additional columns:
    #   age.next = next age
    #   evaldt.next = next evaluation date
    #   value.next = next value
    #   value.inc = value.next-value if value exists, else value
    # All that work is accomplished by tdfIncremental below.
    # Returns a new dataframe sorted by acp, age (sorted in tdfIncremental).
    tda <- tdfIncremental(tro$tdf)
#    if (is.null(tro$dpo)) return(tda)
    
    # Some useful dates.
    # link together the age in tdf with acpbegindt and acpenddt in acpTable
    # Note how:
    #   1.  It is customary to measure the 'age' of an evaluation as
    #       the length of time since the beginning of the exposure period.
    #   2.  The beginning of an exposure period is the instant of time
    #       representing the close of business the day before the first
    #       day of the exposure period.
    #   3.  The 'mondate' S3 class considers a date to be that instant
    #   4.  tro$acpTable$acpbegindt is an R 'Date', so subtracting 1 day
    #       yields the 'Date' of the day before acpbegindt.
    tda$acpLBdt <- as.mondate(tro$acpTable[tda$acp,"acpbegindt"]-1)
    if (attr(tda,"timeunits")=="months") tda$evaldt<-tda$acpLBdt+tda$age
    else tda$evaldt<-tda$acpLBdt+tda$age/12 # "years" = "months"/12
    tda$acpUBdt <- min.mondate.vectors(tda$evaldt,
                        as.mondate(tro$acpTable[tda$acp,"acpenddt"]))
    tda$acpwidth <- tda$acpUBdt - tda$acpLBdt
    if (attr(tda,"timeunits")=="years") tda$acpwidth<-tda$acpwidth/12
    # average date of loss, date-of-loss age
    tda$dol <- mean.mondate.vectors(tda$acpLBdt,tda$acpUBdt)
    tda$dola <- tda$evaldt - tda$dol
    
    tda
    
    }
    
troAdjustBeginningZeros <- function() {
    
    #******************************************************
    #
    # The remainder of this function's code adjusts the beginning value
    # for zeros and NAs. If this is no longer necessary,
    # delete the rest of this code. See balboa test data.
    #
    #******************************************************
    # If dpo is NULL, value.adj will not exist, we'll just return. 
    # Otherwise,
    #   value.adj = value at current age, possibly adjusted for 
    #   value=0 or value=NA if so instructed by parameters in dpo.

    # If value=0 but value.next!=0, then regression thru origin blows up.
    
    # QUESTION: IS VALUE.NEXT=0 A PROBLEM?
     
    # Therefore, adjust those beginning values.
    # Ways to adjust:
    #   1. Back into its expected value given value.next and selected ata
    #         In this case, if value.next is zero, must delete beginning value.
    #   2. Set to NA, which has effect of deleting that development observation
    #   3. Divide value.next by a large number, which would have the effect
    #        of not significantly changing the volume weighted average
    #        In this case, if value.next is non-zero, divide.
    #        Otherwise, delete beginning value.

    if (!is.dpo(dpo)) stop("Improper development pattern object.")
    # Must be 1-1 correspondence between tri and dpo ages
    if (!identical(attr(tdf,"ages"),dpo$age))
        stop("Triangle ages", attr(tdf,"ages"),
            " and development pattern ages", dpo$age," do not match")

    tda$value.adj      <- tda$value
    tda$value.adjusted <- FALSE

    # Adjust for zero beginning values first,
    #   but only if the zero.adj column was stored in dpo
    if (!is.null(dpo$zero.adj)) { # Zero adjustment required
        rn <- which(tda$value==0)
        # if no zeros, move along
        if (length(rn)>0) for (i in rn) {
            rowi <- as.character(tda$age[i])
            if (dpo[rowi,"zero.adj"]=="avg") 
#           tda$beginvalue[i] <- tda$value[i]/dpo[tda$age[i],"ata"]
#           4-23-09: ata's now stored at age at beginning of dev period
            tda$value.adj[i] <- 
                tda$value.next[i]/dpo[rowi,"ata"]
            else
            if (dpo[rowi,"zero.adj"]=="del") 
                # after "deleting" value, will be handled as NA in next section
                tda$value.adj[i] <- NA
            else
            if (is.logical(
                # Hmisc function all.is.numeric will test if argument is
                # numeric and, if so, convert ("vector") at the same time.
                # If not numeric, returns FALSE, thus test of
                # whether large.ata is of mode 'logical'
                    large.ata<-all.is.numeric(dpo[rowi,"zero.adj"],"vector")
                )) 
                stop("Non-numeric zero adjustment value, age",rowi)
            tda$value.adj[i] <- tda$value.next[i]/large.ata
            tda$value.adjusted[i] <- TRUE
            }
        }

    # Now adjust for NA beginning values, but only if
    #   the na.adj column was stored in dpo.
    if (!is.null(dpo$na.adj)) { # NA adjustment requested
        rn <- which(tda$age!=dpo$age[1] & is.na(tda$value.adj))
        # if no NAs, move along
        if (length(rn)>0) for (i in rn) {
            rowi <- as.character(tda$age[i])
            if (dpo[rowi,"na.adj"]=="avg") 
            tda$value.adj[i] <- 
                tda$value.next[i]/dpo[rowi,"ata"]
            else
            if (dpo[rowi,"na.adj"]=="del") tda$value.adj[i] <- NA
            else
            if (is.logical(
                large.ata<-all.is.numeric(dpo[rowi,"na.adj"],"vector")
                ))
                stop("Non-numeric NA adjustment value, age",rowi)
            else tda$value.adj[i] <- tda$value.next[i]/large.ata
            tda$value.adjusted[i] <- TRUE
            }
        }

    # With the adjusted value's, calculate adjusted ata's
    tda$ata.adj<-tda$value.next/tda$value.adj
    tda
    }

as.tdaOlder <- function(tro) { # tdf,dpo=NULL) {
    # Given triangle in df format, return
    # same df with "age-to-age" additional columns:
    #   age.next = next age
    #   evaldt.next = next evaluation date
    #   value.next = next value
    #   value.inc = value.next-value if value exists, else value
    # All that work is accomplished by tdfIncremental below.
    # Returns a new dataframe sorted by acp, age (sorted in tdfIncremental).
    tda <- tdfIncremental(tro$tdf)
#    if (is.null(tro$dpo)) return(tda)
    
    # Some useful dates.
    # link together the age in tdf with acpbegindt and acpenddt in acpTable
    # Note how:
    #   1.  It is customary to measure the 'age' of an evaluation as
    #       the length of time since the beginning of the exposure period.
    #   2.  The beginning of an exposure period is the instant of time
    #       representing the close of business the day before the first
    #       day of the exposure period.
    #   3.  The 'mondate' S3 class considers a date to be that instant
    #   4.  tro$acpTable$acpbegindt is an R 'Date', so subtracting 1 day
    #       yields the 'Date' of the day before acpbegindt.
    tda$acpLBdt <- as.mondate(tro$acpTable[tda$acp,"acpbegindt"]-1)
    if (attr(tda,"timeunits")=="months") tda$evaldt<-tda$acpLBdt+tda$age
    else tda$evaldt<-tda$acpLBdt+tda$age/12 # "years" = "months"/12
    tda$acpUBdt <- min.mondate.vectors(tda$evaldt,
                        as.mondate(tro$acpTable[tda$acp,"acpenddt"]))
    tda$acpwidth <- tda$acpUBdt - tda$acpLBdt
    if (attr(tda,"timeunits")=="years") tda$acpwidth<-tda$acpwidth/12
    tda$dola <- mean.mondate.vectors(tda$acpLBdt,tda$acpUBdt)
    
    # these 3 lines can be removed if adjustments deemed unnecessary,
    # but must change logic in optimalalpha that relies on the adjusted values
    tda$value.adj      <- tda$value
    tda$value.adjusted <- FALSE
    tda$ata.adj <- tda$ata
    
    return(tda)
    
    
    
    
    #******************************************************
    #
    # The remainder of this function's code adjusts the beginning value
    # for zeros and NAs. If this is no longer necessary,
    # delete the rest of this code. See balboa test data.
    #
    #******************************************************
    # If dpo is NULL, value.adj will not exist, we'll just return. 
    # Otherwise,
    #   value.adj = value at current age, possibly adjusted for 
    #   value=0 or value=NA if so instructed by parameters in dpo.

    # If value=0 but value.next!=0, then regression thru origin blows up.
    
    # QUESTION: IS VALUE.NEXT=0 A PROBLEM?
     
    # Therefore, adjust those beginning values.
    # Ways to adjust:
    #   1. Back into its expected value given value.next and selected ata
    #         In this case, if value.next is zero, must delete beginning value.
    #   2. Set to NA, which has effect of deleting that development observation
    #   3. Divide value.next by a large number, which would have the effect
    #        of not significantly changing the volume weighted average
    #        In this case, if value.next is non-zero, divide.
    #        Otherwise, delete beginning value.

    if (!is.dpo(dpo)) stop("Improper development pattern object.")
    # Must be 1-1 correspondence between tri and dpo ages
    if (!identical(attr(tdf,"ages"),dpo$age))
        stop("Triangle ages", attr(tdf,"ages"),
            " and development pattern ages", dpo$age," do not match")

    tda$value.adj      <- tda$value
    tda$value.adjusted <- FALSE

    # Adjust for zero beginning values first,
    #   but only if the zero.adj column was stored in dpo
    if (!is.null(dpo$zero.adj)) { # Zero adjustment required
        rn <- which(tda$value.adj==0)
        # if no zeros, move along
        if (length(rn)>0) for (i in rn) {
            rowi <- as.character(tda$age[i])
            if (dpo[rowi,"zero.adj"]=="avg") 
#           tda$beginvalue[i] <- tda$value[i]/dpo[tda$age[i],"ata"]
#           4-23-09: ata's now stored at age at beginning of dev period
            tda$value.adj[i] <- 
                tda$value.next[i]/dpo[rowi,"ata"]
            else
            if (dpo[rowi,"zero.adj"]=="del") 
                # after "deleting" value, will be handled as NA in next section
                tda$value.adj[i] <- NA
            else
            if (is.logical(
                # Hmisc function all.is.numeric will test if argument is
                # numeric and, if so, convert ("vector") at the same time.
                # If not numeric, returns FALSE, thus test of
                # whether large.ata is of mode 'logical'
                    large.ata<-all.is.numeric(dpo[rowi,"zero.adj"],"vector")
                )) 
                stop("Non-numeric zero adjustment value, age",rowi)
            tda$value.adj[i] <- tda$value.next[i]/large.ata
            tda$value.adjusted[i] <- TRUE
            }
        }

    # Now adjust for NA beginning values, but only if
    #   the na.adj column was stored in dpo.
    if (!is.null(dpo$na.adj)) { # NA adjustment requested
        rn <- which(tda$age!=dpo$age[1] & is.na(tda$value.adj))
        # if no NAs, move along
        if (length(rn)>0) for (i in rn) {
            rowi <- as.character(tda$age[i])
            if (dpo[rowi,"na.adj"]=="avg") 
            tda$value.adj[i] <- 
                tda$value.next[i]/dpo[rowi,"ata"]
            else
            if (dpo[rowi,"na.adj"]=="del") tda$value.adj[i] <- NA
            else
            if (is.logical(
                large.ata<-all.is.numeric(dpo[rowi,"na.adj"],"vector")
                ))
                stop("Non-numeric NA adjustment value, age",rowi)
            else tda$value.adj[i] <- tda$value.next[i]/large.ata
            tda$value.adjusted[i] <- TRUE
            }
        }

    # With the adjusted value's, calculate adjusted ata's
    tda$ata.adj<-tda$value.next/tda$value.adj
    tda
    }

as.tdaOlderer <- function(tdf,dpo=NULL) {
  # Given triangle in df format, return
  # same df with "age-to-age" additional columns:
  #   age.prev = previous age
  #   evaldt.prev = previous evaluation date
  #   value.prev = previous value
  #   value.inc = value-value.prev if value.prev exists, else value
  # All that work is accomplished by tdfIncremental in troutil.r.
  # If dpo is NULL, beginvalue does not exist. Otherwise,
  #   beginvalue = value at previous age, possibly adjusted for
  #     zeros and/or NAs stored in dpo
  # Returns a new dataframe sorted by acp, age.
  tda <- tdfata(tdf)
  if (is.null(dpo)) return(tda)

  # If beginning value = 0, then regression thru origin blows up. 
  # Therefore, adjust those beginning values.
  # Ways to adjust:
  #   1. Back into its expected value given ending value and selected ata
  #         In this case, if ending value is zero, must delete beginning value.
  #   2. Set to NA, which has effect of deleting that development observation
  #   3. Divide ending value by a large number, which would have the effect
  #        of not significantly changing the volume weighted average
  #         In this case, if ending value is non-zero, divide.
  #         Otherwise, delete beginning value.
  if (!is.dpo(dpo)) stop("Improper development pattern object.")
  # Must be 1-1 correspondence between tri and dpo ages
  if ( !identical(attr(tdf,"ages"), dpo$age) )
    stop("Triangle ages",attr(tdf,"ages")," and development pattern ages",
        dpo$age," do not match")

  tda$beginvalue          <- tda$value.prev
  tda$beginvalue.adjusted <- FALSE

  # Adjust for zero beginning values first
  if (!is.null(dpo$zero.adj)) { # Zero adjustment required
    rn <- which(tda$beginvalue==0)
    # if no zeros, move along
    if (length(rn)>0) for (i in rn) {
      if (dpo[tda$age[i],"zero.adj"]=="avg") 
        tda$beginvalue[i] <- tda$value[i]/dpo[tda$age[i],"ata"]
      else
      if (dpo[tda$age[i],"zero.adj"]=="del") tda$beginvalue[i] <- NA
      else
      if (is.logical(
        large.ata<-all.is.numeric(dpo[tda$age[i],"zero.adj"],"vector"))) 
        stop("Non-numeric zero adjustment value")
      else tda$beginvalue[i] <- tda$value[i]/large.ata
      tda$beginvalue.adjusted[i] <- TRUE
      }
    }

  # Now adjust for NA beginning values
  if (!is.null(dpo$na.adj)) { # NA adjustment requested
    rn <- which(tda$age!=dpo$age[1] & is.na(tda$beginvalue))
    # if no NAs, move along
    if (length(rn)>0) for (i in rn) {
      if (dpo[tda$age[i],"na.adj"]=="avg") 
        tda$beginvalue[i] <- tda$value[i]/dpo[tda$age[i],"ata"]
      else
      if (dpo[tda$age[i],"na.adj"]=="del") tda$beginvalue[i] <- NA
      else
      if (is.logical(
        large.ata<-all.is.numeric(dpo[tda$age[i],"na.adj"],"vector"))) 
          stop("Non-numeric NA adjustment value")
      else tda$beginvalue[i] <- tda$value[i]/large.ata
      tda$beginvalue.adjusted[i] <- TRUE
      }
    }

  # With the adjusted beginvalue's, calculate the ata's
  tda$ata.adj<-tda$value/tda$beginvalue

#  attributes(tda) <- c(list(names=names(tda),
#                            row.names=rownames(tda)),
#                       attributes.excluding(tdf,c("names",
#                                                  "row.names",
#                                                  "mcd.row"
#                                                  )
#                                            )
#                       )
  return(tda)

  }

tdfata <- function(x){
    # Calculate and store 'ata' and 'next' fields in tdf x
    if (is.null(x$value.next)) x<-tdfIncremental(x)

    x$ata <- x$value.next/x$value
    agemod<-ceiling(log10(max(x$age))) # # of decimal places needed to print
    attr(x,"agemod")<-agemod
    x$dev.period<-paste(x$age,x$age.next,sep="-")
    x
    }

tdfIncremental <- function(tdf) {
    # Calculate the change in value of a tdf
    # The incremental change in loss, which equals the difference in
    #   cumulative loss between this age and the next age, is stored
    #   with the next age.
    # However, the link ratio, which is the ratio of cumulative loss
    #   at the next age with the loss at this age, is stored with this age.
    # Therefore, for incremental loss it's the prior value that is important
    #   and for ATA's it's the next value that is important.
    # Order the data.frame first by acp and then by age.
    # Then, for each age:
    #   We'll find the prior values by doing a right shift => 
    #       pre-pend a zero and drop the last value.
    #   We'll find the next values by doing a left shift =>
    #        drop the first value and sup-pend a "not available" NA.
    
    # Start tda with a properly sorted ("ordered") version of tdf.
    # Only $acp, $age, $value here, but inherits attributes too.
    tda<-tdf[order(tdf$acp,tdf$age),] 
#    tda$age.prev<-unlist(tapply(tda$age, tda$acp, 
#                                function(x) c(0,x[-length(x)])))
#    tda$dola.prev<-mondate(unlist(tapply(tda$dola, tda$acp, 
#                                function(x) c(0,x[-length(x)]))))
#    tda$evaldt.prev <- 
#        mondate(unlist(tapply(tda$evaldt, tda$acp, 
#                                function(x) c(NA,x[-length(x)]))))
#    tda$value.prev<-unlist(tapply(tda$value, tda$acp, 
#                                function(x) c(NA,x[-length(x)])))
    tda$value.inc<-unlist(tapply(tda$value, tda$acp, 
                                function(x) c(x[1],diff(x))))
    tda$ata<-unlist(tapply(tda$value, tda$acp, 
                                function(x) c(ataRatio(x),NA)))
    tda$dev.period<-unlist(tapply(tda$age, tda$acp, 
          function(x) paste(x,c(x[-1],NA),sep="-")))
#    tda$age.next<-unlist(tapply(tda$age, tda$acp, 
#                                function(x) c(x[-1],NA)))
    tda$value.next<-unlist(tapply(tda$value, tda$acp, 
                                function(x) c(x[-1],NA)))
    tda$age.next<-unlist(tapply(tda$age, tda$acp, 
                                function(x) c(x[-1],NA)))
prn(tda)
prn(tda$evaldt)
prn(tda$acp)
#prn(tapply( (tda$evaldt, tda$acp,function(x) c(x[-1],NA))))
    tda$evaldt.next<-unlist(tapply(tda$evaldt, tda$acp, 
                                function(x) c(x[-1],NA)))
#    tda$ata <- tda$value.next/tda$value
#    agemod<-ceiling(log10(max(tda$age))) # # decimal places needed to print
#    attr(tda,"agemod")<-agemod
#    tda$dev.period<-paste(tda$age,tda$age.next,sep="-")
    tda
    }
    
tdfataOld<-function(x){ # "old" used "prev": stored ata's @ ending (=next) age
  # Calculate and store 'ata' and 'prev' fields in tdf x
  if (is.null(x$value.prev)) x<-tdfIncremental(x)
  i<-!is.na(x$value.prev)
  x$ata[i]<-x$value[i]/x$value.prev[i]
  agemod<-ceiling(log10(max(as.numeric(x$age))))
  attr(x,"agemod")<-agemod
  x$dev.period<-paste(x$age.prev,x$age,sep="-")
  x
  }

troDevelopDiagonal <- function(diagonal,dpo) {
  tro<-list()
  tro$ftr <- troProjectToUltimate(diagonal,dpo)
  # Send results to excel
  cell <- troPrintToExcel(tro) # returns reference to col A cell in next row
  return(tro)
  }
