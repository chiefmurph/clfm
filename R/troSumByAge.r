troSummaryByAge <- function(tdf,LM) {
    # When multiple rows have the same age, we'll calculate the
    #   LDFs and related statistics at the age-summarized level.
    # Group the values on the diagonal tdf by age,
    # project those to ultimate, then calculate the ratios by age
    # from that dataset.
  
    # diagonal values aggregated by age
    d<-tdfSumByAge(tdf)
    # The "accident period" which identifies unique rows
    # will be the age at which the diagonal values were summarized.
    gpddiag <- data.frame(acp     = colnames(d),
                          age     = as.numeric(colnames(d)),
                          value   = c(d),
                          stringsAsFactors=FALSE) 
    ftr <- troProjectToUltimate(gpddiag,LM)
    # Detail rows at ultimate
    m1<-ftr$detail[ftr$detail$age==Inf,c("acp","value","deltar","gammar")]
    # merge with current value and rename columns
    m <- structure(merge(gpddiag,m1,by="acp"),
                   names=c("acp","age","Current","Ultimate",
                           "ParaRisk","ProcRisk"))
    # Put in correct order, oldest age on top
    m <- m[rev(vector.order(m$age)),]
    # Tack on a Sum row 
    m <- rbind(m,
            data.frame(acp="All",
                       age="All",
                       Current=sum(gpddiag$value),
                       Ultimate=sum(m$Ultimate),
                       ParaRisk=ftr$acpsum$deltar[is.infinite(ftr$acpsum$age)],
                       ProcRisk=ftr$acpsum$gammar[is.infinite(ftr$acpsum$age)],
                       stringsAsFactors=FALSE) )
    row.names(m) <- m$acp
    # Calculated values
    m$Ldf <- m$Ultimate/m$Current
    m$TotlRisk <- sqrt(m$ParaRisk^2+m$ProcRisk^2)
    m$ParaR.df <- m$ParaRisk/m$Current
    m$ProcR.df <- m$ProcRisk/m$Current
    m$TotlR.df <- m$TotlRisk/m$Current
    m$ParaR.cv <- m$ParaRisk/m$Ultimate
    m$ProcR.cv <- m$ProcRisk/m$Ultimate
    m$TotlR.cv <- m$TotlRisk/m$Ultimate
    m
  }

trocvs <- function(tro) {
  # Calculate cv's of summarized values in tro$ldfs
  # Detail rows
  d <- merge(tro$diagonal,
             subset(tro$ftr$detail,age==Inf,select=-2), # -2: excl age column
             by="acp")
  d<-d[vector.order(d$acp),]
  # Sum row
  s<-merge(data.frame(acp="Sum",age="Sum",
                        value=sum(tro$diagonal$value),
#                        dol=mean(tro$diagonal$dol),
                        acpbegindt=mean(tro$diagonal$acpbegindt),
                        evaldt=mean(tro$diagonal$evaldt)
                        ),
           subset(tro$ftr$sum,is.infinite(age),select=-2),by="acp")
  # Put together
  m<-rbind(d,s)
  # Rename columns
  names(m)[which(names(m)=="value.x")] <-"Current"
  names(m)[which(names(m)=="value.y")] <-"Ultimate"
  names(m)[which(names(d)=="deltar")] <-"ParaRisk"
  names(m)[which(names(d)=="gammar")] <-"ProcRisk"
  names(m)[which(names(d)=="totalr")]<-"TotlRisk"
  rownames(m)<-m$age
  # Calculate cv's
  m$OS <- m$Ultimate-m$Current
  m$ParaR.cv <- m$ParaRisk/m$Ultimate
  m$ProcR.cv <- m$ProcRisk/m$Ultimate
  m$TotlR.cv <- m$TotlRisk/m$Ultimate
  # Keep a subset of the columns
  return(m[c("acp","age","Current","Ultimate","OS",
             "ParaRisk","ProcRisk","TotlRisk",
             "ParaR.cv","ProcR.cv","TotlR.cv")])
  }

trocvplot <- function(summary.byage) {
  # Plot the cv's of argument summary.byage, where summary.byage 
  # is the output of the troSummaryByAge function above.
#  maxcv<-10 
  maxcv<-1 # Keep large cv's from dominating plot
  m<-summary.byage[c("ParaR.cv","ProcR.cv","TotlR.cv")]
  m<-t(m)
  xlab="Age"
  ylab="Risk/Ultimate"
  x<-m
  x[is.nan(x)]<-0
  x[is.infinite(x)]<-0
  maxm <- max(x)
  ylim=c(0,ifelse(is.nan(maxm),maxcv,min(maxm,maxcv)))
  
  title="CVs of Ultimate Loss by Age"
  par(cex=.7)
  barplot(m,beside=T,names.arg=colnames(m),main=title,xlab=xlab,ylab=ylab,
          axis.lty=1,
          ylim=ylim
          )
  legend("topleft",legend=c("Parameter Risk",
                   "Process Risk",
                   "Total Risk"),inset=.05,
         fill=grey.colors(3))
  }

trocvplot.xl<-function(summary.byage) {
  # Set graphics device to clipboard - metafile
  win.metafile("")
  trocvplot(summary.byage)
  dev.off()
  # Paste graph from clipboard to current position on excel worksheet
  xlObjInvoke("WK","Paste")
  xlSkipRow(38)
  }