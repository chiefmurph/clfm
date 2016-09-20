  ## Author: Dan Murphy
  ## Copyright: Daniel Murphy, chiefmurphy@gmail.com
  ## Date: 10/9/2008
  ##        5/26/09
library(ChainLadder)
library(reshape2)
library(mondate)

#SOURCE.DIR <- "C:\\Users\\Owner\\Documents\\Trinostics\\tro\\"
SOURCE.DIR <- paste(getwd(),"/",sep="")
##########################################################################
# troscript
##########################################################################
source("/Utilities/excelRio.r")

source('troFFM.r')
source('tro.r')
source('troutil.r')

source('dpo.r')

source('troLM.r')

source('troProject.r')
source('predict.ffm.r')
source('linterp.r')

source('troWrapUp.r')
source('troSumByAge.r')
source('troresid.r')

#source('optalpha.r')
#source('troexcel.r')
#source('xl.r')
#library(xlRio)
#source('mondate.r')
#source('trotail.r')
#source('trointerp.r')
#source('troldfs.r')
#source('troGetData.r')


n.terms="3"
min.df=0
ata.err.method="Mack"
zero.int=TRUE
    # Initialize tro, global parameters in FFM environment
    tro <- troFFMInitialize() # 'troFFM.r'
    # Get triangle data, store in tda format; selected factors; 
    #   diagonal to develop.
#    tro <- troGetData(tro)

    readline("Copy triangle from excel to clipboard, hit return when ready ...")
    tro$tri <- pasteFromExcel(header = TRUE, rowheader = TRUE)
    # Convert to tdf format with properties ('triangle data frame').
    tro$tdf <- tri.as.tdf(tro$tri)
    # Get acpTable data from excel
    #tro$acpTable <-  NULL
    #if (is.null(tro$acpTable)) {
        # Try to create table from tri's acp values in its rownames
        yr<-rownames(tro$tri)
        if (!all.is.numeric(yr)) stop ("NULL acpTable")
        tro$acpTable<-data.frame(acp=yr,
            acpbegindt = as.Date(paste(yr,"-1-1",sep="")),
            acpenddt   = as.Date(paste(yr,"-12-31",sep="")),
            stringsAsFactors =FALSE)
        row.names(tro$acpTable) <- yr
    #    }

    readline("Copy 'selected' from excel to clipboard, hit return when ready ...")
    tro$dpo <- pasteFromExcel(header = TRUE, rowheader = TRUE)
    tro$dpo <- as.dpo(tro$dpo, ages = colnames(tro$tri))
    if (is.null(tro$dpo)) tro$dpo <- dpoVWata(tro$tri) # default to vol wtd

    # Expand tdf to include ata's, etc.
    tro$tda <- as.tda(tro)

    tro$diagonal <- tdf.mcd(tro$tda)

#    oa <- calcOptimalAlpha(tro$tda, tro$dpo)
    oa <- CLFMdelta(Triangle = tro$tri, selected = tro$dpo$ata[1:(ncol(tro$tri)-1)])
oa[1] <- -8

x <- tro$tri[as.character(1994:2002),"120"]
y <- tro$tri[as.character(1994:2002),"132"]
alpha <- seq(-10,30,by=.1)
plot(alpha, LRfunction(x, y, alpha), type = "l", main = "Link Ratio Function for 120-132 Development Period")
abline(h=tro$dpo["120", "ata"], lty = "dashed", col = "red")

op <- par(mfrow = c(2, 2))

i <- 2003
x <- tro$tri[as.character(i:(i+8)),"12"]
y <- tro$tri[as.character(i:(i+8)),"24"]
alpha <- seq(-100,100,by=.1)
plot(alpha, LRfunction(x, y, alpha), type = "l", main = "12-24")#, ylim = c(1.55, 1.75))
abline(h=tro$dpo["12", "ata"], lty = "dashed", col = "red")
abline(v=oa["12"], lty = "dashed", col = "blue")
a <- round(oa["12"], 2)
b <- round(tro$dpo["12", "ata"], 3)
text(a, b, labels = paste("(", a, ", ", b, ")", sep = ""), adj = c(1,1))

i <- i - 1
x <- tro$tri[as.character(i:(i+8)),"24"]
y <- tro$tri[as.character(i:(i+8)),"36"]
alpha <- seq(-100,100,by=.1)
plot(alpha, LRfunction(x, y, alpha), type = "l", main = "24-36")#, ylim = c(1.55, 1.75))
abline(h=tro$dpo["24", "ata"], lty = "dashed", col = "red")
abline(v=oa["24"], lty = "dashed", col = "blue")
a <- round(oa["24"], 2)
b <- round(tro$dpo["24", "ata"], 3)
text(a, b, labels = paste("(", a, ", ", b, ")", sep = ""), adj = c(1,1))

i <- i - 1
x <- tro$tri[as.character(i:(i+8)),"36"]
y <- tro$tri[as.character(i:(i+8)),"48"]
alpha <- seq(-100,100,by=.1)
plot(alpha, LRfunction(x, y, alpha), type = "l", main = "36-48")#, ylim = c(1.55, 1.75))
abline(h=tro$dpo["36", "ata"], lty = "dashed", col = "red")
abline(v=oa["36"], lty = "dashed", col = "blue")
a <- round(oa["36"], 2)
b <- round(tro$dpo["36", "ata"], 3)
text(a, b, labels = paste("(", a, ", ", b, ")", sep = ""), adj = c(1,1))

i <- i - 1
x <- tro$tri[as.character(i:(i+8)),"48"]
y <- tro$tri[as.character(i:(i+8)),"60"]
alpha <- seq(-100,100,by=.1)
plot(alpha, LRfunction(x, y, alpha), type = "l", main = "48-60", ylim = c(1.12, 1.16))
abline(h=tro$dpo["48", "ata"], lty = "dashed", col = "red")
a <- 0
b <- round(tro$dpo["48", "ata"], 3)
text(a, b+.001, labels = paste("(NA, ", b, ")", sep = ""), adj = c(NA,0))

par(op)
#Mack method: See below.
# Plug values for oa
#oa[is.na(oa)] <- 1
  #tro:         SUM  77548.261   209.324     454.636   500.510 10964.261 cv=0.04564922
  #ChainLadder: Sum  77548.261   209.324     453.419   499.406 10964.261
#oa[is.na(oa)] <- 4
  #tro:         SUM  77691.710   205.782     460.842   504.700 11107.710 cv=0.04543691
  #ChainLadder: Sum  77691.710   205.782     459.313   503.304 11107.710
#oa[is.na(oa)] <- -4
  #tro:         SUM  77426.755   205.596     532.370   570.690 10842.755 cv=0.0526333
  #ChainLadder: Sum  77426.755   205.596     529.080   567.623 10842.755
#oa[is.na(oa)] <- 0
  #tro:         SUM  77497.315   209.111     462.055   507.171 10913.315 cv=0.04647268
  #ChainLadder: Sum  77497.315   209.111     460.826   506.052 10913.315
    tro$dpo <- c(coef(chainladder(Triangle = tro$tri, delta = oa)), 1.000) # append the tail
    tro$dpo <- as.dpo(tro$dpo, ages = colnames(tro$tri))
    tro$dpo$oa <- c(oa, NA)

    # Calculate link ratio models
    tro$LM <- troLM(tro, min.df = min.df, ata.err.method = ata.err.method, zero.int = zero.int, tol = .0005)

    # Calc future part (tda format) of the triangle using selected factors.
    # Includes both the detail and the all-accident-years-combined summary.
    tro$ftr <- troProjectToUltimate(tro$diagonal, tro$LM, n.terms = n.terms) #,sum.only=sum.only)

    # Wrap up
    tro<- troWrapUp(tro, printxl = FALSE)

  # Future triangle values, then total risk
  # rbind(acast(tro$ftr$detail, acp ~ age, value.var = "value"), Sum = tro$ftr$acpsum$value)
  # rbind(acast(tro$ftr$detail, acp ~ age, value.var = "totalr"), Sum = tro$ftr$acpsum$totalr)

ULT <- rbind(subset(tro$ftr$detail, is.infinite(age)) , subset(tro$ftr$acpsum, subset = is.infinite(age), select = -7))
ULT$IBNR <- ULT$value - c(tro$diagonal$value, sum(tro$diagonal$value))
row.names(ULT) <- ULT$acp
ULT <- ULT[, c("value", "deltar", "gammar", "totalr", "IBNR")]

mcl <- MackChainLadder(tro$tri, alpha = 2 - oa, est.sigma="Mack", mse.method = "Independence")

y <- data.frame(
  ultimate = summary(mcl)$ByOrigin$Ultimate
  , ParamRisk = mcl$Mack.ParameterRisk[,ncol(mcl$Mack.ParameterRisk)]
  , ProcessRisk = mcl$Mack.ProcessRisk[,ncol(mcl$Mack.ProcessRisk)]
  , TotalRisk = summary(mcl)$ByOrigin$Mack.S.E
  , IBNR = summary(mcl)$ByOrigin$IBNR
  )
y <- rbind(y,
data.frame(
  ultimate = summary(mcl)$Totals["Ultimate:", ]
  , ParamRisk = mcl$Total.ParameterRisk[length(mcl$Total.ParameterRisk)]
  , ProcessRisk = mcl$Total.ProcessRisk[length(mcl$Total.ProcessRisk)]
  , TotalRisk = summary(mcl)$Totals["Mack S.E", ]
  , IBNR = summary(mcl)$Totals["IBNR", ]
  , row.names = "Sum")
  )

print(round(as.matrix(ULT), 3), scientific = FALSE, big.mark = ",")
print(round(as.matrix(y), 3), scientific = FALSE, big.mark = ",")


MM <- MackChainLadder(tro$tri, est.sigma="Mack", mse.method = "Independence")
y <- data.frame(
  ultimate = summary(MM)$ByOrigin$Ultimate
  , ParamRisk = MM$Mack.ParameterRisk[,ncol(MM$Mack.ParameterRisk)]
  , ProcessRisk = MM$Mack.ProcessRisk[,ncol(MM$Mack.ProcessRisk)]
  , TotalRisk = summary(MM)$ByOrigin$Mack.S.E
  , IBNR = summary(MM)$ByOrigin$IBNR
  )
y <- rbind(y, 
data.frame(
  ultimate = summary(MM)$Totals["Ultimate:", ]
  , ParamRisk = MM$Total.ParameterRisk[length(MM$Total.ParameterRisk)]
  , ProcessRisk = MM$Total.ProcessRisk[length(MM$Total.ProcessRisk)]
  , TotalRisk = summary(MM)$Totals["Mack S.E", ]
  , IBNR = summary(MM)$Totals["IBNR", ]
  , row.names = "Sum")
  )
#print(round(as.matrix(y), 3), scientific = FALSE, big.mark = ",")
#Mack method: Sum  77382.549   204.852     388.245   438.975 10798.549 cv=0.04065126
