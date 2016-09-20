  ## Author: Dan Murphy
  ## Copyright: Daniel Murphy, chiefmurphy@gmail.com
  ## Date: 10/9/2008


 # 5 Areas in this file
 # 1. troscript
 # 2. trocreate
 # 3. trocomplete
 # 4. trostoch
 # 5. troexcel
 
##########################################################################
# troscript
##########################################################################

source('tro.r')
source('troutil.r')
source('optalpha.r')
source('troexcel.r')
source('dpo.r')
source('troreg.r')
source('xl.r')
source('mondate.r')
source('trotail.r')
source('trointerp.r')
source('troProject.r')
source('troldfs.r')
source('linterp2d.r')
tsr <- function(x) source('Mack1993_Table4.r')

troscript <- function(min.df=0) {
  # Run m<-troscript(0) for Mack's proportional sigma
  
  tro <- new.tro()

  # Get triangle source data from excel.
  tro$tri <- triGetFromExcel(
    "Mack1993_Table4.xls","Triangle")
  if (is.null(tro$tri)) stop ("NULL tri")
  
  # Convert to tdf format with properties ('triangle data frame').
  tro$tdf <- tri.as.tdf(tro$tri)
  
  # Calculate volume weighted factors
  tro$dpo <- dpoVWata(tro$tri,1.05)
  if (is.null(tro$dpo)) stop ("NULL tro")
  
  # Adjust the triangle for zero beginning-period values
#  tro$dpo$zero.adj <- "100"
#  tro$dpo$na.adj <- "100"
#  tro$dpo$zero.adj <- "avg"
#  tro$dpo$zero.adj <- "del"
  tro$tda<-as.tda(tro$tdf,tro$dpo)
  
  # Calculate optimal alpha values per FFM paper, store with dpo.
  tro$dpo <- tdaOptimalAlpha(tro$tda,tro$dpo)

  # Calculate link ratio models
  tro$LM <- tdaReg(tro$tda,tro$dpo)
  # Store the error statistics, etc.
  tro$dpo <- troStoreModelResults(tro$LM,tro$dpo)
  # Calc Mack's heuristic for the df=0 age(s)
  tro$dpo <- tdaMackSigmaHeuristic(tro$dpo,tro$tda)
#return(tro)
  # Store "selected" dpo standard errors sef and sigma (here, Mack)
  tro$dpo <- dpoSelectMackSEs(tro$dpo)
  # For Tail errors, use Mack's idea.
  tro$dpo <- troTailInterpSEsMack(tro$dpo)
  
  # Need losses to develop. We'll use triangle's most current diagonal.
  tro$diagonal <- tdf.mcd(tro$tdf)

  # Calc future part (df format) of the triangle using selected factors
  tro$ftr <- troProjectToUltimate(tro$diagonal,tro$dpo)

  # Summarize cv's and development factors
#  tro$ldfs <- troldfs(tro$diagonal,tro$dpo)
  tro$summary.byage <- troSummaryByAge(tro$diagonal,tro$dpo)
  
  # Send results to excel
  cell <- troPrintToExcel(tro) # returns reference to col A cell in next row
  
  return(tro)
  
  }

  