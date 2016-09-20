setwd("C:\\Users\\Dan\\Desktop\\PaperCLFM2013SpringMeetig")
SOURCE.DIR <- file.path(getwd())
source('troFFM.r')
source('tro.r')
source('dpo.r')
source('optimalAlpha.r')
source('utilsSanDiego.r')
source('troProject.r')
tro <- troFFMInitialize()
triangle <- ChainLadder::RAA
selected <- attr(ChainLadder::ata(cbind(triangle, `Inf` = triangle[,10])), "vwtd")
selected
#tolerance <- .0005 # necessary?
# from troscript_mnually get data final
tro$tri <- triangle
tro$tdf <- ChainLadder:::as.data.frame.triangle(triangle, na.rm = TRUE)
tro$tdf

# Given triangle in df format, tro$tdf, return an "expanded version"
#   with calculated fields. In tro.r.
tro$tda <- as.tda(tro, acp = "origin", age = "dev")

#yr<-rownames(tro$tri)
#if (!all.is.numeric(yr)) stop ("NULL acpTable")
#tro$acpTable<-data.frame(acp=yr,
#                         acpbegindt = as.Date(paste(yr,"-1-1",sep="")),
#                         acpenddt   = as.Date(paste(yr,"-12-31",sep="")),
#                         stringsAsFactors =FALSE)
#row.names(tro$acpTable) <- yr

tro$dpo <- selected
tro$dpo <- as.dpo(tro$dpo, ages = colnames(tro$tri))
tro$dpo
(tro$diagonal <- tdf.mcd(tro$tda))
#oa <- ChainLadder::CLFMdelta(Triangle = tro$tri, head(selected, -1))
tro$dpo$oa <- optimalAlpha(tro$tda, tro$dpo)
print(tro$dpo$oa)
print(tro$dpo)
tro$diagonal <- tdf.mcd(tro$tda)
tro$LM <- troLM(tro)
tro$ftr <- troProjectToUltimate(tro$diagonal, tro$LM)
