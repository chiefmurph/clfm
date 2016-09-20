troTailInterpSEsMack <- function(dpo,Tail=NULL) {
  # Implement Mack suggestion of finding the ata of same size as Tail
  # and imputing its SEs to Tail
  K <- nrow(dpo)
  if (missing(Tail))
    if (is.finite(dpo$endage[K])) stop("Missing Tail.")
    else Tail <- dpo$ata[K]
  if (Tail==1.0) {
    message("Note: ascribing zero std errors for unity tail")
    dpo$sef[K]    <- 0.0
    dpo$sigma[K]  <- 0.0
    dpo$sef2[K]   <- 0.0
    dpo$sigma2[K] <- 0.0
    dpo$oa[K]     <- 0.0
    }
  else {
#   3/10/09 Changed from spline to linear interpolation 
#    tailage <- troTailAgeInterp(dpo,Tail,type="linear")
    # 3/12/09 Changed to new function that does not need
    #   uniroot, approxfun or splinefun
    #   (S+ doesn't have explicit forms of latter 2, as R does)
    tailage <- troTailAgeInterp(dpo,Tail)
    if (attr(tailage,"flag")==1)
        warning("Tail ",Tail," < min ata ",min(dpo$ata[1:Ka]))
    else
    if (attr(tailage,"flag")==2) 
        warning("Tail ",Tail," > max ata ",max(dpo$ata[1:Ka]))
    message("Implementing Mack tail imputation, age ",
            formatC(tailage)," ",attr(dpo,"timeunits"))
    dpo$sef[K]    <- dpoInterpSE(dpo,tailage,what="sef")
    dpo$sigma[K]  <- dpoInterpSE(dpo,tailage,what="sigma")
    dpo$sef2[K]   <- dpo$sef[K]^2  
    dpo$sigma2[K] <- dpo$sigma[K]^2
    dpo$oa[K]     <- dpoInterpSE(dpo,tailage,what="oa")
    }
  return(dpo)
  }

troTailAgeInterp <- function(dpo,Tail=NULL) {
    Ka<-attr(dpo,"Ka")
    if (missing(Tail))
        if (is.finite(dpo$endage[Ka+1])) stop("Missing Tail.")
        else Tail <- dpo$ata[Ka+1]
    solve.interp.scalar(dpo$age[1:Ka],dpo$ata[1:Ka],Tail,rule=2)
    }

