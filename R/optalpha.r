## Author: Dan Murphy
## Copyright: Daniel Murphy, chiefmurphy@gmail.com
## Date: 10/9/2008
#

LRfcnPair <- function(a,B,E) sum(B^(1-a)*E)/sum(B^(2-a))
LRfcnPairD <- function(a,B,E,s) LRfcnPair(a,B,E)-s
LRdist <- function(a,B,E,s) abs(LRfcnPair(a,B,E)-s)

tdaOptimalAlpha<-function(tda, dpo, step.a = 1, tol = .0005, ...){
    # tol: S/B same as in function tdaReg in 'troLM.r'
    #   Calculate the "major" averages: VW (alpha=1), SA (2), RA (3)
    #   If selected factor is within tol of a major average *in that order*
    #       set alpha=indicated integer and we're done.
    # if dpo is pre-loaded with oa that matches the selected factors,
    #   no need to search for optimal alpha. However, all entries of 
    #   dpo$oa must be numeric, otherwise search.
    if (!is.null(dpo$oa)) {
        Ka <- attr(dpo,"Ka") # number of ata's excluding tail to search
        if (is.null(Ka)) stop("Null Ka in tdaOptimalAlpha")
        if (any(!is.na(dpo$oa[1:Ka]))) {
            if (!all(!is.na(dpo$oa[1:Ka]))) { 
                prn(dpo)
                prn(attributes(dpo))
                trowarning("some NA dpo$oa: solving for all")
                }
            else return(dpo$oa)
            }
        }

  # Calculate the optimal alpha's per the FFM paper
  # for  triangle tri and the selected link ratios ata.
  # tdf is the original triangle in df format
  # tdfBegin is the original triangle where the beginning values of 
  #   each development period may have been adjusted for 
  #   missing values or zero values.
  #   beginvalue's can still be zero if (end)value was zero.
  #   We must ignore those (0,0) observations.
  # Wish we could use optimize, but we define optimal to be the least
  #   optimal value (in absolute value), and there's no guarantee that,
  #   in case of a tie, optimize won't return the greater value.
  # Return vector of optimal alphas.
  
  if (!is.dpo(dpo)) stop("Improper development pattern object.")
  if ( !identical(attr(tda,"ages"), dpo$age) )
    stop("Triangle and development pattern ages do not match")

  ata <- dpo$ata
  ages <- dpo$age
  Ka <- attr(dpo,"Ka") # number of ata's excluding tail
  K  <- attr(tda,"K") # number of columns in triangle

  # a is the vector of alpha's to search over for optimality
  # maxa = right endpoint of search interval
  # So search interval = [0,10] on positive side, [-10,0] on negative side
  maxa <- 10
# Change 8/19/11 different positive and negative widths, coincide with PSItable
  maxaPos <- 8
  maxaNeg <- 4
  
  # Start off w/ NAs. NA's at end will be informative
  oavec <- structure(rep(NA,nrow(dpo)),names=rownames(dpo))
  # We'll flag any index where no solution exists.
  attr(oavec, "foundSolution") <- rep(TRUE, nrow(dpo))

  # If ata's are within tolerance of the "major" averages, set oa 
  #     accordingly and can ignore those ages in search for optimal value.
  #     Often, the major averages are close together, so prioritize:
  #     volume weighted, simple average, regression average
  oavec[names(tata <- RAata.tda(tda))][approxeq(tata, dpo[names(tata), "ata"], tol)] <- 0
  oavec[names(tata <- SAata.tda(tda))][approxeq(tata, dpo[names(tata), "ata"], tol)] <- 2
  oavec[names(tata <- VWata.tda(tda))][approxeq(tata, dpo[names(tata), "ata"], tol)] <- 1
  
  # For each ata age ...
  for (k in (1:Ka)[is.na(oavec[1:Ka])]) {
    # After all is said and done, if oa is still NA, it will mean that
    #   no optimal alpha was found.
    oa <- NA 
    # value      = value at age k, 
    # value.next = value at endage k
    # Ignore pairs w/ NA's in one or more cells
    rn <- which(tda$age==dpo$age[k] 
                & !is.na(tda$value) & !is.na(tda$value.next)
                & tda$value!=0
                )
    # If no row numbers identified, move along. Shouldn't happen if 
    # dpo and tda from same triangle.
    if (length(rn)>0) {
      B <- tda$value[rn]
      E <- tda$value.next[rn]
      if (0 %in% B) # Yikes! shouldn't happen
        stop("At age ",k," zeros in B: ",B[1],",",B[2],",",B[3])

      # If more than one observation at age k, do the calculations.
      # Otherwise (when only one link ratio observation) oa=1 by default.
      if (length(rn)>1) {

        # We will look through a three times until found,
        # each time increasing the search resolution of a by factor of ten.
        # 'width' = search resolution
        width <- step.a
        loopcntr <- 1

        while (loopcntr <= 3 & is.na(oa)) {

          a=seq(0,maxa,by=width)

          # positive side first
# Change 8/19/11 different positive and negative widths, coincide with PSItable
          a <- seq(0, maxaPos, by = width)
        
          # Calc (absolute) distance between the selected ata and the average 
          # link ratio for Begin value, End value, and exponents in 'a' ...
          tmp<-sapply(a,LRdist,B,E,ata[k])
          # ... and find the minimum.
          vp<-which.min(tmp)
          posroot <- a[vp]
          # See if uniroot will run within one step of vp
          # Note: it has occurred that the curvature of the function is such
          # that posroot could be at the optimal solution but interpolation
          # between the values one step away moves beyond tolerance.
          # Therefore, if the current root yields an ata within tolerance of
          # the selected ata, then we'll take that value as the solution.
          # Otherwise, we'll use uniroot.
          if (LRdist(posroot,B,E,ata[k])>tol) {
              # "One step" handled differently within interval vs at endpoints
              if (vp>1 & vp<length(a)) { # vp within interval
                # Must be of opposite sign at endpoints for uniroot to run
                if (sign(LRfcnPairD(a[vp-1],B,E,ata[k]))*
                    sign(LRfcnPairD(a[vp+1],B,E,ata[k]))<0) {
                  posroot<-uniroot(LRfcnPairD,c(a[vp-1],a[vp+1]),B,E,ata[k])$root
                  }
                }
            else {
                if (vp==1) {# vp at left endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(a[vp],B,E,ata[k]))*
                      sign(LRfcnPairD(a[vp+1],B,E,ata[k]))<0) {
                    posroot<-uniroot(LRfcnPairD,c(a[vp],a[vp+1]),B,E,ata[k])$root
                    }
                  }
                else { # vp at right endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(a[vp-1],B,E,ata[k]))*
                      sign(LRfcnPairD(a[vp],B,E,ata[k]))<0) {
                    posroot<-uniroot(LRfcnPairD,c(a[vp-1],a[vp]),B,E,ata[k])$root
                    }
                  }
                }

              # if now not within tolerance, give positive root infinite value
              if (LRdist(posroot,B,E,ata[k])>tol) posroot <- Inf
              }

          # negative side next
# Change 8/19/11 different positive and negative widths, coincide with PSItable
          a <- seq(0, maxaNeg, by = width)
          
          tmp<-sapply(-a,LRdist,B,E,ata[k])
          vn<-which.min(tmp)
          negroot <- -a[vn]
          if (LRdist(negroot,B,E,ata[k])>tol) {
              # see if uniroot will run within a step of vn
              if (vn>1 & vn<length(a)) { # within interval
                # must be of opposite sign at endpoints
                if (sign(LRfcnPairD(-a[vn+1],B,E,ata[k]))*
                    sign(LRfcnPairD(-a[vn-1],B,E,ata[k]))<0) {
                  negroot<-uniroot(LRfcnPairD,c(-a[vn+1],-a[vn-1]),B,E,ata[k])$root
                  }
                }
              else {
                if (vn==1) { # at 0, which for negative side is right endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(-a[vn+1],B,E,ata[k]))*
                      sign(LRfcnPairD(-a[vn],B,E,ata[k]))<0) {
                    negroot<-uniroot(LRfcnPairD,c(-a[vn+1],-a[vn]),B,E,ata[k])$root
                    }
                  }
                else { # vn at left endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(-a[vn],B,E,ata[k]))*
                      sign(LRfcnPairD(-a[vn-1],B,E,ata[k]))<0) {
                    negroot<-uniroot(LRfcnPairD,c(-a[vn],-a[vn-1]),B,E,ata[k])$root
                    }
                  }
                }
              # if now not within tolerance, give negative root infinite value
              if (LRdist(negroot,B,E,ata[k])>tol) negroot <- -Inf
              }
          if (posroot < -negroot) oa <- posroot
          else if (is.finite(negroot)) oa <- negroot
          
          # get ready for next try
          loopcntr <- loopcntr+1
          width <- width/10
        
          } # end of three-tries loop
        } # end of 'if more than one observation at age k' logic
      
      # Only one ata observation at age k; oa=1 by default if 
      # selected factor = observed factor, otherwise oa=unity
      else oa <- ifelse(ata[k] == E / B, 1.0, NA)
      if (is.na(oa)) {
        trowarning("No optimal alpha solution for age ",
            dpo$age[k],". Using unity. Error estimates may be unreliable.")
        oa <- 1.0
        attr(oavec, "foundSolution")[k] <- FALSE
        }
      oavec[k] <- oa
      } # end of length(rn)>0 case
    } # end of k loop
    
  # Now that we have all the optimal alpha's, go through them one more time
  # and set them equal to the closest integer if they are within 'tol'
  # of the closest integer.
  roa <- round(oavec,0)
  i <- abs(oavec-roa)<tol & !is.na(oavec)
  oavec[i] <- roa[i]
#oavec[6]<-3.05859579847101 # Manolis' alpha value for 72-84 months
  oavec
  }

tdaOptimalAlphaOld<-function(tda,dpo,step.a=1,tol=.0001){

###Perhaps use nlm or nlminb?

  # Calculate the optimal alpha's per the FFM paper
  # for  triangle tri and the selected link ratios ata.
  # tdf is the original triangle in df format
  # tdfBegin is the original triangle where the beginning values of 
  #   each development period may have been adjusted for 
  #   missing values or zero values.
  #   beginvalue's can still be zero if (end)value was zero.
  #   We must ignore those (0,0) observations.
  # Wish we could use optimize, but we define optimal to be the least
  #   optimal value (in absolute value), and there's no guarantee that,
  #   in case of a tie, optimize won't return the greater value.
  # Return vector of optimal alphas.
  
  if (!is.dpo(dpo)) stop("Improper development pattern object.")
  if ( !identical(attr(tda,"ages"), dpo$age) )
    stop("Triangle and development pattern ages do not match")

  ata <- dpo$ata
  ages <- dpo$age
  Ka <- attr(dpo,"Ka") # number of ata's excluding tail
  K  <- attr(tda,"K") # number of columns in triangle

  # a is the vector of alpha's to search over for optimality
  # maxa = right endpoint of search interval
  # So search interval = [0,10] on positive side, [-10,0] on negative side
  maxa <- 10

  # Start off w/ NAs. Then NA's at end can be informative
  dpo$oa <- NA # fills all rows of dpo

  # For each ata age ...
  for (k in 1:Ka) {
    # After all is said and done, if oa is still NA, it will mean that
    #   no optimal alpha was found.
    oa <- NA 
    # beginvalue = value at age k, end value = value at endage k
    # Ignore pairs w/ NA's in one or more cells
    rn <- which(tda$age==dpo$endage[k] 
                & !is.na(tda$beginvalue) & !is.na(tda$value)
                & tda$beginvalue!=0
                )
    # If no row numbers identified, move along. Shouldn't happen if 
    # dpo and tda from same triangle.
    if (length(rn)>0) {
      B <- tda$beginvalue[rn]
      E <- tda$value[rn]
      if (0 %in% B) # Yikes! shouldn't happen
        stop("At age ",k," zeros in B: ",B[1],",",B[2],",",B[3])
      # If more than one observation at age k, do the calculations.
      # Otherwise (when only one link ratio observation) oa=1 by default.
      if (length(rn)>1) {

        # We will look through a three times until found,
        # each time increasing the search resolution of a by factor of ten.
        # 'width' = search resolution
        width <- step.a
        loopcntr <- 1

        while (loopcntr <= 3 & is.na(oa)) {

          a=seq(0,maxa,by=width)

          # positive side first
        
          # Calc (absolute) distance between the selected ata and the average 
          # link ratio for Begin value, End value, and exponents in 'a' ...
          tmp<-sapply(a,LRdist,B,E,ata[k])
          # ... and find the minimum.
          vp<-which.min(tmp)
          posroot <- a[vp]
#if (k==1) {prn(vp);prn(posroot)}
          # See if uniroot will run within one step of vp
          # Note: it has occurred that the curvature of the function is such
          # that posroot could be at the optimal solution but interpolation
          # between the values one step away moves beyond tolerance.
          # Therefore, if the current root yields an ata within tolerance of
          # the selected ata, then we'll take that value as the solution.
          # Otherwise, we'll use uniroot.
          if (LRdist(posroot,B,E,ata[k])>tol) {
              # "One step" handled differently within interval vs at endpoints
              if (vp>1 & vp<length(a)) { # vp within interval
                # Must be of opposite sign at endpoints for uniroot to run
                if (sign(LRfcnPairD(a[vp-1],B,E,ata[k]))*
                    sign(LRfcnPairD(a[vp+1],B,E,ata[k]))<0) {
                  posroot<-uniroot(LRfcnPairD,c(a[vp-1],a[vp+1]),B,E,ata[k])$root
#if (k==1) {prn("after uniroot");prn(posroot)}
                  }
                }
            else {
                if (vp==1) {# vp at left endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(a[vp],B,E,ata[k]))*
                      sign(LRfcnPairD(a[vp+1],B,E,ata[k]))<0) {
                    posroot<-uniroot(LRfcnPairD,c(a[vp],a[vp+1]),B,E,ata[k])$root
                    }
                  }
                else { # vp at right endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(a[vp-1],B,E,ata[k]))*
                      sign(LRfcnPairD(a[vp],B,E,ata[k]))<0) {
                    posroot<-uniroot(LRfcnPairD,c(a[vp-1],a[vp]),B,E,ata[k])$root
                    }
                  }
                }

#if (k==1) prn(LRdist(posroot,B,E,ata[k]))
              # if now not within tolerance, give positive root infinite value
              if (LRdist(posroot,B,E,ata[k])>tol) posroot <- Inf
              }

          # negative side next
          tmp<-sapply(-a,LRdist,B,E,ata[k])
          vn<-which.min(tmp)
          negroot <- -a[vn]
          if (LRdist(negroot,B,E,ata[k])>tol) {
              # see if uniroot will run within a step of vn
              if (vn>1 & vn<length(a)) { # within interval
                # must be of opposite sign at endpoints
                if (sign(LRfcnPairD(-a[vn+1],B,E,ata[k]))*
                    sign(LRfcnPairD(-a[vn-1],B,E,ata[k]))<0) {
                  negroot<-uniroot(LRfcnPairD,c(-a[vn+1],-a[vn-1]),B,E,ata[k])$root
                  }
                }
              else {
                if (vn==1) { # at 0, which for negative side is right endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(-a[vn+1],B,E,ata[k]))*
                      sign(LRfcnPairD(-a[vn],B,E,ata[k]))<0) {
                    negroot<-uniroot(LRfcnPairD,c(-a[vn+1],-a[vn]),B,E,ata[k])$root
                    }
                  }
                else { # vn at left endpoint
                  # must be of opposite sign at endpoints
                  if (sign(LRfcnPairD(-a[vn],B,E,ata[k]))*
                      sign(LRfcnPairD(-a[vn-1],B,E,ata[k]))<0) {
                    negroot<-uniroot(LRfcnPairD,c(-a[vn],-a[vn-1]),B,E,ata[k])$root
                    }
                  }
                }
              # if now not within tolerance, give negative root infinite value
              if (LRdist(negroot,B,E,ata[k])>tol) negroot <- -Inf
              }
          if (posroot < -negroot) oa <- posroot
          else if (is.finite(negroot)) oa <- negroot
          
          # get ready for next try
          loopcntr <- loopcntr+1
          width <- width/10
        
          } # end of three-tries loop
        } # end of 'if more than one observation at age k' logic
      
      # Only one ata observation at age k; oa=1 by default
      else oa<-1.0
      if (is.na(oa)) {
        warning("No optimal alpha solution for period ",k,". Using unity.")
        oa<-1.0
        }
      dpo$oa[k] <- oa
      } # end of length(rn)>0 case
    } # end of k loop
    
  # Now that we have all the optimal alpha's, go through them one more time
  # and set them equal to the closest integer if they are within 'tol'
  # of the closest integer.
  roa <- round(dpo$oa,0)
  i <- abs(dpo$oa-roa)<tol & !is.na(dpo$oa)
  dpo$oa[i]<-roa[i]

  return(dpo)

  }

triAdjustForZeros <- function(tri,dpo) {



    ##############################
    #
    # Is this function necessary anymore?
    #
    ##############################
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
  triadj <- tri
  for (k in (length(tri)-1):1)
    triadj[,k] <- sapply(1:nrow(tri), function(x) 
      if (is.na(tri[[x,k]])) NA
      else
      if (tri[[x,k]]!=0) tri[[x,k]]
      else
      if (dpo$zero.ata[k]=="avg") 
        ifelse(tri[[x,k+1]]==0, NA, tri[[x,k+1]]/dpo$ata[k])
      else
      if (dpo$zero.ata[k]=="del") NA
      else
      if (is.na(tri[[x,k+1]])) NA
      else
      if (tri[[x,k+1]]!= 0) tri[[x,k+1]]/as.numeric(dpo$zero.ata[k])
      else NA
      )
  return(triadj)
  }  

