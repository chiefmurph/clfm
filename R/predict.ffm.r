predict.ffm.longversion <- function(LM,newdata,n.terms="3",
                        ata.f="spline",tail.f="exponentialdecay") {
    # Un-implemented as of 9/18/09
    # This long version does all the checking:
    # LM is a single FFM model or a list of FFM models.
    # newdata is a dataframe of values, etc.
    # If LM is a single model, then its age must match newdata's age
    # If LM is a list of models, then we'll match entries in newdata by age
    #   to the correct model in the list.
    # If newdata does not have an acp column, we'll assign it one with
    #   the value '1'.
    # If newdata does not have deltar or gammar columns, then we'll
    #   assign it columns of zeroes.
    }

predict.ffm <- function(object,newdata,LM,n.terms="3",
                        ata.f="spline",tail.f="exponentialdecay") {
    if (nrow(newdata)==0) return(NULL)
    # Description:
    # Function to project the losses to the next age using
    #   the ffm model for the current age.
    # Arguments:
    # object    Object of class inheriting from "ffm". 
    #           Holds the ffm model for one age of the chain ladder method.
    #           For Mack method, the result of a weighted linear regression.
    #           Expected to have components
    #               age
    #               endage
    #               ata
    #               sef
    #               sigma
    #               alpha
    # newdata   Data frame containing losses as of object's age and deltar
    #           and gammar error estimates corresponding to those losses.
    #           Expected to have components
    #               acp
    #               age
    #               value
    #               deltar
    #               gammar
    #           If the losses are actual observations, deltar and gammar = 0.
    #           If the losses are predictions, deltar and gammar 
    #           are likely non-zero.
    # ata.f, tail.f: functions to use for interpolating link ratios and 
    #           errors -- i.e., ata, sef, sigma -- for newdata$age<>object$age
    #           ata.f used when newdata$age is bounded within the set of 
    #               LM$age's
    #           tail.f used when newdata$age > oldest LM$age
    #           
    # Details:
    # For rows of newdata not "on age" this function uses interpolation 
    #   to modify the object's statistics (coef, std. error, and sigma) 
    #   before carrying out the "fit" and "error" calculations.
    # alpha is not interpolated (variance exponent is constant throughout
    #   the development period.
    intrparms <- interpata(newdata$age, object, LM, ata.f, tail.f)
    data.frame(acp=newdata$acp,
               age=object$endage,
               value=intrparms$ata*newdata$value,
               deltar = d <- ffm.deltar.formula(intrparms$ata,
                                         intrparms$sef,
                                         newdata$value,
                                         newdata$deltar,
                                         n.terms=n.terms),
               gammar = g <- ffm.gammar.formula(intrparms$ata,
                                         intrparms$sigma,
                                         object$alpha,
                                         newdata$value,
                                         newdata$gammar),
               totalr=sqrt(d^2+g^2),
               stringsAsFactors=FALSE)
    }

interpata <- function(age, object, LM,
                        ata.f=c("spline","linear"),
                        tail.f=c("exponentialdecay","inversepower")) {
    # Technically, it's interpolation if age x <= tailage,
    #   otherwise it's extrapolation.

    ata <- rep(object$ata,length(age))
    sef <- rep(object$sef,length(age))
    sigma <- rep(object$sigma,length(age))

    # Interpo/extrapolation only necessary if x is not "on-age".
    i <- age!=object$age
    if (sum(i)>0) { # some TRUEs; ie, some ages "off-age"
        # Don't have to interpo/extrapolate all age's, only those "off-age".
        #   Thus the [i]'s in logic below.
        if (is.unity(object$ata)) {
            # When ata=1, interpo/extrapolated ata's are also 1 (set above). 
            # Also, expect sef=sigma=0, but if not,
            # interpo/extrapolate the sef and the sigma of this unity ata
            sefn <- object$sef!=0
            sigman <- object$sigma!=0
            if (sefn | sigman) {
                if (is.finite(object$endage)) {
                    # interpolate in proportion to where age lies
                    #   within object's development period
                    lambda <- (age[i]-object$age) /
                              (object$endage-object$age)
                    if (sefn) sef[i] <- lambda*object$sef
                    if (sigman) sigma[i]<-lambda*object$sigma
                    }
                else { 
                    # Extrapolate sef and sigma using exponential decay.
                    # Find the fitted values as of object's age.
                    # For each newdata row, find the extrapolated sef's and
                    #   sigma's and adjust in proportion to the ratio 
                    #   of object's values to its fitted values.
                    X <- sapply(LM,function(x) x$age)
                    # sef
                    if (sefn) {
                        Y <- log(sapply(LM,function(x) x$sef))
                        lmss <- lm(Y~X, data=data.frame(Y,X))
                        oess <- exp(predict(lmss,
                                    newdata=data.frame(X=object$age)))
                        ratio <- object$sef/oess
                        sef[i] <- ratio*exp(predict(lmss,
                                    newdata=data.frame(X=age[i])))
                        }
                    # sigma
                    if (sigman) {
                        Y <- log(sapply(LM,function(x) x$sigma))
                        lmss <- lm(Y~X, data=data.frame(Y,X))
                        oess <- exp(predict(lmss,
                                    newdata=data.frame(X=object$age)))
                        ratio <- object$sigma/oess
                        if (sigman) sigma[i] <- ratio*exp(predict(lmss,
                                            newdata=data.frame(X=age[i])))
                        }
                    }
                }
            }
        else { # object$ata != 1
            # We'll interpo/extrapolate ata and factor object's sef and sigma
            #   in proportion to the ratio of the development portion of
            #   the interpo/extrapolated ata to the development portion of
            #   object's ata.
            # We'll interpo/extrapolate ata by interpo/extrapolating the ldf
            #   pattern and getting the age-to-age factor by dividing by
            #   the ldf as of object's endage.
            # Pull age's, ata's out of LM, calculate the ldf pattern.
            ldfages <- sapply(LM,function(x) x$age)
            ldf <- rev(cumprod(rev(sapply(LM,function(x) x$ata))))
            # If object not a tail, use cubic spline or linear interpolation
            if (is.finite(object$endage)) {
                f <- if (match.arg(ata.f)=="spline")
                         splinefun(ldfages,ldf,method="natural")
                     else approxfun(ldfages,ldf,method="linear")
                # If beyond age to which this ata develops losses, ata=1
                #   (i.e., no further development).
                # Otherwise, ata = interpolated ldf / next "on-age" ldf, 
                #   where next ldf is the ldf vector w/o 1st element and 
                #   with 1.000 @ infinity.
                ata[i] <- ifelse(age[i]>object$endage, 1,
                    f(age[i])/c(ldf[-1],"Inf"=1.000)[as.character(object$endage)])
                }
            else { # next age is infinity; i.e., object is the tail.
                # Therefore,  age is on or beyond the tail and must extrapolate. 
                # Two methods to date are 
                #   "exponential decay" and "inverse power"
                #   Discussed in Sherman, PCAS 1984.
                #   Both work on the "development portion" of the ldf (ie, ldf-1)
                xs1 <- ldf>1
                if (sum(xs1)<2)
                    warning("Too few positive ata's to extrapolate beyond tail; using tail")
                else {
                    tail.f <- match.arg(tail.f)
                    if (tail.f[1]=="exponentialdecay") {
                        Y <- log(ldf[xs1]-1)
                        X <- ldfages[xs1]
                        lma <- lm(Y~X, data=data.frame(X,Y))
                        ratio <- (object$ata-1)/exp(predict(lma,newdata=data.frame(X=object$age)))
                        ata[i] <- ratio * 
                            exp(predict(lma,newdata=data.frame(X=age[i]))) + 1
                        }
                    else { # inverse power
                        Y <- log(ldf[xs1]-1)
                        X <- log(ldfages[xs1])
                        lma <- lm(Y~X, data=data.frame(X,Y))
                        ratio <- (object$ata-1)/exp(predict(lma,newdata=data.frame(X=log(object$age))))
                        ata[i] <- ratio * 
                            exp(predict(lma,newdata=data.frame(X=log(age[i])))) + 1
                        }
                    }
                }
            # no need to worry about zero in denominator since not unity ata
            lambda <- (ata[i]-1)/(object$ata-1)
            sef[i] <- lambda*object$sef
            sigma[i] <- lambda*object$sigma
            }
        }
    list(ata=ata, sef=sef, sigma=sigma)
    }

# Note that fields 'f' and 'sef' from the model and 'value' and 'deltar' 
#   from the the data.frame can be scalars or vectors.
ffm.deltar.formula <- function(f,sef,value,deltar,n.terms="3") {
    if (n.terms=="3") sqrt((sef*value)^2+(f*deltar)^2+(sef*deltar)^2)
    else sqrt((sef*value)^2+(f*deltar)^2)
    }
    
ffm.gammar.formula <- function(ata,sigma,alpha,value,gammar) {
# The following prn statements were inserted to check Manolis' 1.362 psi calculation in the paper
# Another one is at the end of this function.
#prn(alpha)
#if (alpha == 1.158) prn(ata)
#if (alpha == 1.158) prn(sigma)
#if (alpha == 1.158) prn(value)
#if (alpha == 1.158) prn(gammar)
    z <- numeric(length(value))
    idx <- !is.zero(value)
    if (sum(idx)) {
      psi <- ffm.psi.function(gammar[idx]/value[idx], alpha)
        z[idx] <- sqrt(value[idx]^alpha*ffm.psi.function(gammar[idx]/value[idx],alpha)*sigma[idx]^2 +
                       (ata[idx]*gammar[idx])^2)
        # 
        if (any(value[idx]^alpha*ffm.psi.function(gammar[idx]/value[idx],alpha)*sigma[idx]^2 +
                           (ata[idx]*gammar[idx])^2<0)) {
            x <- value[idx]^alpha*ffm.psi.function(gammar[idx]/value[idx],alpha)*sigma[idx]^2 +
                           (ata[idx]*gammar[idx])^2
            jdx <- x<0
            prn(length(x[jdx]))
            prn(x[jdx])
            prn(value[idx][jdx])
            prn(gammar[idx][jdx])
            prn(head(ata))
            prn(length(sigma))
            prn(head(sigma))
            prn(alpha)
            prn(head(value))
            prn(head(gammar))
            }
       }
    # where begin value = 0, mult prior gammar by ata
    z[!idx] <- ata[!idx]*gammar[!idx]
#if (alpha == 1.158) prn(psi)
    z
    }

ffm.psi.function <- function(cv,alpha) {
    # alpha is a scalar
    # cv can be a vector
    # For S+, might have to formulate differently.
    # Logic of integer alpha, positive and negative alpha is all wrapped up
    #   here in the function defined within sapply
    sapply(cv,function(x) {
        if (alpha>=0)
          # If alpha is close to an integer, return the
          #   exact value the for raw moment
          if (is.wholenumber(alpha))
              structure(normal.rawmoment.ratio(round(alpha),x),flag=0)
          # otherwise, linear interpolation between integer pairs
          else
            structure(linterp.scalar(
                c(x1<-floor(alpha),x2<-ceiling(alpha)),
                c(normal.rawmoment.ratio(x1,x),normal.rawmoment.ratio(x2,x)),
                  alpha),
              flag=0)
        else {
          # Negative alpha's need special processing.
          # PSI is a table read from disk. The table was generated by
          #   running simulations of normal random variables z with mean 1
          #   and standard error equal to various values of cv.
          #   For a sequence of negative values of alpha, the ratio
          #   of the mean of z^alpha to mean(z)^alpha
          #   was stored in the PSI table at that (alpha,cv pair).
          #   Refer to the R script PSI.r.
          #
          PSItable <- getPSI()
          # PSItable is a matrix indexed by alpha, cv.
          #   Run linear interpolation on both dimensions.
          psi <- linterp2d.scalar(alpha,x,PSItable)
          flag <- attr(psi,"flag")
          if (flag>0) {
              trowarning(linterp2d.outofbounds.msg(PSItable,flag,"alpha",alpha,"cv",x))
              }
          psi
          }
        } )
    }

normal.rawmoment.ratio <- function(n,cv)
    sum(sapply(seq(0,n-n%%2,by=2), function(j)
            prod(if (j>0) n:(n-(j-1)) else 1)*2^(-j/2)/factorial(j/2)*cv^j))

