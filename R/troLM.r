trowarning <- function(...) warning(...)
create.models.open<-function(data.cum, closed, opt.alpha=TRUE, zero.a=ifelse(opt.alpha,FALSE,TRUE), alpha=ifelse(opt.alpha,0,1), min.df=0) {
    # Returns a list of lists of models:
    #   For each age, calculates the models that develop open claims at
    #       that age to ultimate
    lapply(1:(ncol(data.cum)-1), function(j) {
           create.models(array(data.cum[!closed[,j], j:ncol(data.cum)], 
                               dim = c(sum(!closed[,j]),ncol(data.cum) - j + 1),
                               dimnames=list(dimnames(data.cum)[[1]][!closed[,j]],dimnames(data.cum)[[2]][j:ncol(data.cum)])
                              ),
                         array(closed[!closed[,j],j:ncol(data.cum)],
                               dim = c(sum(!closed[,j]),ncol(data.cum) - j + 1)
                              ),
                         opt.alpha=opt.alpha, zero.a=zero.a, 
                         alpha=alpha, min.df=min.df)
           }
          )
    }

create.models<-function(data.cum, closed, open.only=FALSE, opt.alpha=TRUE, zero.a=ifelse(opt.alpha,FALSE,TRUE), alpha=ifelse(opt.alpha,0,1), min.df=0) {
#
# if min.df>0, will combine periods for small models (df<min) into a 
#   simultaneous regression.
# 
    # Function to minimize is the one that calculates the AIC of a model over
    #   all possible values of alpha; model, qr=F to speed execution
    f.no.a<-function(alpha,X,Y) AIC(lm(Y~X+0,w=X^-alpha,model=FALSE,qr=FALSE))
    f.w.a<-function(alpha,X,Y) AIC(lm(Y~X,w=X^-alpha,model=FALSE,qr=FALSE))
    ages <- structure(as.numeric(colnames(data.cum)),names=colnames(data.cum))
    models<-vector("list",length(ages)-1)
    names(models) <- names(ages[-1])
#    data.incr <- t(apply(data.cum,1,diff)) # regress incremental on cumulative
    # regress incrementals=column differences on cumulatives
    data.incr <- t(diff(t(data.cum)))
    for (i in 1:(length(ages)-1)) {
        models[[i]]$i <- i
        models[[i]]$age <- ages[i+1] # age of Y; "target" age
        models[[i]]$period <- paste(names(ages[i]),names(ages[i+1]),sep="-")
        if (nrow(data.cum)==0L) {
            models[[i]]$nobs <- 0L
            models[[i]]$situation <- "no data"
            }
        else {
            models[[i]]$situation <- "lm" # assume enough data to run lm
            exclude <- is.na(data.cum[,i]) | is.na(data.incr[,i])
            if (open.only) exclude <- exclude | closed[,i]
            X<-data.cum[!exclude,i]
            Y<-data.incr[!exclude,i]
            models[[i]]$nobs <- length(Y)
            models[[i]]$X <- X
            models[[i]]$Y <- Y
            all.closed <- sum(!closed[!exclude,i],na.rm=TRUE)==0
            if (all.closed) models[[i]]$situation <- "all.closed"
            else
            if (all(is.na(Y))) models[[i]]$situation <- "all.na"
            else
            if (models[[i]]$nobs==1L) models[[i]]$situation <- "only 1 obs"
            else
            if (models[[i]]$nobs==2L) models[[i]]$situation <- "only 2 obs"
            else
            if (sum(closed[!exclude,i])/length(Y)>.9995)
                models[[i]]$situation <- "essentially.all.closed"
            else
            if (vector.approx.equal(Y[!is.na(Y)])) models[[i]]$situation <- "all.approx.equal"
            }
        # situation identified. can run lm?
        if (models[[i]]$situation=="lm") {
            if (opt.alpha) {
                soln.no.a<-optimize(f.no.a,X=X,Y=Y,lower=-2,upper=8) # find optimal alpha
                soln.w.a<-optimize(f.w.a,X=X,Y=Y,lower=-2,upper=8) # find optimal alpha
                }
            else {
                soln.no.a<-list(minimum=alpha,objective=f.no.a(alpha,X,Y))
                soln.w.a<-list(minimum=alpha,objective=f.w.a(alpha,X,Y))
                }
            if (zero.a | soln.no.a$objective<=soln.w.a$objective) {
                no.a <- TRUE
                alpha<-soln.no.a$minimum
                aic<-soln.no.a$objective
                # Now run lm again and keep the qr 
                LM<-lm(Y~X+0,w=X^-alpha,model=TRUE,qr=TRUE)
                }
            else {
                no.a <- FALSE
                alpha<-soln.w.a$minimum
                aic<-soln.w.a$objective
                LM<-lm(Y~X,w=X^-alpha,model=TRUE,qr=TRUE)
                }
            sLM <- summary(LM)
            models[[i]]$LM <- LM
            models[[i]]$alpha = alpha
            models[[i]]$aic = aic
            if (no.a) models[[i]]$coefs <- matrix(c(0,0,sLM$coef[,1:2]),nrow=2,byrow=TRUE,dimnames=list(c("(Intercept)","X"),c("Estimate","Std. Error")))
            else models[[i]]$coefs <- sLM$coef[,1:2]
            models[[i]]$sigma <- sLM$sigma
            models[[i]]$cov.mat <- sLM$cov.unscaled
            models[[i]]$df <- sLM$df[2]
            models[[i]]$residuals <- sLM$residuals
            models[[i]]$r.squared <- sLM$r.squared
            models[[i]]$adj.r.squared <- sLM$adj.r.squared
            models[[i]]$fstatistic <- sLM$fstatistic
            }
        else { # identity model
            models[[i]]$LM <- NA
            models[[i]]$alpha = NA
            models[[i]]$aic = NA
            models[[i]]$coefs <- matrix(c(0,1,0,0),nrow=2,dimnames=list(c("(Intercept)","X"),c("Estimate","Std. Error")))
            models[[i]]$sigma <- 0
            models[[i]]$cov.mat <- NA
            models[[i]]$df <- NA
            models[[i]]$residuals <- NA
            models[[i]]$r.squared <- NA
            models[[i]]$adj.r.squared <- NA
            models[[i]]$fstatistic <- NA
            }
        }
    if (min.df>0) models<-min.df.func(models,min.df=min.df,opt.alpha=opt.alpha,alpha=alpha)
    models
    }
min.df.func <- function(mods,min.df=3,opt.alpha=TRUE,alpha=ifelse(opt.alpha,0,1)) {
    f.no.a<-function(alpha,X,Y) AIC(lm(Y~X+0,w=rowSums(X)^-alpha,model=FALSE,qr=FALSE))
    if (min.df<3) stop("min.df < 3 not allowed. Execution terminated.")
    ata.err.method <- "min.df"
    dfs <- sapply(mods, function(x) ifelse(is.na(x$df),0,x$df))
    nobs <- sapply(mods, function(x) ifelse(is.na(x$nobs),0,x$nobs))
    no.lm <- sapply(mods, function(x) x$situation!="lm")
    cum.dfs<-rev(cumsum(rev(dfs)))
    w<-which(cum.dfs>min.df)
    if (length(w)==0L) return(mods) # no aggregation of models meets min.df
    if (tail(w,1) == length(mods)) return(mods) # no agg necessary, last model large enough
    w <- tail(w,1):length(mods)
    w <- w[nobs[w]>0L]
    # all but last model followed by zeroes
#    X <- do.call(rbind, lapply(1:length(w),function(i) {
#        X<-matrix(0,ncol=length(w),nrow=nobs[w[i]])
#        X[,i]<-mods[[w[i]]]$X
#        X
#        }))
    X <- do.call(rbind, lapply(w,function(i) {
        X<-matrix(0,nrow=nobs[i],ncol=length(w),dimnames=list(NULL,w))
        X[,as.character(i)]<-mods[[i]]$X
        X
        }))
    Y <- do.call(c, lapply(mods[w], function(x) x$Y))
    # Now run the multiple regression
    if (opt.alpha) {
        soln.no.a<-optimize(f.no.a,X=X,Y=Y,lower=-2,upper=8) # find optimal alpha
        }
    else {
        soln.no.a<-list(minimum=alpha,objective=f.no.a(alpha,X,Y))
        }
    no.a <- TRUE
    alpha<-soln.no.a$minimum
    aic<-soln.no.a$objective
    # Now run lm again and keep the qr 
    LM<-lm(Y~X+0,w=rowSums(X)^-alpha,model=TRUE,qr=TRUE)
    sLM <- summary(LM)
    # Store the results back in the applicable models
    nobs<-nobs[w]
    cumnobs<-cumsum(nobs)
    for (i in 1:length(w)) {
        mods[[w[i]]]$situation <- "df>min"
        mods[[w[i]]]$LM <- LM
        mods[[w[i]]]$alpha = alpha
        mods[[w[i]]]$aic = aic
        mods[[w[i]]]$coefs <- matrix(c(0,0,sLM$coef[i,1:2]),nrow=2,byrow=TRUE,dimnames=list(c("(Intercept)","X"),c("Estimate","Std. Error")))
        mods[[w[i]]]$sigma <- sLM$sigma
        mods[[w[i]]]$cov.mat <- sLM$cov.unscaled
        mods[[w[i]]]$df <- sLM$df[2]
        mods[[w[i]]]$residuals <- sLM$residuals[(cumnobs[i]-nobs[i]+1):cumnobs[i]]
        mods[[w[i]]]$r.squared <- sLM$r.squared
        mods[[w[i]]]$adj.r.squared <- sLM$adj.r.squared
        mods[[w[i]]]$fstatistic <- sLM$fstatistic
        }
    mods
    }

summary.models <- function(mods) {
    ages <- attr(mods,"ages")
    do.call(rbind, lapply(mods, function(x) 
            data.frame(period=x$period,situation=x$situation,
                       a=x$coefs[1,1],se.a=x$coefs[1,2],b=x$coefs[2,1],se.b=x$coefs[2,2],
                       sigma=x$sigma,nobs=x$nobs,alpha=x$alpha,df=x$df,
                       r2.adj=x$adj.r.squared,
                       aic=x$aic)
        ))
    }
plot.models <- function(mods) {
    library(lattice)
#    ndx <- which(sapply(mods,function(x) x$situation=="lm"))
    ndx <- which(sapply(mods,function(x) !is.na(x$df)))
    w<-do.call(rbind,lapply(mods[ndx],function(x) data.frame(residuals=x$residuals,age=x$age,period=x$period)))
    x11()
    densityplot(~residuals|period,data=w,scales="free",main=deparse(substitute(mods)))
    }

troLM <- function(tro,min.df=0,
                      ata.err.method=c("Mack","loglinear"),
                      tail.err.method=c("Mack","loglinear"),
                      unity.tail.method=c("Zero","Mack","loglinear"),
                      zero.int=TRUE,
                      zero.rm=FALSE, ...) {
    # zero.int=TRUE -> no intercept in models
    # zero.rm: remove zeroes, NA's
    #   if TRUE, exclude the observations from the regression estimation
    #   if FALSE, leave the obs in, but change the beginning value
    #               to a infinitesimally small number, which will
    #               make the by-claim detail ultimate = ay agg ultimate
    # Calculate optimal alpha values per FFM paper
    #   unless already "pre-loaded"
    oa <- tro$dpo$oa
    if (is.null(oa)) oa <- tdaOptimalAlpha(tro$tda,tro$dpo)
    # Calculate link ratio models consistent with those oa's
    tdareg <-
    tdaReg(tro$tda,tro$dpo,oa,min.df=min.df,
                              ata.err.method=ata.err.method,
                              tail.err.method=tail.err.method,
                              unity.tail.method=unity.tail.method,
                              zero.int=zero.int,
                              zero.rm=zero.rm, ...)
    }

troStoreModelResults <- function(LM,dpo) {
    # 'sf': std error of the link ratio estimate f
    # 'sr': std error of the residuals ("sigma")
    for (x in names(LM)) {
        lmsum <- summary(LM[[x]])
        dpo[x,"f"]<-lmsum$coef[,"Estimate"]
        dpo[x,"sf"]<-lmsum$coef[,"Std. Error"]
        dpo[x,"sr"]<-lmsum$sigma
        dpo[x,"df"] <-LM[[x]]$df
        dpo[x,"formula"] <-
            paste(as.character(formula(LM[[x]]))[c(2,1,3)],collapse="")
        dpo[x,"AIC"] <- AIC(LM[[x]])
        }
    return(dpo)
    }
  
tdaReg <- function(tda,dpo,oa,min.df=0,
                              ata.err.method=c("Mack","loglinear"),
                              tail.err.method=c("Mack","loglinear"),
                              unity.tail.method=c("Zero","Mack","loglinear"),
                              zero.int=TRUE,
                              zero.rm=FALSE,
                              tol=.00001) {
    # tol: S/B same as in function tdaOptimalAlpha in 'optalpha.r'

    # Given tda and optimal alphas in oa, run regression models
    if (!is.agevector(attr(tda,"ages")))
        stop("tda 'ages' attribute must be an age vector")
    ata.err.method <- match.arg(ata.err.method)
    tail.err.method <- match.arg(tail.err.method)
    unity.tail.method <- match.arg(unity.tail.method)
    # Need 'X' for "development period" string
    xnames<-paste("X",dpo$age,sep="")
    ynames<-paste("X",dpo$endage,sep="")
.f<-function(k) {
    alpha <- oa[k]
    dat <- tdareg.gatherdat(tda, dpo$age[k], dpo$endage[k], alpha, zero.rm)
    nobs<- nrow(dat)
    if (nobs>0) {
        # Assign variate values to R variables with the names 'xname'.
        assign(xnames[k],dat$X)
        assign(ynames[k],dat$Y)
        if (zero.int)
            # This is how we get formulas like "X24 ~ X12 + 0"
             frmla <- formula(paste(ynames[k],"~",xnames[k],"+0",sep=""))
        else frmla <- formula(paste(ynames[k],"~",xnames[k],     sep=""))
        if (all(approx.equal(dat$X,dat$Y))) 
            LM <- list (ata=1,sef=0,sigma=0,df=nobs)
        else { 
            if (alpha!=0) {
                # The weights will be a simple vector of values
                X.to.neg.alpha <- dat$X^-alpha
                # this form shows lm(formula = frmla, weights = X.to.neg.alpha)
                LM <- lm(frmla,weights=X.to.neg.alpha)
                }
            else LM <- lm(frmla)
            slm <- summary(LM)
            LM$ata <- slm$coef[,"Estimate"]
            LM$sef <- slm$coef[,"Std. Error"]
            LM$sigma <- slm$sigma
            }
        LM$age <- dpo$age[k]
        LM$endage <- dpo$endage[k]
        LM$alpha <- alpha
        LM$xname <- xnames[k]
        LM$yname <- ynames[k]
        LM$nobs <- nobs
        # Use these "key fields" to match up LM results with tda
        LM$keyfields <- dat[-1:-2] # all but X and Y
        class(LM) <- c("ffm",oldClass(LM))
        }
    else { # no observations for this selected ata
        ffmInfoMsg("No obs for period",
                   paste(dpo$age[k],dpo$endage[k],sep="-"),
                   "; ata=", dpo$ata[k])
        LM <- list(age=dpo$age[k],
                   endage=dpo$endage[k],
                   ata=dpo$ata[k],
                   sef=0,
                   sigma=0,
                   df=0,
        # 3/26/10 changed mean to 0
                   alpha=0)
#                   alpha=mean(oa,na.rm=TRUE))
        LM$nobs  <- 0
        LM$xname <- xnames[k]
        LM$yname <- ynames[k]
        class(LM)<-"ffm"
        }
    if (!approxeq(LM$ata,dpo$ata[k],tol=tol)) {
#        trowarning("tdaReg: calculated ata ",LM$ata," at age ",dpo$age[k]," set to selected value ", dpo$ata[k])
        trowarning("tdaReg: calculated ata ",LM$ata," at age ",dpo$age[k]," set to selected value ", dpo$ata[k], "; alpha=", round(alpha, 3))
#prn(LM$ata)
#prn(dpo$ata[k])
#prn(alpha)
#prn(dat)
        LM$ata <- dpo$ata[k]
        }
    LM
    }
    nLM <- nrow(dpo)
    LM<-lapply(1:nLM,.f)
    # To reference models by age, eg, L[["12"]]:
    names(LM)<-dpo$age
    class(LM)<-c("ffmlist","list")
#return(LM)
    
    # 1. Make sure all models in LM are complete
#CheckLM
    # Does the beginage of a model's development period = endage of
    # the previous model's development period?
    if (nLM>1) for (i in 2:nLM) {
        if (LM[[i]]$age!=LM[[i-1]]$endage) 
            stop("Mismatched ages in model i=",i,": age[i]=",
                    LM[[i]]$age,", endage[i-1]=",LM[[i-1]]$endage)
        }
    # 2. If the last model is not a tail, append a unity tail.
    if (is.finite(LM[[nLM]]$endage)) {
        ffmInfoMsg("Last model not a tail, appending errorless unity tail.")
        tailage <- LM[[nLM]]$endage
        nLM <- nLM+1
        LM[[as.character(tailage)]] <- list(age=tailage,
                 endage=Inf,
                 ata=1.000,
                 sef=0,
                 sigma=0,
                 df=0,
                 alpha=0)
        LM[[nLM]]$nobs  <- 0
        LM[[nLM]]$xname <- paste("X",tailage,sep="")
        LM[[nLM]]$yname <- "XInf"
        class(LM[[nLM]])<-"ffm"
        }

    # 3. Do all models have ata's?
    for (i in 1:nLM)
        if (is.na(LM[[i]]$ata)) stop("model",i," ata is NA")

    # 4. Do all non-tail models have sef's, sigma's, and alpha's?

    # 4a. Check ata's first (i<nLM)
    #       idx will be T when something missing
    idx <- sapply(1:(nLM-1), function(i)
       is.null(LM[[i]]$sef) | is.null(LM[[i]]$sigma) | is.null(LM[[i]]$alpha) |
       is.nan(LM[[i]]$sef)  | is.nan(LM[[i]]$sigma)  | is.nan(LM[[i]]$alpha)
                 )
    if (sum(idx)>0) { # something missing somewhere
        # For missing item, estimate using one of the recognized methods.
        ata.err.method <- match.arg(ata.err.method)
        # min.df >0 takes priority
        if (min.df>0) { #LM<-troMultipleReg(LM,min.df=min.df)
            ata.err.method <- "min.df"
            #troMultipleReg <- function(LM,min.df=3) {
            # We've finished looping through ata's.
            # From least to most mature age, loop through LM until find first
            #     model with degrees of freedom < min.df.
            # Group that age and all more mature ages into one multiple regression.
            # So as to not have to worry about the model that matches the df being
            #     one for which the Mack heuristic was run, we'll arbitrarily
            #     not allow min.df to be anything less than 3.
            if (min.df<3) stop("min.df < 3 not allowed. Execution terminated.")

            # indices of LM's models whose degrees of freedom < min.df
            small.models <- which((sapply(LM,function(x) x$df)<min.df)&
                                  (y<-sapply(LM,function(x) x$nobs)>0))
            # If no df's < min.df or only one small model, nothing to do.
            if (length(small.models)<2) return(LM)

            X <- Y <- W <- numeric(0)
            totobs <- sum(sapply(LM[small.models],function(x) x$nobs))
            # all but last model followed by zeroes
            for (i in head(small.models,-1)) {
                # Build a matrix that looks like
                # y x 0 0 0
                # y x 0 0 0
                # y x 0 0 0
                # y x 0 0 0
                # y 0 x 0 0
                # y 0 x 0 0
                # y 0 x 0 0
                # y 0 0 x 0
                # y 0 0 x 0
                # y 0 0 0 x
                # Notice that each set of x's is followed by 10 zeros
                if (LM[[i]]$nobs>0) {
                    Y <- c(Y,unlist(LM[[i]]$model[1]))
                    X <- c(X,unlist(LM[[i]]$model[2]),rep(0,totobs))
                    W <- c(W,unlist(LM[[i]]$model[3]))
                    }
                }
            i <- tail(small.models,1)
            Y <- c(Y, unlist(LM[[i]]$model[1]))
            X <- c(X, unlist(LM[[i]]$model[2]))
            W <- c(W, unlist(LM[[i]]$model[3]))
            dim(X) <-c(totobs, length(small.models))
            # Now run the multiple regression
            slm<-summary(lm(Y~X+0, weights=W))
            # Store the results back in the applicable models
            for (i in 1:length(small.models)) {
                LM[[small.models[i]]]$ata <- slm$coef[i,"Estimate"]
                LM[[small.models[i]]]$sef <- slm$coef[i,"Std. Error"]
                LM[[small.models[i]]]$sigma <- slm$sigma
                LM[[small.models[i]]]$df <- slm$df[2]
                }
            #LM
            }
        else
        if (ata.err.method=="Mack") {
            for (i in which(idx)) {
                a <-LM.MackHeuristic.ATA(LM,i)
                LM[[i]][names(a)] <- a
                }
            }
        else
        if (ata.err.method=="loglinear") {
            # First, cannot extrapolate from < 2 datapoints
            if ((nLM-1-sum(idx))< 2) stop(nLM-1-sum(idx),
                "models are insufficient for loglinear sigma extrapolation")
            # Second, all models must have the same alpha
            ndx <- which(idx)
            if (!vector.all.equal(sapply(LM[ndx], function(x) x$alpha)))
                stop('Alphas must be identical for loglinear sigma method.')
            alpha <- LM[ndx][[1]]$alpha # this is the identical value
            # If all ata models are complete, we're done. Otherwise, ...
            Y<-sapply(LM[ndx], function(x) x$sigma)
            X<-sapply(LM[ndx], function(x) x$ata)
            sigma.lm<-lm(log(Y)~X)
            # predict sigma for all non-tail periods
            ndx <- which(!idx)
            predY.log <- predict(sigma.lm,
                                 newdata=data.frame(X=sapply(LM[ndx],
                                 function(x) x$ata)))
            sigma.pred <- exp(predY.log)
            # save sigma, etc for idx periods
            for (i in ndx) {
                if (!is.null(LM[[i]]$sigma))
                    trowarning("Overwriting sigma at age",LM[[i]]$age)
                LM[[i]]$sigma<-sigma.pred[i]
                sef <- LM[[i]]$sigma /
                            sqrt(sum(LM[[i]]$model[LM[[i]]$xname]^(2-alpha)))
                if (!is.null(LM[[i]]$sef))
                    trowarning("Overwriting sef at age",LM[[i]]$age)
                LM[[i]]$sigma<-sef
                if (!is.null(LM[[i]]$df))
                    trowarning("Overwriting df at age",LM[[i]]$age)
                LM[[i]]$df<-1 # the smallest positive value
                if (!is.null(LM[[i]]$alpha))
                    trowarning("Overwriting alpha at age",LM[[i]]$age)
                LM[[i]]$alpha <- alpha
                attr(LM[[i]],"ata.err.method") <- "loglinear"
                }
            ffmInfoMsg("Implementing loglinear ata error method at age(s)",
                        paste(ndx,collapse=","))
            }
        else stop("incomplete models, no recognized method specified")
        }

    # 4b. Does tail model have sef & sigma?
    #   If not, estimate them using one of the recognized methods.
    if (is.na(LM[[nLM]]$sef) | is.na(LM[[nLM]]$sigma)) {
        if (troUnityTail(LM[[nLM]]$ata)) { # see function in troutil.r
            unity.tail.method <- match.arg(unity.tail.method)
            if (unity.tail.method=="Zero") {
                LM[[nLM]]$sef <- 0
                LM[[nLM]]$sigma <- 0
                LM[[nLM]]$df <- 0
                LM[[nLM]]$alpha <- 0
                }
            else
            if (unity.tail.method=="Mack") { # use ata heuristic on tail
                a <-LM.MackHeuristic.ATA(LM,nLM)
                LM[[nLM]][names(a)] <- a
                }
            else # use loglinear ata method on tail
                ssa <- LM.LogLinear.ATA(LM,nLM)
            }
        else {
            tail.err.method <- match.arg(tail.err.method)
            if (ata.err.method=="Mack")
                ssa <- LM.MackHeuristic.Tail(LM)
            else
                stop(ata.err.method,"not implemented!")
            }
        LM[[nLM]]$sef <- ssa[1]
        LM[[nLM]]$sigma <- ssa[2]
        LM[[nLM]]$df <- ssa[3]
        LM[[nLM]]$alpha <- ssa[4]
        }
    attr(LM,"min.df")<-min.df
    attr(LM,"ata.err.method")<-ata.err.method
    attr(LM,"tail.err.method")<-tail.err.method
    attr(LM,"unity.tail.method")<-unity.tail.method
    attr(LM,"ages") <- sapply(LM,function(x) x$age)
    attr(LM,"endages") <- sapply(LM,function(x) x$endage)

    LM

    }

tdareg.gatherdat <- function(tda, age.begin, age.end, alpha, zero.rm) {
    with(tda, {
    # For an (X,Y) pair to qualify as an observed datapoint, the age of X
    #   must = age.begin and the Y value cannot be NA.
    # We'll allow X=NA for the time being for the cases that
    #   Y is not NA and we want X==0 and X==NA to be equivalent.
    ai <- age==age.begin & !is.na(value.next)
    if (sum(ai)==0) return(data.frame()) # empty d.f -> no obs
    if (zero.rm) { # remove zeroes, NA's
#        zi <- value[ai]!=0 & !is.na(value[ai])
# Change 8/19/11 approximately zero
        zi <- !approx.equal(value[ai], 0) & !is.na(value[ai])
        if (alpha != 0) { # cannot have negative weights in wtd regression
            ni <- value[ai][zi] < 0
            nni <- sum(ni)
            if (nni>0) trowarning("Discarding ",nni," observations < 0 at age ",
                                  age.begin)
            dat <- data.frame(X=value[ai][zi][!ni],
                               Y=value.next[ai][zi][!ni],
                               acp=acp[ai][zi][!ni],
                               age=age.end)#,
#                               evaldt=evaldt[ai][zi][!ni])
            }
        else # regression will be unweighted
            dat <- data.frame(X=value[ai][zi],
                               Y=value.next[ai][zi])
        }
    else { # Don't remove zeroes, treat NA's as zeroes.
        if (alpha!=0) {
            ni <- value[ai]<0
            ni[is.na(ni)] <- FALSE # NAs will become zeroes
            nni<-sum(ni)

            if (nni>0) trowarning("Discarding ",nni," observations < 0 at age ",
                                  age.begin)
            dat  <-data.frame(X=value[ai][!ni],
                               Y=value.next[ai][!ni],
                               acp=acp[ai][!ni],
                               age=age.end
                               )
            # zeroes become a very small positive number otherwise weights
            #   in wtd regression would be infinite
# perhaps could make the small.number = a penny
            dat$X[dat$X==0] <- small.number <- .Machine$double.eps*.5
            dat$X[is.na(dat$X)] <- small.number
            dat$Y[is.na(dat$Y)] <- 0
            }
        else {
            dat <- data.frame(X=value[ai],
                               Y=value.next[ai],
                               acp=acp[ai],
                               age=age.end
#                               , evaldt=evaldt[ai]
                               )
            dat$X[is.na(dat$X)] <- 0
            dat$Y[is.na(dat$Y)] <- 0
            }
        }
    dat
    })
    }

LM.MackHeuristic.ATA <- function(LM,k){
    if (k<3) 
        stop("Not enough prior periods to implement Mack's heuristic at age",k)
    sigma <- sqrt(min(LM[[k-1]]$sigma^4/LM[[k-2]]$sigma^2,
                      min(LM[[k-2]]$sigma^2, LM[[k-1]]$sigma^2))
                  )    
    # That implements Mack's idea for sigma. The formula for sef is:
    #     sef^2 = sigma^2 / denominator
    # where denominator = sum of beginning losses ^ (2-alpha)
    #       in the calculation of this particular ata.
    # Those beginning loss values are the values in LM[[.]]$model in the
    #   column whose name is stored in LM[[.]]$xname
    sef <- sigma/sqrt(sum(LM[[k]]$model[LM[[k]]$xname]^(2-LM[[k]]$alpha)))
    # For alpha and df -- not aspects of Mack's original method --
    #   we will use linear interpolation on the prior two periods, 
    #   where the proportionality factor is based on the sigma's.
    lambda <- (sigma-LM[[k-2]]$sigma)/(LM[[k-1]]$sigma-LM[[k-2]]$sigma)
    df <- LM[[k-2]]$df + lambda * (LM[[k-1]]$df-LM[[k-2]]$df)
    # changed alpha from interpolation to prior alpha 3/26/10
    alpha <- LM[[k-1]]$alpha
#    alpha <- LM[[k-2]]$alpha + lambda * (LM[[k-1]]$alpha-LM[[k-2]]$alpha)
#    ffmInfoMsg(paste("Implementing Mack ata error heuristics at age ",
    # changed from ffmInfoMsg to trowarning because the order of the message
    # showing up on the console was reverse of actual order
    trowarning(paste("Implementing Mack ata error heuristics at age ",
                     LM[[k]]$age, ", alpha = prior alpha = ", alpha, sep=""))
    list(sef=sef,
         sigma=sigma,
         df=df,
         alpha=alpha)
    }

LM.MackHeuristic.Tail <- function(LM,k) {
    ata <- sapply(LM[1:(k-1)],function(x) x$ata)
    age <- sapply(LM[1:(k-1)],function(x) x$age)
    DF <- sapply(LM[1:(k-1)],function(x) x$df)
    # If Tail is larger than any ata, use statistics from largest ata.
    # If there are ties for the largest ata, use the one with the largest
    #   degrees of freedom.
    # If there are ties after that, use last choice.
    if (LM[[k]]$ata > max(ata)) { # tail larger than all ata's
        candidates <- which(LM[[k]]$ata==max(ata))
        maxDF <- max(DF[candidates],na.rm=TRUE)
        if (is.na(maxDF)) stop("all df's NA")
        i <- which(DF[candidates]==maxDF)
        i <- tail(candidates[i],1)
        ffmInfoMsg("Tail",LM[[k]]$ata,
            " > all ata's. Errors proportional to largest ata errors, age=",
            LM[[i]]$age)
        }
    else
    # If Tail is smaller than any ata, use statistics from smallest ata.
    # If there are ties for the smallest ata, use the one with largest df.
    if (LM[[k]]$ata < min(ata)) { # tail smaller than all ata's
        candidates <- which(LM[[k]]$ata==max(ata))
        maxDF <- max(DF[candidates],na.rm=TRUE)
        if (is.na(maxDF)) stop("all df's NA")
        if (maxDF<1) stop("all df's NA")
        i <- which(DF[candidates]==maxDF)
        i <- tail(candidates[i],1)
        trowarning("Tail",LM[[k]]$ata,
            " < all ata's. Errors proportional to smallest ata errors, age=",
            LM[[i]]$age)
        }
    else { # tail bounded by some ata's
        # Find the age (method="constant",rule=1) whose ata (potentially interpolated) = Tail
        # If more than one such age, take the last one ("ties=max")
        g <- approxfun(ata,age,method="constant",rule=1,ties=max)
        i <- g(LM[[k]]$ata) # give g an ata, it returns an age
        if (is.na(i)) stop ("i is na")
        ffmInfoMsg("Implementing Mack tail error heuristics with ata age",
            LM[[i]]$age)
        }
    # Errors will be proportional to the ratio of the excess portion 
    #   of the tail and the factor it was matched with
    lambda <- (LM[[k]]$ata-1)/(LM[[i]]$ata-1)
    c(LM[[i]]$sef*lambda,
      LM[[i]]$sigma*lambda,
      LM[[i]]$df,
      LM[[i]]$alpha)
    }

troMultipleReg <- function(LM,min.df=3) {
    # We've finished looping through ata's.
    # From least to most mature age, loop through LM until find first 
    #     model with degrees of freedom < min.df.
    # Group that age and all more mature ages into one multiple regression.
    # So as to not have to worry about the model that matches the df being 
    #     one for which the Mack heuristic was run, we'll arbitrarily
    #     not allow min.df to be anything less than 3.
    if (min.df<3) stop("min.df < 3 not allowed. Execution terminated.")
  
    # indices of LM's models whose degrees of freedom < min.df
    small.models <- which((sapply(LM,function(x) x$df)<min.df)&
                          (y<-sapply(LM,function(x) x$nobs)>0))
    # If no df's < min.df or only one small model, nothing to do.
    if (length(small.models)<2) return(LM) 
    
    X <- Y <- W <- numeric(0)
    totobs <- sum(sapply(LM[small.models],function(x) x$nobs))
    # all but last model followed by zeroes
    for (i in head(small.models,-1)) {
        # Build a matrix that looks like
        # y x 0 0 0
        # y x 0 0 0
        # y x 0 0 0
        # y x 0 0 0
        # y 0 x 0 0
        # y 0 x 0 0
        # y 0 x 0 0
        # y 0 0 x 0
        # y 0 0 x 0
        # y 0 0 0 x
        # Notice that each set of x's is followed by 10 zeros
        if (LM[[i]]$nobs>0) {
            Y <- c(Y,unlist(LM[[i]]$model[1]))
            X <- c(X,unlist(LM[[i]]$model[2]),rep(0,totobs))
            W <- c(W,unlist(LM[[i]]$model[3]))
            }
        }
    i <- tail(small.models,1)
    Y <- c(Y,unlist(LM[[i]]$model[1]))
    X <- c(X,unlist(LM[[i]]$model[2]))
    W <- c(W,unlist(LM[[i]]$model[3]))
    dim(X) <-c(totobs,length(small.models))
    # Now run the multiple regression
    slm<-summary(lm(Y~X+0,weights=W))
    # Store the results back in the applicable models
    for (i in 1:length(small.models)) {
        LM[[small.models[i]]]$ata <- slm$coef[i,"Estimate"]
        LM[[small.models[i]]]$sef <- slm$coef[i,"Std. Error"]
        LM[[small.models[i]]]$sigma <- slm$sigma
        }
    LM
    }
  
troSummarize.LM <- function(LM) {
#    require(gdata) # for nobs
    .f <- function(x) 
        if (inherits(x,"lm"))
            data.frame(age   = x$age,
                       endage= x$endage,
                       ata   = x$ata,
                       sef   = x$sef, 
                       sigma = x$sigma, 
                       alpha = x$alpha,
                       df    = x$df,
                       formula = paste(as.character(formula(x))[c(2,1,3)],collapse=""),
                       AIC   = AIC(x),
#                       BIC   = AIC(x,k=nobs(x)),
                       BIC   = AIC(x,k=x$nobs),
                stringsAsFactors=FALSE)
        else
            data.frame(age   = x$age,
                       endage= x$endage,
                       ata   = x$ata,
                       sef   = x$sef, 
                       sigma = x$sigma, 
                       alpha = x$alpha,
                       df    = x$df,
                       formula = NA,
                       AIC   = NA,
                       BIC   = NA,
                stringsAsFactors=FALSE)
      
    if (inherits(LM,"list")) {
        DF <- data.frame()
        for (x in LM) DF<-rbind(DF,.f(x))
        }
    else DF<-.f(LM)
    row.names(DF)<-DF$age
    DF
    }

troLM.endage.to.beginage <- function(endage,LM) {
#    ages<-attr(LM,"ages")
#    endages<-attr(LM,"endages")
#    ages[match(endage,endages)]
    attr(LM,"ages")[match(endage,attr(LM,"endages"))]
    }
