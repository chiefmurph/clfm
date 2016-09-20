####################
# PSI.r
# Makes the PSI table used in troProject under the FFM model.
#

# Simply source this file.
# The four commands at the bottom call the TNP.build function
# and the PSI.build function and save the tables to disk.

# The purpose of the TNP.build function is to find the (mu,std) pairs where, 
#   for each pair of parameters, the positive normal random variable 
#   with those parameters will have mean mu=1 and cv in cvseq.
# Those (mu,std) pairs (which are the same as the (1,cv) pairs) are saved
#   to disk.
#
# The purpose of the PSI.build function is to create the PSI table.
#   The PSI table is a matrix whose rows are indexed by alpha, a sequence
#   of real numbers between -4 and 8 at .25 increments, whose columns are
#   indexed by cv, a sequence of numbers between 0 and 2 at .025 increments,
#   and whose value is the ratio of E(X^alpha)/E(X)^alpha for a
#   positive normal random variable X with mean 1 and the indexed cv, the 
#   (mu,sigma) parameters of which were previously calculated by TNP.build.


# Not used. R-Devel-oper suggestion for how to generate a vector x long of 
# samples of normal random variables truncated below and above by a and b.
#rnorm.trunc <- function(x,mu=0,std=1,a=-Inf,b=Inf) {
#    n <- length(x)
#    if (n==1) n <- x
#    totsamp <- rnorm(n,mu,std)
#    outside <- totsamp<a | totsamp>b
#    while (any(outside)) {
#        totsamp[outside] <- rnorm(sum(outside),mu,std)
#        outside <- totsamp<a | totsamp>b
#        }
#    totsamp
#    }

# Generate positive normal random variables.
# Contrasted with the above function, there is no guarantee that
#   the truncated sample will have length N.
# This version truncates the sample of N normal random variables at zero
#   but does not truncate above.
# It also truncates below at a point where the inverse would be more
#   than 'cutoff' standard deviations above the mean. Default cutoff=3.72
#   corresponds to the 99.99%-ile of the normal distribution.
# It is expected that the function will be called with positive mu; in fact,
#   most lkely called with mu=1 and std=coefficient of variation.
rnorm.pos <- function(mu,std,N=10000,cutoff=3.72) {
    if (std<0) return(NA)
    if (std==0) return(rep(mu,N))
    y<-rnorm(N,mu,std)
    # Perform rejection sampling. See discussion at
    # http://www.biostat.wustl.edu/archives/html/s-news/2001-04/msg00033.html
    y<-y[y>0]
    # Also omit small values whose reciprocals will be "overly large" where
    # "overly large" means more than cutoff std's above mu.
    if (is.finite(cutoff)) {
        radius <- cutoff*std
        if (mu-radius<=0) {
            inverse.too.large <- y < 1/(mu+radius)
            y<-y[!inverse.too.large]
            }
        }
    y
    }

# This function is used by optim to determine how close the
# calculated mu (parms[1]) and sigma (parms[2]) are to the
# targetmu and targetsd, based on the sum of the absolute differences.
abs.diffparms <- function(parms,targetmu,targetsd,N=10000) {
    mu <- parms[1]
    std <- parms[2]
    y <- rnorm.pos(mu,std,N=N)
    if (length(y)==0) {prn(mu);prn(std);stop("y empty")}
    else
    if (is.na(y[1])) Inf
    else
    abs(mean(y)-targetmu)+abs(sd(y)-targetsd)
    }
    
norm.pos.musigma <- function(targetsd,targetmu=1,N=10000) {
    # Given targetmu and targetsd (probably the cv)
    # find the values of mu and sd such that a truncated normal distribution
    # (truncated at zero, i.e., must be positive) has mean = target mu
    # and standard deviation = target sd.
    S<-optim(c(targetmu,targetsd),
             abs.diffparms, # Function to be minimized (above)
             gr=NULL,       # No gradient function provided
             targetmu=targetmu,    # Other arguments ...
             targetsd=targetsd, # ... to be passed ...
             N=N)           # ... to abs.diffparms function.
#prn(length(rnorm.pos(S$par[1],S$par[2],N=N))/N)
    S$par
    }

norm.pos.musigmaconstr <- function(targetsd,targetmu=1,N=10000) {
    # Not used.
    # Same as above, but with constrained optimization.
    # Unconstrained version above is about twice as fast.
    a<-Sys.time()
    S<-constrOptim(c(targetmu,targetsd),
                abs.diffparms,
                gr=NULL,
                ui=matrix(c(1,0,0,1),2,2), # ui %*% c(mu,sd)>= ci=c(0,0)
                ci=c(0,0),                 # mu, sd must be positive
                targetmu=targetmu,
                targetsd=targetsd,
                N=N)
    b<-Sys.time()
    print(b-a)
    S$par
    }

# Function to find the (mu,std) pairs where, for each pair of parameters,
#   the positive normal random variable with those parameters will
#   have mean mu-1 and cv in cvseq.

# 3/17/09: Final builds of TNP and PSI were run as follows:
# > TNP <- TNP.build(N=1000000)
# > save(TNP,file='TNP.RData')
# > PSI <- PSI.build() # no arguments -> inputs taken from TNP.RData
# > save(PSI,file='PSI.Rdata')
TNP.build <- function(cvseq=seq(0,0.5,by=.025),N=10000) {
    set.seed(123)
    a<-Sys.time()
    cat("Starting at ",format(a),"\n")
    flush.console()
    TNP<-matrix(nrow=2,ncol=length(cvseq),
        dimnames=list(c("mu","sd"),formatC(cvseq,digits=3,format="f")))
    attr(TNP,"cvseq")<-cvseq
    attr(TNP,"N")<-N
    for (i in 1:length(cvseq)) {
        TNP[,i] <- norm.pos.musigma(cvseq[i],targetmu=1,N=N)
        save(TNP,file='TNP.RData')
        cat("cv= ",cvseq[i],", elapsed time=",format(Sys.time()-a),"\n")
        flush.console()
        }
    cat("Ending at ",format(Sys.time()),"\n")
    flush.console()
    TNP
    }

# Function to build the PSI matrix used in the calculation
#   of process risk under the FFM model. A PSI value is the
#   ratio E(X^n)/[E(X)]^n.
PSI.build <- function(TNP=NULL,alphaseq=seq(-4,8,by=.25)) {
    # ffm(p.16): formula=E(x^a)*sigma2+f^2*gammar2
    # E(X^n) = (E^n(X))*PSI(alpha,cv). See
    #     http://en.wikipedia.org/wiki/Normal_distribution#Moments
    # For non-integer alpha, we do linear interpolation on the 
    #   integer moments.

    # alphaseq -> row
    # cvseq -> col

    a<-Sys.time()
    cat("Starting at ",format(a),"\n")
    flush.console()
    set.seed(31609)
    if (missing(TNP)) load('TNP.RData')
    cvseq <- attr(TNP,"cvseq")
    N <- attr(TNP,"N")
    prn(cvseq)
    prn(N)
    PSI<-matrix(nrow=length(alphaseq),ncol=length(cvseq),
                dimnames=list(formatC(alphaseq,digits=3,format="f"),
                              formatC(cvseq,digits=3,format="f")    ))
    attr(PSI,"alphaseq")<-alphaseq
    attr(PSI,"cvseq")<-cvseq
    # Two attributes used for lookups by linterp2d
    attr(PSI,"rowvec") <- alphaseq
    attr(PSI,"colvec") <- cvseq
    attr(PSI,"N")<-N
    for (i in 1:length(cvseq)) {
        #mu <- TNP["mu",i]
        #std <- TNP["sd",i]
        # When cvseq!=attr(TNP,"cvseq"), the indices of cvseq and the
        # columns of TNP do not match.
        # Therefore, above 'i' needs to be replaced by finding the index
        # of cvseq in attr(TNP,"cvseq"). BUT .3 does not match!
        # I.e., the command .3 %in% seq(0,.8,by=.025) returns FALSE
        # Even just .3 %in% seq(.2,.3) returns FALSE
        # WORKAROUND: Match the column names.
        std<-formatC(cvseq[i],digits=3,format="f")
        mu <- TNP["mu",std]
        std <- TNP["sd",std]
        y <- rnorm.pos(mu,std,N)
        negvals <- sum(y<0)
        # 3.72 = qnorm(.9999), the 1-in-10,000 year event
        smallvals <- sum(y<1/(mu+3.72*std))
        y<-y[y>0] 
        mny <- mean(y)
        y <- y/mean(y)
        sdy <- sd(y)
        # y scaled to mean=1, so (mean(y))^a = 1, so omitted in next line
        PSI[,i] <- sapply(alphaseq, function(a) mean(y^a))
#        save(PSI,file='PSI.RData')
        # Hit Ctrl-W to toggle console to update when written to 
        # (default is only when input required)
        cat("cv= ",      cvseq[i],
            "TNP col.ndx=",match(std,colnames(TNP)),
            ", negvals=",negvals,
            ", smallvals=",smallvals,
            ", mny=",    mny,
            ", sdy=",    sdy,
            ", elapsed time=",format(Sys.time()-a),"\n")
        flush.console()
        }
    cat("Ending at ",format(Sys.time()),"\n")
    flush.console()
    PSI
    
    }
 TNP <- TNP.build(cvseq=seq(0,2,by=.025),N=100000)
 save(TNP,file='TNP.RData')
 PSI <- PSI.build() # no arguments -> inputs taken from TNP.RData
 save(PSI,file='PSI.Rdata')
