troresid <- function(tro) {
    require(MASS)
    # Loop through the list of linear models in tro's LM component,
    #   pulling out the residuals and fitted values and "keyfields" which
    #   as of 2/1/2010 hold acp,age, and evaldt.
    # Index by accident period ('acp') and end age of the development period.
    P <- do.call(rbind,
        lapply(tro$LM,function(x) {
            if (x$nobs>0) {
                if (inherits(x,"lm")) { 
                    if (x$nobs-x$rank-1>0)
                        cbind("residuals"=residuals(x),
                              "standardized.residuals"=stdres(x),
                              "studentized.residuals"=studres(x),
                              "fitted.values"=fitted.values(x),
                              x$keyfields)
                    else
                        cbind("residuals"=residuals(x),
                              "standardized.residuals"=NA,
                              "studentized.residuals"=NA,
                              "fitted.values"=fitted.values(x),
                              x$keyfields)
                    }
                else
                if (is.unity(x$ata)) {
                    y <- merge(x$keyfields,tro$tda)
                    cbind("residuals"=0,
                          "standardized.residuals"=0,
                          "studentized.residuals"=0,
                          "fitted.values"=y$value,
                          x$keyfields)
                    }
                else {
                    y <- merge(x$keyfields,tro$tda)
                    cbind("residuals"=NA,
                          "standardized.residuals"=NA,
                          "studentized.residuals"=NA,
                          "fitted.values"=NA,
                          x$keyfields)
                    }
                }
            })
        )
    rownames(P) <- paste(P$acp,P$age,sep="-")
    P
    }
troresid.homo <- function(tro) {
    r<-troresid(tro)
    r<-cbind(matrix(as.numeric(unlist(strsplit(rownames(r),"-"))),ncol=2,byrow=T,dimnames=list(NULL,c("acp","age"))),r)
    tdf.as.tri(r,what="residuals")
    r<-merge(r,tro$tda,by.x=c("acp","age"),by.y=c("acp","age.next"))
    r<-r[-5]
    r<-merge(r,tro$dpo,by.x="age",by.y="age")
    r$resid.homo<-r$residuals/(r$value^(r$oa/2))
    lag=1
    for (i in (lag+1):7) print(as.numeric(cor.test(trir[(lag+1):(10-i+lag),i-lag],trir[1:(10-i),i])[c("estimate","p.value")]))
    r
    }
troPlotResiduals <- function(tro,main=NULL) {
    # We will plot the residuals as functions of
    #   1. accident period  ("AY dimension")
    #   2. age ("age dimension")
    #   3. evaluation date ("CY dimension")
    #   4. fitted value (observe heteroscedasticity)
    # The y-values of the plots are the P$residuals.
    # The x-values will be the appropriate corresponding fields in
    #   the tro$tda data.frame according to matching rownames.
    par(mfrow=c(3,2),         # 6 plots on one page
        oma = c(0, 0, 5, 0))  # outer margins necessary for page title
    #
    plot(tro$P$acp,
        tro$P$residual,ylab="residuals",
        xlab="Accident Period",
        main="By Accident Period")
    #
    plot(as.factor(tro$P$age),
        tro$P$residual,ylab="residuals",
        xlab=paste("(",attr(tro$tda,"timeunits"),")",sep=""),
        main="By Projected Age")
    # by evaluation date
    if (!is.null(tro$tda$evaldt.next)) {
        par(las=2) # x-axis date labels printed vertically
        plot(as.factor(as.Date(tro$tda[rownames(tro$P),"evaldt.next"])),
            tro$P$residual,ylab="residuals",
            #        xlab="Date",
            main="By Evaluation Date of Projected Loss")
        par(las=0) # back to default
        }
    #
    plot(tro$P$fitted.values,
        tro$P$residual,ylab="residuals",
        xlab="Projected Value",
        main="By Value of Projected Loss")
    # Normality test
    # Only report on statistics for the periods for which "lm" was run, which
    #   are those models for which AIC was run.
    s<-troSummarize.LM(tro$LM)
    P <- tro$P[tro$P$age %in% s[!is.na(s$AIC),"endage"],]
    qqnorm(P$residuals)
    qqline(P$residuals)
    # skip to next plotting frame and display the test's p-value
    plot.new()
    N<-min(length(P$residuals),5000) # 5000=shapiro.test max sample size
    shap.p <- shapiro.test(sample(P$residuals,N))
    text(.5,.9, shap.p$method)
    shap.p.value <- round(shap.p$p.value,5)
    text(.5,.7, paste("p.value", shap.p.value, sep="="))
    text(.5,.5,"At .05 significance level:")
    text(.5,.3, paste(ifelse(shap.p.value<.05,"Should reject",
                                              "Cannot reject"),
                      "the hypothesis"))
    text(.5,.1,"that residuals are normally distributed.")

    # Finally, the overall title of the page of plots
    if (is.null(main))
        mtext("Chain Ladder Residuals (Actual-Projected Values)", 
                outer = TRUE, cex = 1.5)
    else
        mtext(paste(main, "Chain Ladder Residuals (Actual-Projected Values)", sep="\n"),
                outer = TRUE, cex = 1.5)
    par(mfrow=c(1,1))
    }
troPlotStandardizedResiduals <- function(tro, main=NULL) {
    # We will plot the residuals as functions of
    #   1. accident period  ("AY dimension")
    #   2. age ("age dimension")
    #   3. evaluation date ("CY dimension")
    #   4. fitted value (observe heteroscedasticity)
    # The y-values of the plots are the P$residuals.
    # The x-values will be the appropriate corresponding fields in
    #   the tro$tda data.frame according to matching rownames.
    par(mfrow=c(3,2),         # 6 plots on one page
        oma = c(0, 0, 5, 0))  # outer margins necessary for page title
    #
    plot(tro$P$acp,
        tro$P$standardized.residuals,ylab="standardized.residuals",
        xlab="Accident Period",
        main="By Accident Period")
    #
    plot(as.factor(tro$P$age),
        tro$P$standardized.residuals,ylab="standardized.residuals",
        xlab=paste("(",attr(tro$tda,"timeunits"),")",sep=""),
        main="By Projected Age")
    # by evaluation date
    if (!is.null(tro$tda$evaldt.next)) {
        par(las=2) # x-axis date labels printed vertically
        plot(as.factor(as.Date(tro$tda[rownames(tro$P),"evaldt.next"])),
            tro$P$standardized.residuals,ylab="standardized.residuals",
            #        xlab="Date",
            main="By Evaluation Date of Projected Loss")
        par(las=0) # back to default
        }
    #
    plot(tro$P$fitted.values,
        tro$P$standardized.residuals,ylab="standardized.residuals",
        xlab="Projected Value",
        main="By Value of Projected Loss")
    # Normality test
    # Only report on statistics for the periods for which "lm" was run, which
    #   are those models for which AIC was run.
    s<-troSummarize.LM(tro$LM)
    P <- tro$P[tro$P$age %in% s[!is.na(s$AIC),"endage"],]
    qqnorm(P$standardized.residuals)
    qqline(P$standardized.residuals)
    # skip to next plotting frame and display the test's p-value
    plot.new()
    N<-min(length(P$standardized.residuals),5000) # 5000=shapiro.test max sample size
    shap.p <- shapiro.test(sample(P$standardized.residuals,N))
    text(.5,.9, shap.p$method)
    shap.p.value <- round(shap.p$p.value,5)
    text(.5,.7, paste("p.value", shap.p.value, sep="="))
    text(.5,.5,"At .05 significance level:")
    text(.5,.3, paste(ifelse(shap.p.value<.05,"Should reject",
                                              "Cannot reject"),
                      "the hypothesis"))
    text(.5,.1,"that residuals are normally distributed.")

    # Finally, the overall title of the page of plots    
    if (is.null(main))
        mtext("Chain Ladder Residuals (Actual-Projected Values)", 
                outer = TRUE, cex = 1.5)
    else
        mtext(paste(main, "Chain Ladder Standardized Residuals", sep="\n"),
                outer = TRUE, cex = 1.5)
    par(mfrow=c(1,1))
    }


