troWrapUp <- function(tro,printxl=TRUE) {

    # Summarize cv's and development factors
    tro$summary.byage <- troSummaryByAge(tro$diagonal,tro$LM)

    # Calculate residuals for predicted cumulative losses
    tro$P <- troresid(tro)

    # Send results to excel
    if (printxl) cell <- troPrintToExcel(tro)

    tro

    }

troPlotLogLinearATA <- function(LM) {
    # In the situation where sigma etc for models with df<2 are replaced
    #   by values based on log-linear regression of sigma on age,
    #   display what's going on with two side-by-side graphs, the first
    #   being the log(sigma) values, regression line, and predictions
    #   and the second being the same on the original sigma (y-axis) scale.
    nLM <- length(LM)
    idx <- attr(LM,"loglinear ata err modeled models")
    Y<-sapply(LM[names(idx[!idx])], function(x) x$sigma)
    Y.log <- log(Y)
    X<-sapply(LM[names(idx[!idx])], function(x) x$age)
    sigma.lm<-lm(Y.log~X)
    # predict sigma for all non-tail periods
    predY.log <- predict(sigma.lm,newdata=data.frame(X=sapply(LM[1:(nLM-1)], function(x) x$age)))
    sigma.pred <- exp(predY.log)
    par(mfrow=c(1,2),
        oma = c(0, 0, 3, 0))  # outer margins necessary for page title
    xlim <- c(0,max(attr(LM,"ages")))
    ylim <- range(c(Y.log,predY.log))
    plot(log(Y)~X,xlim=xlim,ylim=ylim,xlab="Age",ylab="log(sigma)")
    abline(sigma.lm,lty="dashed")
    points(as.numeric(names(predY.log[idx])),predY.log[idx],col="red",pch=21,bg="red")
    plot(Y~X,xlim=xlim,ylim=c(0,max(c(Y,sigma.pred))),xlab="Age",ylab="sigma")
    lines(as.numeric(names(predY.log)),sigma.pred,lty="dashed")
    points(as.numeric(names(predY.log[idx])),exp(predY.log[idx]),col="red",pch=21,bg="red")
    # Finally, the overall title of the page of plots    
    mtext("log-linear sigma method for ata's w/ 1 data point", 
            outer = TRUE, cex = 1.5)

    }
