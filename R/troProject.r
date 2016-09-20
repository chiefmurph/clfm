####################
#
troProjectToUltimate <- function(diagonal,LM,n.terms="3") { 
    if (missing(diagonal) | is.null(diagonal)) stop("diagonal missing")
    if (0==(N<-nrow(diagonal))) return(NULL)
    nLM <- length(LM)
    # projection of detailed data is straightforward
    detail<-project.diagtdf(diagonal,LM,n.terms=n.terms)
#    diagendage<-sapply(diagonal$age, function(x) LMendage[LMendage>x][1])
    # Projection of the sum across acp's is straightforward for 
    #   the expected value (sum of the detail expected value) and for
    #   the process risk (sqrt of the sum of the detail process variance).
    #   Parameter risk is handled by adding together all the diagonal values
    #   that will be projected to the next endage and adding in the prior 
    #   projected value that will also be projected to the same endage. Then
    #   the deltar (parameter risk) formula will act on that beginning value.

    # Create a "development pattern" table 
    dp <- data.frame(begin.age  =sapply(LM,function(x) x$age),
                     end.age    =sapply(LM,function(x) x$endage),
                     observed.begin.value=0,
                     ata        =sapply(LM,function(x) x$ata))

    # Aggregate the projected detail values by endage
    acpsum <- aggregate(data.frame(value=detail$value,gammar2=detail$gammar^2),by=list(age=detail$age),sum)
    # Aggregate the diagonal values by endage = age to which they will
    #   be projected.
    #   The function in sapply finds the first end.age after diagonal's age
    diagonalsum <- aggregate(data.frame(value.x=diagonal$value),by=list(age=sapply(diagonal$age, function(x) dp$end.age[dp$end.age>x][1])),sum)
    # merge them together
    acpsum <- merge(acpsum, diagonalsum, by="age")
    rownames(acpsum)<-acpsum$age
    acpsum$acp<-"SUM"
    # For each age in acpsum, identify the beginning age of the development
    #   pattern 
    acpsum$begin.age <- troLM.endage.to.beginage(acpsum$age,LM)
    # Save the development parameters in acpsum
    acpsum$ata <- sapply(as.character(acpsum$begin.age),function(x) LM[[x]]$ata)
    # For each age in acpsum, identify the beginning value "x" such that
    #   the total in acpsum$value "y" is the projection of "x"
    acpsum$begin.value <- acpsum$value/acpsum$ata
    # The above replaces the three lines below, which are algebraically equivalent
#    acpsum[as.character(y$end.age),"observed.begin.value"]<-y$value
#    acpsum$projected.begin.value = c(0,acpsum$value[-nrow(acpsum)])
#    acpsum$begin.value = acpsum$observed.begin.value+acpsum$projected.begin.value
    # deltar will be filled in by 'predict' below
    acpsum$deltar <- 0
    # but gammar2 is the sum of the detail gammar2's when acp's are independent
    acpsum$gammar<-sqrt(acpsum$gammar2)
#    if (all(on.age.LM(diagonal,LM))) {
#        # See COMMENTS below
#        # "prior" values to be added to next diagonal before projection
#        deltar.prior <- 0
#        for (i in 1:nrow(acpsum)) {
#            newdata<-data.frame(acp="SUM",
#                                age=acpsum$begin.age[i],
#                                value=acpsum$begin.value[i],
#                                deltar=deltar.prior,
#                                gammar=0)
#            deltar.prior<-acpsum$deltar[i] <- 
#                predict(LM[[i]],newdata,LM,n.terms=n.terms)$deltar
#            }
#        }
#    else {
#        # See COMMENTS below
#        acpsum$begin.value <- acpsum$value / dp$ata
    deltar.prior <- 0
    gammar.prior <- 0
    for (i in 1:nrow(acpsum)) {
        newdata<-data.frame(acp="SUM",
                            age=acpsum$begin.age[i],
                            value=acpsum$begin.value[i],
                            deltar=deltar.prior, gammar=gammar.prior)
        x<-predict(LM[[as.character(acpsum$begin.age)[i]]],
                   newdata, LM, n.terms=n.terms)
        # x$deltar assumes development across the "full period" (= default
        #   when all diagonal values are on-age). 
        # When not all diagonal values are on-age, then the implied ata
        #   will be different from the period's ata. In that case, deltar 
        #   will be proportional to the "full-period" value in the same 
        #   proportion as the implied ata bears to the actual ata.
        if (acpsum$value[i]!=0) proportion <- x$value/acpsum$value[i]
        else proportion <- 1
        deltar.prior <- acpsum$deltar[i] <- x$deltar*proportion
#        gammar.prior<-acpsum$gammar[i]
        }

    acpsum$totalr<-sqrt(acpsum$gammar2+acpsum$deltar^2)
    
    list(detail=detail,
         acpsum=acpsum[c("acp","age","value","deltar","gammar","totalr","ata")])
    }
    #########################################################
    #     COMMENTS
    #########################################################
    # Given a list of models, LM, by age, project 'diagonal' to the endage 
    #     of the last model in the list.
    # Properties of models within LM:
    # * endage of a model equals the age of the next model.
    # * ata is the link ratio
    # * sef is the standard error of that link ratio
    # * sigma is the estimate of the residual standard error
    # * alpha is the exponent of the value of loss, X, at the beginning of
    #       the development period such that the variance of the error term
    #       is proportional to X^alpha. This is the FFM assumption.
    #       Combined with sigma, the variance of the error term is
    #       X^alpha*sigma^2
    # All models "should" have an ata field. If there is only one 
    #   observation in the triangle for that development period, 
    #   the model may not have a numerical value for sef or sigma.
    #
    # * each model is of a class for which a 'predict' function exists,
    #       can accept a newdata data frame with components 
    #           acp
    #           age
    #           value
    #           deltar
    #           gammar
    #       and returns a data frame with the same components
    #
    # Value:
    #   A list with a detail section and an acpsum section
    #   The detail section will have an entry for every acp and
    #       every future age.
    #   The acpsum section will have an entry for every future age, 
    #       with projection from the value which is the sum of 
    #       the diagonal value's at that age and the projecte values 
    #       from the prior age.
    #   When all diagonal elements are on-age, that calculation is
    #       straightforward. 
    #   When a diagonal element is off-age, then its first projection 
    #       will be from an interpolated ata and so that diagonal element's
    #       ata, deltar, and gammar estimates will not line up with the
    #       estimates corresponding to the on-age projection from the 
    #       prior age. For the off-age case, we add up the projected values
    #       at each age and "undevelop" them to get the expected diagonal
    #       value at the prior age. We call this a "shortcut" and know that
    #       it only approximates the error estimates for the accident
    #       period total. Note: The Shortcut would yield the same estimates
    #       as the non-Shortcut, on-age algorithm, and they would be exact,
    #       but we do not do that below.
    
    #   Note: Even in the off-age case, gammar for the total can be 
    #       calculated exactly as the square root of the sum of the gammar^2
    #       values because the assumption is that acp's are independent.
    #       When there's only one acp per age, it's not an issue.
    #       Only when there's more than one acp per age are we worried
    #       that the estimate of deltar and gammar for the acp-total will
    #       be off. So we will use the exact value for gammar and adjust
    #       the calculated deltar value by the ratio of the exact gammar
    #       value to the calculated gammar value.
    #       Note also that when using the weighted average link ratios
    #       (alpha=1) the calculated deltar values for the total equal the
    #       exact values (the sqrt of the sum of the squares).
    #

    ########################################################################
    #   ACPSUM SECTION
    #
    # When all diagonal values are on-age, the algorithm would be as 
    #   described in Murphy (PCAS 1994):
    #   1. For each projected age, add up 
    #       a. the predicted values at that age, and 
    #       b. the diagonal values at that age.
    #   1-Shortcut:
    #      Note that we can arrive at the sum of the predicted values and
    #       diagonal values at an age by dividing the sum of the projected
    #       values at the next age by the ata factor that projects
    #       losses to that next age.
    #   2. With that total value and with the calculated value of
    #       deltar and gammar at that age (which start at zero as of
    #       the youngest age as of which no projections have yet been made), 
    #       we predict the value of the acpsum at the next age. 
    #   3A. That predicted sum as of next age will be identical to 
    #       the sum of the predicted detail as of the next age, and the
    #       deltar and gammar values of the predicted sum will be the
    #       appropriate corresponding risk estimates.
    #   3B. Note that the gammar values calculated this way will coincide
    #       with the value you would get simply by summing the gammar^2
    #       detail values by age, and then taking the square root, i.e.,
    #           sqrt(aggregate(detail$gammar^2,list(detail$age),sum)$x)
    #    
    # I suspect that the square of the gammar estimates resulting from 
    #   the Shortcut will equal the sum of the squares of the gammar
    #   estimates of the acp detail (3B above). I've found that to be the 
    #   case with my test, but have yet to prove it holds in all cases.
    #

project.diagtdf<-function(tdf,LM,n.terms="3") {
    # Programmed to avoid repeated use of rbind, which will slow execution
    #   when running model on thousands of detailed data rows.
    # The technique of creating a list of results of calculations, then
    #   using 'do.call(rbind,...)' at the end is a well-documented
    #   technique in R. 
    #   The problem with successively building 'detail' within
    #   the loop as detail <- rbind(detail,predict(...))
    #   is that each rbind copies to a new 'detail' data frame then destroys
    #   the old 'detail' data frame, resulting in slow execution, especially 
    #   after many iterations with larger and larger 'detail' data frames.
    detail <- structure(vector("list",nLM<-length(LM)),names=names(LM))
    agec.prior <- '' # age of prior projections
    # vector of endage's for each entry in LM list
    LMendage <- sapply(LM,function(x) x$endage)
    # First LM index where tdf's age < LM's endage.
    agendx <- sapply(tdf$age, function(x) which(x<LMendage)[1]) 
    for (k in 1:nLM) {
        i <- agendx==k
        # Possible for not to have a diag value to develop at this age
        if (sum(i)>0) newdata<-rbind(data.frame(tdf[i,1:3],
                                                deltar=0,gammar=0,totalr=0),
                                     detail[[agec.prior]]) #=NULL 1st loop iter
        else newdata <- data.frame(detail[[agec.prior]])
#if (k==18){
#prn(k)
#prn(head(i))
#prn(sum(i))
#prn(head(detail[[agec.prior]]))
#prn(nrow(newdata))
#prn(head(newdata))
#}
        if (nrow(newdata)>0)
            detail[[k]] <- predict(LM[[k]],newdata,LM,n.terms=n.terms)
        agec.prior <- as.character(LM[[k]]$age)
        }
    do.call(rbind,detail)
    }

on.age.LM <- function(tdf,LM) 
    sapply(tdf$age,function(x) x %in% sapply(LM,function(x) x$age))
