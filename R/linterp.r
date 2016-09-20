linterp <- function(x,y,xvalue,what=c("first","last"),rule=1)
    # x and y are vectors, xvalue is an x value, find the y value 
    sapply(xvalue, 
        function(z) linterp.scalar(x,y,z,what=what,rule=rule))

linterp.scalar <- function(x,y,xvalue,what=c("first","last"),rule=1) {
    what <- match.arg(what) 
    # if xvalue out of bounds, return y's at x endpts (rule=2) or NAs (rule=1).
    # See documentation for 'approx"
    i<-boundingndx(x,xvalue,what=what)
    if (length(i)>1) #bounded
        approx(x[i],y[i],xout=xvalue)$y
    else
    if (rule==1) NA
    else
    if (rule==2)
    if (xvalue<min(x)) y[which.min(x)]
    else y[which.max(x)]
    else stop("Programmer error: rule ",rule," not recognized")
    }

linterp2d.scalar <- function(xval,yval,v,xvec=NULL,yvec=NULL) {
    # x corresponds to a value in the row-direction
    # y corresponds to a value in the column direction
    # v is a matrix of values for which you want to get
    #   an interpolated value for the (x,y) pair.
    # xvec is an ordered vector of x values corresponding to the rows of v
    # yvec is an ordered vector of y values corresponding to the columns of v
    #
    # Decoding flag attribute:
    #   !any(flag): xval and yval in bounds
    #   flag[1]: xval<min x, result assumes xval=min x
    #   flag[2]: xval>max x, result assumes xval=max x
    #   flag[3]: yval<min y, result assumes yval=min y
    #   flag[4]: yval>max y, result assumes yval=max y
    # Ideas outlined at http://wiki.tcl.tk/11339
    if (missing(xvec)) xvec <- attr(v,"rowvec") 
    if (missing(yvec)) yvec <- attr(v,"colvec") 
    Lx <- length(xvec)
    Ly <- length(yvec)
    if (Lx!=nrow(v)) stop("xvec, v mismatch")
    if (Ly!=ncol(v)) stop("yvec, v mismatch")
    # Check x and y against their respective ranges.

    flag <- 0 # assume no problem
    xlbf <- xval<xvec[1];  if (xlbf) {xval<-xvec[1]; flag <- 1}
    xubf <- xval>xvec[Lx]; if (xubf) {xval<-xvec[Lx]; flag <- 2}
    ylbf <- yval<yvec[1];  if (ylbf) {yval<-yvec[1]; flag <- flag + 4}
    yubf <- yval>yvec[Ly]; if (yubf) {yval<-yvec[Ly]; flag <- flag + 8}

#    flag <- logical(4) # assume no problem
#    flag[1]<-xlbf <- xval<xvec[1];  if (xlbf) xval<-xvec[1]
#    flag[2]<-xubf <- xval>xvec[Lx]; if (xubf) xval<-xvec[Lx]
#    flag[3]<-ylbf <- yval<yvec[1];  if (ylbf) yval<-yvec[1]
#    flag[4]<-yubf <- yval>yvec[Ly]; if (yubf) yval<-yvec[Ly]
    # Get the row (x) and column(y) indices bounding xval and yval.
    # Note: boundndx's calculation returns (1,1) for values <= min and
    #       returns (n,n) (n=length(xvec or yvec) for values > max;
    #       otherwise boundndx assumes open-left/closed-right intervals,
    #       returning (i-1,i) where xvec[i-1]<xval<=xvec[i].
    xndx <- boundingndx(xvec,xval)
    yndx <- boundingndx(yvec,yval)
    a <- ifelse(xndx[1]==xndx[2], 0, 
        (v[xndx[2],yndx[1]]-v[xndx[1],yndx[1]])/(xvec[xndx[2]]-xvec[xndx[1]]))
    b <- ifelse(yndx[1]==yndx[2], 0, 
        (v[xndx[1],yndx[2]]-v[xndx[1],yndx[1]])/(yvec[yndx[2]]-yvec[yndx[1]]))
    c <- ifelse(
                xndx[1]==xndx[2] | yndx[1]==yndx[2], 
                0,
                (v[xndx[1],yndx[1]] + v[xndx[2],yndx[2]] - 
                              v[xndx[1],yndx[2]] - v[xndx[2],yndx[1]]) /
                ((xvec[xndx[2]]-xvec[xndx[1]])*(yvec[yndx[2]]-yvec[yndx[1]]))
                )
    structure(
      v[xndx[1],yndx[1]] + a*(xval-xvec[xndx[1]]) + b*(yval-yvec[yndx[1]]) +
                           c*(xval-xvec[xndx[1]])*(yval-yvec[yndx[1]]),
      "flag"=flag)
    }

linterp2d.outofbounds.msg<-function(M,flag,rowname,rowval,colname,colval) {
    msg <-""
    if (flag==0) return(msg)
    rowflag<-flag%%4
    colflag<-flag%/%4
    if (rowflag==1) 
        msg<-paste(rowname, " '", round(rowval,options()$digits), "' < ",
             ifelse(is.null(rownames(M)),"minimum allowable value",
                 paste(rownames(M)[1],"- using", rownames(M)[1],sep=" ")),
             sep="")
    else
    if (rowflag==2) 
        msg <- paste(rowname, " '", round(rowval,options()$digits), "' > ", 
             ifelse(is.null(rownames(M)),"maximum allowable value",
                 paste(rownames(M)[nrow(M)],"- using", rownames(M)[nrow(M)])),
             sep="")
    if (colflag==1) msg <- paste(msg, if (msg!="") ", and ",
                     colname, " '", round(colval,options()$digits), "' < ", 
             ifelse(is.null(colnames(M)),"minimum allowable value",
                 paste(colnames(M)[1],"- using", colnames(M)[1],sep=" ")),
             sep="")
    else
    if (colflag==2) msg <- paste(msg, if (msg!="") ", and ",
                     colname, " '", round(colval,options()$digits), "' > ",
             ifelse(is.null(colnames(M)),"maximum allowable value",
                 paste(colnames(M)[ncol(M)],"- using", colnames(M)[ncol(M)])),
             sep="")
    msg
    }

solve.interp.scalar <- function(x,y,yvalue,what=c("first","last"),rule=1) {
    i<-boundingndx.scalar(y,yvalue,what=what,rule=rule)
    structure(if (is.na(i[1])) NA
        else
            x[i[1]] + if(length(y)>1 & y[i[2]]!=y[i[1]])
                (yvalue-y[i[1]])*(x[i[2]]-x[i[1]])/(y[i[2]]-y[i[1]])
                else 0,
            flag=attr(i,"flag"))
    }

bound<-function(xvec,x,what=c("first","last"),rule=1) {
    a<-boundingndx(xvec,x,what=what,rule=rule)
    structure(xvec[a], dim=if (length(x)>1) c(2,length(x)) else NULL,
                       flag=attr(a,"flag")            )
    }

bound.scalar <- function(xvec,x,what=c("first","last"),rule=1) {
    idx <- boundingndx.scalar(xvec,x,what=what,rule=rule)
    structure(xvec[idx],flag=attr(idx,"flag"))
    }

boundingndx <- function(vec,value,what=c("first","last"),rule=1) {
    # Can't use findInterval because it requires vec be 
    # sorted non-decreasingly.
    a<-lapply(value, 
        function(x) boundingndx.scalar(vec,x,what=what,rule=rule))
    structure(unlist(a),
        dim=if (length(value)>1) c(2,length(value)) else NULL,
        flag=unname(unlist(sapply(a,attributes)))            )
    }

boundingndx.scalar<-function(vec,value,what=c("first","last"),rule=1) {
    # vec is a vector of values
    # Find the pair of indices of vec such that value
    # lies between vec's values, inclusive of endpoints.
    # If pair found:
    #   If value coincides with the first entry in vec, returns c(1,2)
    #   If value coincides with the last entry in vec, returns c(n-1,n)
    #       where n=length(vec).
    #   If value coincides with an interior entry of vec, say vec[i], returns
    #       c(i-1,i) if what="first", c(i,i+1) if what="last".
    # If no such pair exists, returns c(NA,NA) if rule=1. Otherwise (rule=2)
    #   if value > max(vec) returns c(i,i) where i=first/last index of max(vec), 
    #   else returns c(i,i) where i=first/last index of min(vec)
    #   ("rule" logic reflects logic found in 'approx'). 
    # If length(vec)=1, returns 1 if value=vec. If value!=vec returns 
    #   NA if rule=1, returns 1 if rule=2; flag set accordingly.
    # If length(vec)=0, returns NA, flag=3.
    # flag= 0: bound found
    #       1: value < min(vec)
    #       2: value > max(vec)
    #       3: indeterminable because length(vec)=0
    #       4: NA's in vec
    #
    if (NA %in% vec) return(structure(NA,flag=4))
    what <- match.arg(what)
    if (is.na(match(rule,1:2))) 
        stop("Programmer error: rule ",rule," unrecognized.")
    if (length(vec)>1) {
        # Here's all the work. If fails, length(x)=0.
        x<-which(vec[-length(vec)]<=value & value<=vec[-1])
        y<-which(vec[-length(vec)]>=value & value>=vec[-1])
        # Probably most likely case is that it's bounded.
        if (length(x)>0 | length(y)>0) # bounded
            idx <- structure(if (what=="first") c(m<-min(x,y),m+1) 
                             else c(m<-max(x,y),m+1),
                             flag=0)
        else #not bounded
        if (rule==1) idx <- structure(c(NA,NA),flag=ifelse(value<min(vec),1,2))
        else # rule=2; return the first/last index corresponding to min/max
        if (value<min(vec)) 
            idx <- structure(
                if (what=="first") rep(which.min(vec),2) 
                else rep(length(vec)-which.min(rev(vec))+1,2),
                flag=1)
        else 
            idx <- structure(
                if (what=="first") rep(which.max(vec),2)
                else rep(length(vec)-which.max(rev(vec))+1,2),
                flag=2)
        }
    else 
    if (length(vec)>0)
        if (value<vec) 
            idx <- structure(ifelse (rule==1, NA, 1), flag=1)
        else 
        if (value>vec) 
            idx <- structure(ifelse (rule==1, NA, 1), flag=2)
        else 
            idx <- structure(1, flag=0)
    else idx <- structure(NA, flag=3)
    idx
    }

linterp.md <- function () {
    # See also mapprox, discussed at this site:
    # http://www.biostat.wustl.edu/archives/html/s-news/1998-10/msg00034.html
    # Note: 1998
    #
#    > I have a data set allocated on a rectangular grid, say z(x,y) where
#> x_seq(from=xb,to=xe,length=xl);y_seq(from=yb,to=ye,length=yl)
#> Is there a command in S+ to perform interpolation of z into an arbitrary 
#> point 
#> X0,Y0? (Simon Rosenfeld)
#
#Simon--
#
#the routine 'mapprox' below does linear interpolation at a set of m-dimensional 
#points onto an m-dimensional grid of values.
#If the m-dimensional grid has dim( n1, n2, ..., nm), then the points for 
#interpolation are by default assumed to fall in the range
#(1:n1), (1:n2), ..., (1:nm) (per dimension), but you can set the scale of each 
#dimension of the grid using the x.range parameter.

#It's probably wise to experiment with this routine-- I haven't used it directly 
#for quite a while. No guarantees!

#Hope this helps
#Mark Bravington
#m.v.bravington@cefas.co.uk

#Example:

#> # Test interpolation on a fine grid
#> yy <- outer( seq( 0, 1, 0.013), seq( 0, 1, 0.013), 
#               function( x, y) sin( x*y*pi))
#> mapprox(rbind(c(0.2,0.7),c(0.9,0.4),c(0.2,0.3)), yy, 
#           x.range=rbind(c(0,0),c(1,1)))
#$x:
#     [,1] [,2] 
#[1,]  0.2  0.7
#[2,]  0.9  0.4
#[3,]  0.2  0.3

#$y:
#[1] 0.4162331 0.8928420 0.1829592

    # copied from above site, but UNTESTED!
mapprox <- function(x.out,y,x.range=rbind(1,dim(y))) {
  # Linear interpolation of grid values onto points in m dimesions. 
  # Points outside "x.range" are pushed back to side/edge of grid
  # "y" is an m-dimensional array of size n1 x n2 x ... x n.m
  # "x.range" (if supplied) is a matrix of dimension "c(2,m)"
  # "x.out" is a matrix of dimension "c(n.out,m)"

  n.m <- length(dim(y))
  n.out <- nrow(x.out)
  r <- rep(n.out,n.m)

  x <- (x.out-rep(x.range[1,],r)) * 
         rep(dim(y)-1,r) /
         rep(x.range[2,]-x.range[1,],r)
  f.out <- x - floor(x)
  x <- 1 + floor(x)
    
  below <- x < 1; x[below] <- 1; f.out[below] <- 0
  above <- x > rep(dim(y),r); x[above] <- rep(dim(y)-1,r)[above];
      f.out[above] <- 1 

  dx <- as.matrix(do.call("expand.grid",rep(list(c(0,1)),n.m)))
  two.to.n.m <- prod(rep(2,n.m))
  nx <- array(x,c(dim(x),two.to.n.m))
  nx <- nx + aperm(array(dx,c(dim(dx),n.out)),c(3,2,1))
  nx <- matrix(aperm(nx,c(3,1,2)),n.out*two.to.n.m,n.m) 
  y.out <- matrix(y[nx],two.to.n.m,n.out)

  fprod <- matrix(0,n.m*n.out*two.to.n.m,3)
  fprod[,1] <- rep(1:n.out,rep(two.to.n.m,n.out))
  fprod[,2] <-rep(1:n.m,rep(n.out*two.to.n.m,n.m))
  fprod[,3] <-1+dx[rep(1:two.to.n.m,n.out),]

  lff <- array(log(c(1-f.out,f.out)),c(dim(f.out),2))
  fy <- lff[fprod]
  dim(fy) <- c(n.out * two.to.n.m,n.m)
  fy <- fy %*% rep.int(1,n.m)
  fy <- exp(fy)
  dim(fy) <- c(two.to.n.m,n.out)
  
  y.out <- y.out*fy
  y.out <- rep.int(1,two.to.n.m) %*% y.out 

  return(x=x.out,y=as.vector(y.out))
}

    }