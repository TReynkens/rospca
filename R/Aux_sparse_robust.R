# Auxiliary functions for ROBPCA and ROSPCA; function for diagnostic plot.

#Compile functions to increase speed, included in R base for versions >= R 2.13.0
setCompilerOptions("suppressAll",TRUE)
enableJIT(1) #just-in-time (JIT) compilation

#################################################################################
#Auxiliary functions for ROBPCA and ROSPCA
#################################################################################

#univariate MCD estimator based on Matlab code from LIBRA:
#https://wis.kuleuven.be/stat/robust/LIBRA
#x is the data vector and h the size of the subsets
unimcd <- function (x, h) {
  nobs <- length(x)
  len <- nobs-h+1
  xorig <- x
  
  if (len==1) {
    tmcd <- mean(x)
    smcd <- sd(x)
    weights <- rep(1,nobs)
    
  } else {

    #Sort data and keep indices
    sorted <- sort(x,index.return=TRUE)
    x <- sorted$x
    I <- sorted$ix
    
    
    #####
    #Look for h-subset with smallest variance
    #We do this by looking at contiguous h-subsets which allows us 
    #to compute the variances (and means) recursively
    
    #Sum of h-subsets
    ax <- numeric(len)
    ax[1] <- sum(x[1:h])
    for (samp in 2:len) {
      ax[samp] <- ax[samp-1]-x[samp-1]+x[samp+h-1]
    }
    
    #Sample variances of h-subsets multiplied by h-1
    ax2 <- ax^2/h
    sq <- numeric(len)
    sq[1] <- sum(x[1:h]^2)-ax2[1]
    for (samp in 2:len) {
      sq[samp] <- sq[samp-1]-x[samp-1]^2+x[samp+h-1]^2-ax2[samp]+ax2[samp-1];
    }
    
    #Subset(s) with minimal variance
    sqmin <- min(sq)
    ii <- which(sq==sqmin)
    ndup <- length(ii)
    optSets <- ax[ii]
    
    #initial mean is mean of h-subset
    initmean <- optSets[floor((ndup+1)/2)]/h 
    #initial variance is variance of h-subset
    initcov <- sqmin/(h-1) 
    
    #####
    # Consistency factor
    
    res <- (xorig-initmean)^2/initcov
    sortres <- sort(res)
    factor <- sortres[h]/qchisq(h/nobs,1)
    initcov <- factor*initcov
    
    #####
    #Weights and reweighting
    
    #raw squared robust distances
    res <- (xorig-initmean)^2/initcov 
    quantile <- qchisq(0.975,1)
    
    #raw weights
    weights <- (res<=quantile) 
    
    #reweighting procedure
    tmcd <- sum(xorig*weights)/sum(weights)
    smcd <- sqrt(sum((xorig-tmcd)^2*weights)/(sum(weights)-1))
  }
  
  return(list(tmcd=tmcd,smcd=smcd,weights=weights))

}




#Standardise the matrix X using the location estimator f_m
#and the scale estimator f_s (both are applied to the columns of X).
#In the default case the data are standardised classically.
#Can handle data matrices!!!!
matStand <- function (X, f_c = mean, f_s = sd, center = TRUE, scale = TRUE) {
  
  p <- ncol(X)
  
  if (center) {
    center <- apply(X,2,f_c)
  } else {
    center <- rep(0,p)
  }
  
  if (scale) {
    scale <- apply(X,2,f_s)
    if (any(scale<10^(-16))) {
      stop("The scale measure of at least one variable is 0.")
    }
  } else {
    scale <- rep(1,p)
  }
  
  data <- sweep(X,2,center,"-")
  
  return(list(data=sweep(data,2,scale,"/"),c=center,s=scale))
  
}

#Standardise the matrix X using the vectors center and scale
matStandVec <- function (X, center, scale) {
  return(sweep(sweep(X,2,center,"-"),2,scale,"/"))
}

#Better to use sweep(X,2,x,"-") instead of X-repvec(x,nrow(X))
# #Vector x will be repeated n times (each row of the returned matrix is the vector x)
# repvec <- function (x, n){
#   sapply(x,FUN=rep,times=n)
# }

#Norm of a vector
vecnorm <- function (x) {
  
  return(sqrt(sum(x^2)))
}


#Cutoff for the ODs, also returns ODs
coOD <- function (od, h = length(od), skew = FALSE) {
  
  #Problems can occur when the ODs are very small
  if (all(od<10^(-14))) {
    co.od <- 0; od <- rep(0,length(od))
  } else {
    if (skew) {
      mcod <- medcouple(od)
      if (mcod>0) {
        co <- quantile(od,0.75) + 1.5*exp(3*mcod)*IQR(od) 
      } else {
        co <- quantile(od,0.75) + 1.5*IQR(od) 
      } 
      #cutoff is largest odh smaller than co
      co.od <- max(od[which(od<=co)])
    } else {
      mcd <- unimcd(od^(2/3),h=h)
      co.od <- (mcd$tmcd+mcd$smcd*qnorm(0.975))^(3/2)
    }
  }
  return(list(od=od,co.od=co.od,wieghts=weights))
}

# Function to compute score distances and orthogonal distances, it also returns the cutoff values.
# X is the data matrix; Tn and P are the scores and loadings matrix of X using a certain PCA method.
# l contains the eigenvalues and mu is an estimate for the centre.
# h is an option for unimcd (size of subsets).
# skew indicates if cutoffs for skew distributions are used
distPCA <- function (X, Tn, P, l, mu = rep(0,ncol(as.matrix(X))), h = NULL, skew = FALSE, co.od = TRUE) {
  X <- as.matrix(X)
  Tn <- as.matrix(Tn)
  d <- dim(Tn)
  n <- d[1]; k <- d[2]
  if (is.null(h)) h <- min(ceiling(n*0.75)+1,n)

  mu <- as.vector(mu) #Problems can occur with dimensions
  
  #Score distance, scores are mean-centred
  if (any(l<10^(-14))) {
    #Zero variation in at least one PCA direction
    #Neglect PCA direction(s) with zero variation
    ix <- which(l<10^(-14))
    k2 <- k-length(ix)
    #Specify size of diagonal matrix to avoid problems when k2=1
    sd <- sqrt(mahalanobis(Tn[,-ix],rep(0,k2), diag(1/l[-ix],k2,k2), inverted=TRUE))
  } else {
    #Specify size of diagonal matrix to avoid problems when k=1
    sd <- sqrt(mahalanobis(Tn,rep(0,k), diag(1/l,k,k), inverted=TRUE))
    k2 <- k
  }
  #Cutoff value
  if (skew) {
    #Cutoff is defined similarly as cutoff for ODs
    co.sd <- coOD(sd,skew=TRUE)$co.od
  }else{
    co.sd <- sqrt(qchisq(0.975,df=k2))
  }
  
  #Orthogonal distance
  od <- apply(sweep(X,2,mu,"-")-Tn%*%t(P),1,vecnorm)
  
  #Cutoff
   if(co.od) {
    tmp <- coOD(od=od,h=h,skew=skew)
  } else {
    tmp <- list()
    tmp$od <- od
    tmp$co.od <- NULL
  }
  
  return(list(sd=sd,od=tmp$od,cutoff.sd=co.sd,cutoff.od=tmp$co.od))
}


#################################################################################
#Diagnostic plot
#################################################################################

#Label outliers on diagnostic plot
#Code by Valentin Todorov for rrcov package
labelDD <- function (x, y, id.n.sd = 3, id.n.od = 3, off = 0.02) {
  xrange <- par("usr")
  xrange <- xrange[2] - xrange[1]
  if (id.n.sd > 0 && id.n.od > 0) {
    n <- length(x)
    ind.sd <- sort(x, index.return = TRUE)$ix
    ind.sd <- ind.sd[(n - id.n.sd + 1):n]
    ind.od <- sort(y, index.return = TRUE)$ix
    ind.od <- ind.od[(n - id.n.od + 1):n]
    lab <- ind.od
    if (is.character(names(y))) {
      lab <- names(y[ind.od])
    }
    text(x[ind.od] - off * xrange, y[ind.od], lab)
    lab <- ind.sd
    if (is.character(names(x))) {
      lab <- names(x[ind.sd])   
    }
    text(x[ind.sd] - off * xrange, y[ind.sd], lab)
  }
}

#Make a diagnostic plot.
#If label_out=TRUE, the outliers will be labelled on the plot.
diagPlot <- function (res,title="Robust PCA",col="black",pch=16,labelOut=TRUE,id=3){
  
  sd <- res$sd
  od <- res$od
  co.sd <- res$cutoff.sd
  co.od <- res$cutoff.od
  

  if (is.null(pch)) {
    pch <- ifelse(is.null(col),1,16)
  }
  if (is.null(col)) {
    col <- "black"
  }
 
  plot(sd,od,xlab="Score distance",ylab="Orthogonal distance",main=title,pch=pch,col=col,
       xlim=c(0,max(sd)*1.1),ylim=c(0,max(od)*1.1))


  abline(v=co.sd)
  abline(h=co.od)
  
  #label outliers on plot
  if (labelOut) {
    labelDD(sd,od,id.n.sd=id,id.n.od=id)
  }     
}
