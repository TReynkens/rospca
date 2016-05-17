

robpca <- function (x, k = 0, kmax = 10, alpha = 0.75, h = NULL, mcd = FALSE, ndir = "all", skew = FALSE, ...)
{

  ## k    -   Number of principal components to compute. If \code{k} is missing,
  ##              or \code{k = 0}, the algorithm itself will determine the number of
  ##              components by finding such \code{k} that \code{l_k/l_1 >= 10.E-3} and
  ##              \code{sum_1:k l_j/sum_1:r l_j >= 0.8}. It is preferable to
  ##              investigate the scree plot in order to choose the number
  ##              of components and the run again. Default is \code{k=0}.
  ##
  ## kmax -   Maximal number of principal components to compute Default is \code{kmax=10}.
  ##              If \code{k} is provided, \code{kmax} does not need to be specified,
  ##              unless \code{k} is larger than 10.
  ## alpha    This parameter measures the fraction of outliers the algorithm should
  ##              resist. In MCD alpha controls the size of the subsets over which the determinant
  ##              is minimized, i.e. alpha*n observations are used for computing the determinant.
  ##              Allowed values are between 0.5 and 1 and the default is 0.5.
  ## h    -  (n-h+1) measures the number of outliers the algorithm should 
  ##               resist. Any value between n/2 and n may be specified. (default = NULL)
  ##               Alpha and h may not both be specified.
  ## mcd  -   Logical - when the number of variables is sufficiently small,
  ##              the loadings are computed as the eigenvectors of the MCD covariance matrix,
  ##              hence the function \code{\link{CovMcd}()} is automatically called. The number of
  ##              principal components is then taken as k = rank(x). Default is \code{mcd=TRUE}.
  ##              If \code{mcd=FALSE}, the ROBPCA algorithm is always applied.
  ## ndir:        Number of directions used when computing outlyingness (or adjusted outlyingness when skew=TRUE). 
  ##              When "all" (default), we use n choose 2 directions which are all directions through 2 data points.
  ##              Giving ndir > n choose 2 as an input also results in using all directions through 2 data points.  
  ## trace    whether to print intermediate results. Default is \code{trace = FALSE}}
  ## skew indicates if the version for skew data is applied
 
  ##  Example:
  ##  data(hbk)
  ##  pca <- robpca(hbk)
  
  ## We return a list with following values
  ##
  ## loadings     -   Sparse robust loadings (eigenvectors)
  ## eigenvalues  -   Robust eigenvalues
  ## scores       -   Scores matrix
  ## center       -   Centre of the data
  ## k            -   Number of (chosen) principal components
  ## h            -   The quantile h used throughout the algorithm
  ## alpha        -   The robustness parameter alpha used throughout the algorithm
  ## n.obs        -   Number of observations
  ## sd           -   Robust score distances within the robust PCA subspace
  ## od           -   Orthogonal distances to the robust PCA subspace
  ## cutoff.od    -   Cutoff value for the robust score distance
  ## cutoff.sd    -   Cutoff value for the orthogonal distance
  ## flag.sd      -   The observations whose score distance is larger than cutoff.sd
  ##                  receive a flag.sd equal to zero.
  ## flag.sd      -   The observations whose orthogonal distance is larger than cutoff.sd
  ##                  can be considered as outliers and receive a flag equal to zero.
  ## flag.all     -   The observations whose score distance is larger than cutoff.sd
  ##                  or whose orthogonal distance is larger than cutoff.od
  ##                  can be considered as outliers and receive a flag equal to zero.
  ##                  The regular observations receive flag 1.
  ## H0           -   Indices of the h least outlying points  
  ## H1           -   1 for observations that are kept after reweighting step
  # Based on code for ROBPCA in rrcov package (PcaHubert function).
  # Author: Tom Reynkens (tomreynkens@hotmail.com)
  
  #The number of subsets used in covMCD (MCD version) is set to 1000.
  #Or use all directions in outlyingness (using mrfDepth) when no mcd is applied.
  #We also changed the computation of distances and flags and allow for an input element h.
  #Using the skew option, one can also use the skewed version of ROBPCA.
  
  
  scale <- FALSE; signflip <- TRUE; trace <- FALSE
  
  cl <- match.call()
  
  if(missing(x))
    stop("You have to provide at least some data")
  
  data <- as.matrix(x)
  n <- nrow(data)
  p <- ncol(data)

  ## ___Step 1___: Reduce the data space to the affine subspace spanned by the n observations
  ##  Apply svd() to the mean-centered data matrix. If n > p we use the kernel approach -
  ##  the decomposition is based on computing the eignevalues and eigenvectors of(X-m)(X-m)'
  #Xsvd <- kernelEVD(data, scale=scale, signflip=signflip)
  #Xsvd <- classPC(data, scale=scale, signflip=signflip)
  Xsvd <- classPC(data, scale=scale, signflip=signflip)
  
  if(Xsvd$rank == 0)
    stop("All data points collapse!")

  ## VT::27.08.2010: introduce 'scale' parameter; return the scale in the value object
  ##
  myscale  <-  vector('numeric', p) + 1
  if(scale)
    myscale <- Xsvd$scale
  
  ##
  ## verify and set the input parameters: alpha, k and kmax
  ## determine h based on alpha and kmax, using the function h.alpha.n()
  ##

  ## VT::06.11.2012 - kmax <= floor(n/2) is too restrictive
  ##    kmax <- max(min(floor(kmax), floor(n/2), Xsvd$rank),1)
  kmax <- max(min(floor(kmax), Xsvd$rank),1)
  
  
  if((k <- floor(k)) < 0)
    k <- 0
  else if(k > kmax) {
    warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
    k <- kmax
  }
  
  if(missing(alpha))
  {
    default.alpha <- alpha
    if(is.null(h)){
      h <- min(h.alpha.n(alpha, n, kmax), n)
    }

    alpha <- h/n
    if(k == 0) {
      if(h < floor((n+kmax+1)/2)) {
        h <- floor((n+kmax+1)/2)
        alpha <- h/n
        warning(paste("h should be larger than (n+kmax+1)/2. It is set to its minimum value ", h, ".", sep=""))
      }
    }
    else {
      
      if(h < floor((n+k+1)/2)) {
        h <- floor((n+k+1)/2)
        alpha <- h/n
        warning(paste("h should be larger than (n+k+1)/2. It is set to its minimum value ", h, ".", sep=""))
      }
    }
    if(h > n) {
      
      alpha <- default.alpha
      if(k == 0)
        h <- h.alpha.n(alpha, n, kmax)
      else
        h <- h.alpha.n(alpha, n, k)
      warning(paste("h should be smaller than n = ", n, ". It is set to its default value ", h, ".", sep=""))
    }
  }else
  {
    
    if(alpha < 0.5 | alpha > 1)
      stop("Alpha is out of range: should be between 1/2 and 1")

    hdef=h
    if(k == 0)
      h <- h.alpha.n(alpha, n, kmax)
    else
      h <- h.alpha.n(alpha, n, k)
    
    if(!is.null(hdef))
      warning(paste("Both alpha and h are given, h is set to its default value ",h,".",sep=""))
  }
  
#   X <- Xsvd$scores
  X <- scale(data, center=Xsvd$center, scale=Xsvd$scale) %*% Xsvd$loadings
  center <- Xsvd$center
  rot <- Xsvd$loadings
  
  ##
  ## ___Step 2___: Either calculate the standard PCA on the MCD covariance matrix (p<<n)
  ##  or apply the ROBPCA algorithm. If mcd=FALSE or skew=TRUE, always apply ROBPCA.
  ##
  if (mcd & skew) {
    warning("The skew-adjusted version of ROBPCA cannot be applied with the MCD adaptation, \"mcd\" is set to FALSE.")
  }
  if(ncol(X) <= min(floor(n/5), kmax) & mcd & !skew)    # p << n => apply MCD
  {
    if(trace)
      cat("\nApplying MCD.\n")
    
    ## If k was not specified, set it equal to the number of columns in X
    ##
    if(k != 0)
      k <- min(k, ncol(X))
    else {
      k <- ncol(X)
      if(trace)
        cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")
    }
    
    X.mcd <- CovMcd(as.data.frame(X),alpha=alpha,nsamp=1000,seed=set.seed(251120134))
    X.mcd.svd <- svd(getCov(X.mcd))
    #scores <- (X - repmat(getCenter(X.mcd), nrow(X), 1)) %*% X.mcd.svd$u
    scores <- (X - matrix(rep(getCenter(X.mcd), times=nrow(X)), nrow=nrow(X), byrow=TRUE)) %*% X.mcd.svd$u
    
    center <- as.vector(center + getCenter(X.mcd) %*% t(rot))
    H0 <- X.mcd@best
    
    eigenvalues <- X.mcd.svd$d[1:k]
    loadings <- Xsvd$loadings %*% X.mcd.svd$u[,1:k]
    scores <- as.matrix(scores[,1:k])
    if(is.list(dimnames(data)) && !is.null(dimnames(data)[[1]]))
    {
      dimnames(scores)[[1]] <- dimnames(data)[[1]]
    } else {
      dimnames(scores)[[1]] <- 1:n
    }
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))
    dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))
    
    # Dummy indexset to avoid errors
    indexset <- H0
    
  }
  else                                        # p > n or mcd=FALSE => apply the ROBPCA algorithm
  {
    if(trace)
      cat("\nApplying the projection method of Hubert.\n")
    
  
    if (ndir=="all" & n>500) {
      warning("Computing all directions for this value of n can take very long, we
            recommend to set ndir to 5000.")
    }
    
    if (skew) {
      # Adjusted outlyingness from mrfDepth package, fast enough to compute all directions
      outl_obj <- adjOutlyingness(X, options = list(type="Rotation", ndir=ndir))

    } else {
      # Outlyingness from mrfDepth package, fast enough to compute all directions
      outl_obj <- outlyingness(X, options = list(type="Rotation", h=h, ndir=ndir, stand="unimcd"))
    }

    outl <- outl_obj$outlyingnessX
    if(!is.numeric(outl)) {
      stop("Something went wrong with the outlyingness computations.")
    }
    
    if(any(!is.finite(outl)) | any(is.nan(outl))) {
      stop("The obtained outlyingness is not finite or NaN.")
    }
    
    if(min(outl)<10^(-16)) {
      stop("The obtained outlyingness is not strictly positive.")
    }
    
    if(length(outl)!=nrow(X)) {
      stop("The obtained outlyingness does not have the correct length.")
    }
    
    is <- order(outl)  #Order indices by outlyingness   
    H0 <- (1:n) %in% is[1:h]
    Xh <- X[H0, ]                      # the h data points with smallest outlyingness
    #Xh.svd <- classSVD(Xh)
    Xh.svd <- classPC(Xh)
    
    kmax <- min(Xh.svd$rank, kmax)
    if(trace)
      cat("\nEigenvalues: ", Xh.svd$eigenvalues, "\n")
    
    ##
    ## Find the number of PC 'k'
    ## Use the test l_k/l_1 >= 10.E-3, i.e. the ratio of
    ## the k-th eigenvalue to the first eigenvalue (sorted decreasingly) is larger than
    ## 10.E/3 and the fraction of the cumulative dispersion is larger or equal 80%
    ##
    if(k == 0)
    {
      test <- which(Xh.svd$eigenvalues/Xh.svd$eigenvalues[1] <= 1.E-3)
      k <- if(length(test) != 0)  min(min(Xh.svd$rank, test[1]), kmax)
      else                   min(Xh.svd$rank, kmax)
      
      cumulative <- cumsum(Xh.svd$eigenvalues[1:k])/sum(Xh.svd$eigenvalues)
      if(cumulative[k] > 0.8) {
        k <- which(cumulative >= 0.8)[1]
      }
      if(trace)
        cat(paste("The number of principal components is set by the algorithm. It is set to ", k, ".\n", sep=""))
    }
    
    if(trace)
      cat("\nXsvd$rank, Xh.svd$rank, k and kmax: ", Xsvd$rank, Xh.svd$rank, k, kmax,"\n")
    
    ## perform extra reweighting step
    if(k != Xsvd$rank)
    {
      ## VT::27.08.2010 - bug report from Stephen Milborrow: if n is small relative to p
      ## k can be < Xsvd$rank but larger than Xh.svd$rank - the number of observations in
      ## Xh.svd is roughly half of these in Xsvd
      k <- min(Xh.svd$rank, k)
      
      #XRc <- X - repmat(Xh.svd$center, n, 1)
      XRc <- X - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)
      
      Xtilde <- XRc %*% Xh.svd$loadings[,1:k] %*% t(Xh.svd$loadings[,1:k])
      Rdiff <- XRc - Xtilde
      odh <- apply(Rdiff, 1, vecnorm)
      
      #Cutoff for the ODs
      tmp <- coOD(od=odh,h=h,skew=skew)
      odh <- tmp$od
      cutoffodh <- tmp$co.od
      
      indexset <- (odh <= cutoffodh)
      #Xh.svd <- classSVD(X[indexset,])
      Xh.svd <- classPC(X[indexset,])
      k <- min(Xh.svd$rank, k)
    } else {
      indexset <- H0
    }
    
    ## Project the data points on the subspace spanned by the first k0 eigenvectors
    center <- center + Xh.svd$center %*% t(rot)
    rot <- rot %*% Xh.svd$loadings
    #X2 <- (X - repmat(Xh.svd$center, n, 1)) %*% Xh.svd$loadings
    X2 <- (X - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)) %*% Xh.svd$loadings
    
    X2 <- as.matrix(X2[ ,1:k])
    rot <- as.matrix(rot[ ,1:k])
    
    if (skew) {
      
      #Adjusted outlyingness in k-dimensional subspace
      outl2 <- adjOutlyingness(X2, options = list(type="Rotation", ndir=ndir))$outlyingnessX
      H2 <- order(outl2)  #Order indices by outlyingness    
      Xh2 <- X2[H2[1:h], ] # the h data points with smallest outlyingness
      ee <- eigen(cov(Xh2))
      P6 <- ee$vectors
      X2center <- colMeans(Xh2)
      #MCD is replaced by mean and (sample) covariance of h points with smallest adjusted outlyingness
      #in the k-dimensional subspace.
      center <- as.vector(center + X2center %*% t(rot))
      eigenvalues <- ee$values
      loadings <- rot %*% P6
      scores <- sweep(X2,2,X2center,"-") %*% P6
      
      
    }else{
      
      
      ## Perform now MCD on X2 in order to obtain a robust scatter matrix: find h data points whose
      ##  covariance matrix has minimal determinant
      ##
      ## First apply C-step with the h points selected in the first step, i.e. those that
      ## determine the covariance matrix after the C-steps have converged.
      ##
      mah <- mahalanobis(X2, center=rep(0, ncol(X2)), cov=diag(Xh.svd$eigenvalues[1:k], nrow=k))
      oldobj <- prod(Xh.svd$eigenvalues[1:k])

      niter <- 100
      for(j in 1:niter) {
        if(trace)
          cat("\nIter=",j, " h=", h, " k=", k, " obj=", oldobj, "\n")

        Xh <- X2[order(mah)[1:h], ]
        Xh <- as.matrix(Xh)

        #Xh.svd <- classSVD(as.matrix(Xh))
        # Problems with classPC from robustbase when 1 column,
        # fix manually
        if(ncol(Xh)==1) {
          center <- mean(Xh)
          nn <- nrow(Xh)
          x <- scale(Xh, center=TRUE, scale=FALSE)
          svd <- svd(x/sqrt(nn-1), nu=0)
          loadings <- as.numeric(svd$v)
          if(loadings<0) loadings <- -loadings
          eigenvalues <- (svd$d)^2
          Xh.svd <- list(loadings=loadings, eigenvalues=eigenvalues, rank=1,
                         center=as.matrix(attr(x, "scaled:center")), scale=FALSE)
        } else {
          Xh.svd <- classPC(as.matrix(Xh))
        }

        obj <- prod(Xh.svd$eigenvalues)
        #X2 <- (X2 - repmat(Xh.svd$center, n, 1)) %*% Xh.svd$loadings
        X2 <- (X2 - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)) %*% Xh.svd$loadings
        
        
        center <- center + Xh.svd$center %*% t(rot)
        rot <- rot %*% Xh.svd$loadings
        #mah <- mahalanobis(X2, center=zeros(1, ncol(X2)), cov=diag(Xh.svd$eigenvalues, nrow=length(Xh.svd$eigenvalues)))
        mah <- mahalanobis(X2, center=matrix(0,1, ncol(X2)), cov=diag(Xh.svd$eigenvalues, nrow=length(Xh.svd$eigenvalues)))
        
        if(Xh.svd$rank == k & abs(oldobj - obj) < 1.E-12)
          break
        
        oldobj <- obj
        if(Xh.svd$rank < k) {
          j <- 1
          k <- Xh.svd$rank
        }
      }
      
      ## Perform now MCD on X2
      X2mcd <- CovMcd(X2, nsamp=250, alpha=alpha)
      if(trace)
        cat("\nMCD crit=",X2mcd@crit," and C-Step obj function=",obj," Abs difference=", abs(X2mcd@crit-obj), "\n")
      
      ## VT::14.12.2009 - if there is even a slight difference between mcd$crit and obj
      ## and it is on the negative side, the following reweighting step will be triggered,
      ## which could lead to unwanted difference in the results. Therefore compare with
      ## a tolerance 1E-16.
      eps <- 1e-16
      if(X2mcd@crit < obj + eps)
      {
        X2cov <- getCov(X2mcd)
        X2center <- getCenter(X2mcd)
        if(trace)
          cat("\nFinal step - PC of MCD cov used.\n")
      }else
      {
        consistencyfactor <- median(mah)/qchisq(0.5,k)
        mah <- mah/consistencyfactor
        weights <- ifelse(mah <= qchisq(0.975, k), TRUE, FALSE)
        
        ##          VT::27.08.2010 - not necessary, cov.wt is doing it ptoperly
        ##            wcov <- .wcov(X2, weights)
        wcov <- cov.wt(x=X2, wt=weights, method="ML")
        X2center <- wcov$center
        X2cov <- wcov$cov
        if(trace)
          cat("\nFinal step - PC of a reweighted cov used.\n")
      }
      
      ee <- eigen(X2cov)
      P6 <- ee$vectors
      
      center <- as.vector(center + X2center %*% t(rot))
      eigenvalues <- ee$values
      loadings <- rot %*% P6
      #scores <- (X2 - repmat(X2center, n, 1)) %*% P6
      scores <- (X2 - matrix(rep(X2center, times=n), nrow=n, byrow=TRUE)) %*% P6
      
    }
  }

  #Change names of columns and/or rows
  if(is.list(dimnames(data)) && !is.null(dimnames(data)[[1]]))
  {
    dimnames(scores)[[1]] <- dimnames(data)[[1]]
  } else {
    dimnames(scores)[[1]] <- 1:n
  }
  dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))
  dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))
  
  #Compute distances and cutoffs
  A <- distPCA(X=data, Tn=scores, P=loadings, l=eigenvalues, mu=center, h=h, skew=skew)
  if (skew) {
    #SDs are now adjusted outylingnesses
    sd <- outl2
    #Cutoff is defined similarly as cutoff for ODs
    cutoff.sd <- coOD(sd,skew=TRUE)$co.od
  } else {
    sd <- A$sd
    cutoff.sd <- A$cutoff.sd
  }
  od <- A$od
  cutoff.od <- A$cutoff.od
  flag.sd <- (sd<=cutoff.sd) #FALSE if outlying in PCA subspace
  flag.od <- (od<=cutoff.od) #FALSE if outlying for OD
  flag.all <- (flag.od*flag.sd)==1 #FALSE if outlier (all 3 types)
  
  #Results
  res <- list(loadings=loadings, eigenvalues=eigenvalues, scores=scores,  
              center=center, k=k, H0=H0, H1=indexset, alpha=alpha, h=h,
              sd=sd, od=od, cutoff.sd=cutoff.sd, cutoff.od=cutoff.od,
              flag.sd=flag.sd, flag.od=flag.od, flag.all=flag.all)
  return(res)
}
