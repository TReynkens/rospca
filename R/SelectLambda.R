
#########################################################
# Selection of lambda using AIC or BIC 
#########################################################

# Number of non-zero elements of a matrix
nnz <- function (X, tol = 10^(-5)) {
  
  return(sum(abs(X)>tol))
}


# BIC criterion based on h1 smallest ODs
# Idea from Maronna (2005): "Principal components and orthogonal regression based on robust scales"
# and Allen & Maletic-Savatic (2011): "Sparse non-negative generalized PCA 
#                                      with applications to metabolomics"
#
# res is the output of a sparse (robust) PCA method
# res0 is the output of the unconstrained (robust) PCA method (so with lambda=0)
# h1 is the number of ODs to look at

BICod <- function (res, h1) {
  
  odh <- sort(res$od)[1:h1]
  df <- nnz(res$loadings)
  a <- nrow(res$loadings)*h1
  return(log(sum(odh^2)/a)+log(a)*df/a)
  
}


# Total variance of matrix X using the square of the scale measure fs
matvar <- function (X, fs = sd) {
  
  return(sum(apply(X, 2, fs)^2))
}

# BIC criterion 
# Xcen is the centred data matrix and res is output from a PCA method
# V0 is the residual variance of the method with lambda=0
# fs is the used scale measure
Bic <- function (Xcen, res, V0, fs = sd, aic = FALSE) {
  
  P <- res$loadings
  V <- matvar(Xcen-Xcen%*%P%*%t(P), fs=fs)
  df <- nnz(P)
  n <- nrow(Xcen)
  if (aic) {
    return(V/V0+2/n*df) 
  } else {
    return(V/V0+log(n)/n*df) 
  }
  
}

# Selection of lambda using information criterion (IC)
selectLambda <- function (X, k, kmax = 10, method = "ROSPCA", lmin = 0, lmax = 2, lstep = 0.02,
                         alpha = 0.75, stand = TRUE, skew = FALSE, multicore = FALSE, mc.cores = NULL, P = NULL,
                         ndir = "all") {
  X <- as.matrix(X)

  Lambda <- seq(lmin, lmax, by=lstep)
  l <- length(Lambda)
  
  # Faster version because part of calculations are independent of lambda
  if (method %in% c("ROSPCA","ROSPCAg") ) {
    temp <- rospca_part1(X, k=k, alpha=alpha, stand=stand, skew=skew, kmax=kmax, ndir=ndir)
    h1 <- sum(temp$H1)
    f_pca <- function (Y, lambda = 0) {
      return(rospca_part2(Y, H0=temp$H0, H1=temp$H1, k=temp$k, alpha=temp$alpha, h=temp$h, lambda=lambda,
                          stand=stand, skew=skew))}
    
    f_IC <- function (x) {BICod(x, h1=h1)}; type="BICod"
    
  } else if (method %in% c("SPCAg", "SCoTLASS")) {
    h1 <- nrow(X)
    f_pca <- select_method(method=method, k=k, stand=stand)
    f_IC <- function (x) {BICod(x, h1=h1)}; type="BICod"
    
  } else if (method=="SRPCA") {
    h1 <- floor(alpha*nrow(X))
    f_pca <- select_method(method=method, k=k, stand=stand)
    res0 <- f_pca(X, lambda = 0)
    P <- res0$loadings
    Xcen <- matStand(X, f_c=median, scale=FALSE)$data
    V0 <- matvar(Xcen-Xcen%*%P%*%t(P), fs=qn)
    
    f_IC <- function (x) {Bic(x, fs=qn, V0=V0, Xcen=Xcen)}; type="BIC"
    
  } else {
    stop("Invalid method, use \"ROSPCAg\", \"SPCAg\", \"SCoTLASS\" or \"SRPCA\".")
  }
  
  if (is.null(P)) {
    f <- function (l, Y) {
      return(f_IC(f_pca(Y, lambda=l)))
    }
  } else {
    # Check if correct input
    if (ncol(X)!=nrow(P) | ncol(P)!=k) {
      stop("P has the wrong size.")
    }
    # Include angle
    f <- function (l, Y) {
      res <- f_pca(Y, lambda=l)
      return(c(f_IC(res), angle(P, res$loadings)))
    }
    
  }
  
  
  if (Sys.info()[['sysname']]=="Windows" | !multicore) {
    results <- sapply(Lambda, FUN=f, Y=X)
  } else {
    if (is.null(mc.cores)) mc.cores=detectCores()-1
    # do.call(cbind, list) transforms the results into matrix where every column contains the elements of a list element
    results <- do.call(cbind, mclapply(Lambda, FUN=f, Y=X, mc.cores=mc.cores))
  }
  
  
  if (is.null(P)) {
    ic <- results
  } else {
    ic <- results[1,]
  }
  
  ic <- as.numeric(ic)
  
  # optimal index is index corresponding to smallest IC
  opt_ind <- which.min(ic)
  lambda_opt <- Lambda[opt_ind]
  min.ic <- ic[opt_ind]
  
  res <- f_pca(X, lambda=lambda_opt)
  
  lst <- list(opt.lambda=lambda_opt, min.IC=min.ic, Lambda=Lambda, IC=ic, loadings=res$loadings, fit=res, type=type)
  
  if (is.null(P)) {
    lst$measure <- NULL
  } else {
    lst$measure <- results[2,]
  }
  
  return(lst)
}


# sl is the output from selectLambda
# indicate indicates if the optimal lambda value according to the IC is indicated on the IC curve
selectPlot <- function (sl, indicate = TRUE, main = NULL) {
  
  type <- sl$type
  plot(sl$Lambda, as.numeric(sl$IC), type="b", xlab=bquote(lambda), ylab=type, main=main)
  if (indicate) {
    lambda <- sl$opt.lambda
    abline(v=lambda, lty=2, lwd=2)
    mtext(lambda, at=lambda, side=1)
  }
}

