
###############################################
# General function for simulations             #
###############################################

# Auxiliary output functions
# res is a result of class sPCAgrid
# X is the data frame
# k is the number of principal components
aux_output <- function(res, X, k) {
  
  eigenvalues <- (res$sdev^2)[1:k]
  distances <- distPCA(X=X, Tn=res$scores, P=res$loadings, l=eigenvalues)
  return(list(loadings=res$loadings, eigenvalues=eigenvalues, Xst=X, scores=res$scores, od=distances$od, sd=distances$sd))

}

# Selection of method, the function f_pca will always at least return the loadings matrix, 
# the eigenvalues and Xst.
# lambda_def is the default value of lambda for the Grid based sparse methods
# para is the para option for the sparse methods that use spca
# mcd is an option for ROBPCA
select_method <- function (method = "CPCA", k = 2, alpha = 0.75, lambda_def = 0, mcd = FALSE, stand = TRUE) {
  lambda <- NULL
  if (method=="SCoTLASS") { method <- "SPCAg"}
  if (method=="ROSPCAg") { method <- "ROSPCA"}
  
  f_pca <- switch(method,"CPCA" = {
    function (X) { 
      if (stand) {Xst <- matStand(X, f_c=mean, f_s=sd)$data
      } else {Xst <- X}
      A <- princomp(Xst)
      return(list(loadings=A$loadings[,1:k], eigenvalues=(A$sdev[1:k])^2, Xst=Xst, scores=A$scores))}
   
  }, "RPCA"={
    function (X) { 
      if (stand) {Xst=matStand(X, f_c=median, f_s=qn)$data
      } else {Xst <- X}
      A <- sPCAgrid(Xst, k=k, method="qn", lambda=0, scale=NULL, center=rep(0, ncol(Xst)), glo.scatter=1, maxiter=75)
      return(aux_output(res=A, X=Xst, k=k))}   
      
  }, "ROBPCA"={
    function (X) { 
      if (stand) {Xst=matStand(X, f_c=median, f_s=qn)$data
      } else {Xst <- X}
      return(robpca(Xst, k=k, alpha=alpha, mcd=mcd))}
    
  }, "SPCAg"={
    function (X, lambda=lambda_def) { 
      if (stand) {Xst=matStand(X, f_c=mean, f_s=sd)$data
      } else {Xst <- X}
      A <- sPCAgrid(Xst, k=k, method="sd", lambda=lambda, scale=NULL, center=rep(0, ncol(Xst)), 
                    glo.scatter=1, maxiter=75)
      return(aux_output(res=A, X=Xst, k=k))}
    
  }, "SRPCA"={
    function (X, lambda=lambda_def) { 
      if (stand) {Xst=matStand(X, f_c=median, f_s=qn)$data
      } else {Xst <- X}
      A <- sPCAgrid(Xst, k=k, method="qn", lambda=lambda, scale=NULL, center=rep(0, ncol(Xst)), 
                    glo.scatter=1, maxiter=75)
      return(aux_output(res=A, X=Xst, k=k))}
    
  }, "ROSPCA"={
    # We can directly return the output from rospca because it is in a list format.
    # The previous methods have an S4 object as output or are not suitable for direct return.
    function (X, lambda=lambda_def) { 
      return(rospca(X, k=k, lambda=lambda, stand=stand, alpha=alpha, grid=TRUE))}
    

  }, stop("Invalid method name, use \"CPCA\", \"RPCA\", \"ROBPCA\", \"SPCAg\" (or \"SCoTLASS\"), \"SRPCA\",
        \"ROSPCA\" (or \"ROSPCAg\")."))
    return(f_pca)
}

# Generate correlation matrix where the blocks have off-diagonal elements equal to the elements of the vector a. 
# p is the number of dimensions
# bLength is the size of a block
# The last block is equal to the identity matrix and can have a different size (size=p-bLength*length(a))
genMat  <-  function (a = c(0.9,0.5,0), p = 10, bLength = 4) { 
  
  k <- length(a)-1
  if (p<=bLength*k) {stop("Dimension needs to be at least bLength*length(a)")}
  
  if (any(abs(a)>1)) {
    stop("Off-diagonal elements need to be smaller than 1 in absolute value!")}
    
  X <- diag(p) # diagonal zero matrix of size p
  
  for (t in 0:(k-1)) {
    for (i in (t*bLength+1):((t+1)*bLength)) {
      for( j in (t*bLength+1):((t+1)*bLength)) {
        if (j != i) {X[i,j]=a[t+1]}
      }
    }
  }
  # Off-diagonal elements of last group
  if (a[k+1]!=0) {
    for (i in (k*bLength+1):p) {
      for (j in (k*bLength+1):p) {
        if (j != i) {X[i,j]=a[k+1]}
      }
    }
  }
  return(X)
}

# Generate data 
# a and bLength are options for the genMat function
dataGen <- function (m = 100, n = 100, p = 10, a = c(0.9,0.5,0), bLength = 4, SD = c(10,5,2), eps = 0, seed = TRUE) {
  
  data <- list() # List of m data matrices
  ind <- list() # List of vectors containing the indices of the outliers
  
  if (p<8) {stop("Dimension needs to be at least 8")}
  
  # Correlation matrix
  R <- genMat(p=p, a=a, bLength=bLength)
  
  k <- length(a)-1
  if (k<2) {stop("a should have length at least 3.")}
  if (length(SD)!=length(a)) {stop("SD have the same length as a")}
  
  sdx <- rep(0, p)
  for (i in 1:k) {
    # The ith block (with length bLength) gets st. dev SD[i]
    sdx[(i*bLength)-(0:(bLength-1))] <- rep(SD[i], bLength)
  }
  # Last block (with length p-(k*bLength)) gets st. dev SD[k+1]
  sdx[(k*bLength+1):p] <- rep(SD[k+1], p-k*bLength)       
  SDx <- diag(sdx)
  # Covariance matrix
  Sigma <- SDx%*%R%*%SDx

  # Convert epsilon when given in %
  if (eps>=1) {eps <- eps/100}
  
  # Mean of outliers
  mu_out <- 25*c(c(0,-4,4,2,0,4,-4,2), rep(c(3,-3), length.out=p-8))
  # Number of outliers
  n_out <- floor(eps*n)
  # Variance of outliers
  sigma2 <- 20
  
  for (i in 1:m) {
    if (seed) {set.seed(23071990+i)}
    
    X <- rmvnorm(n=n, mean=rep(0, p), sigma=Sigma)
    # Add noise (from a p-variate standard normal distribution)
    X <- X+rmvnorm(n=n, mean=rep(0, p), sigma=diag(p))
     
    if (eps>0) {
      # Outliers (select indices of rows to replace)
      index <- sample(1:n, n_out, replace=FALSE)
      X_out <- rmvnorm(n=n_out, mu_out, sigma=sigma2*diag(p))
    
      X[index,] <- X_out
     } else {index=0}
    
    ind[[i]] <- index
    data[[i]] <- X
  }
  return(list(data=data, ind=ind, R=R, Sigma=Sigma))
}



###############################################
# Angles                                       #
###############################################

# Last principal angle between A and B
# We use Bjorck & Golub, Numerical methods for computing angles between linear subspaces,
#  Math. Comp. 27 (1973), pp. 579-594.
# More stable to compute than maxsub 
#  (eigenvalues of matrix in maxsub can be outside of [0,1]).
# Based on matlab function subspace.
angle <- function (A, B) { 
  
  A <- pracma::orth(A)
  B <- pracma::orth(B)
  # Check rank and swap
  if (ncol(A) < ncol(B)) {
    tmp <- A
    A <- B
    B <- tmp
  }
  # Compute the projection according to Equation 13 in Bjorck and Golub (1973).
  B <- B - A%*%(t(A)%*%B)
  # Make sure its magnitude is less than 1.
  theta <- asin(min(1, max(svd(B)$d)))
  return(theta/(pi/2))
}


# The maximal angle is defined as the arccos of the smallest eigenvalue of the matrix
# Z as in Krzanowski (1979).
maxsub <- function (P, PP) {
  
  P <- P[,1:ncol(PP)]
  Z <- t(PP)%*%P%*%t(P)%*%PP
  
  # Eigenvalues can be different from 1 if equal due to numerical problems
  if (sum(abs(Z-diag(ncol(PP))))<10^{-14}) {  
    return(0)
  } else {
    # Take absolute value because some eigenvalues can be very close to 0
    # but negative due to numerical problems. In theory all eigenvalues are between 0 and 1
    return(  acos( sqrt(min(abs(eigen(Z)$values))) )/(pi/2)  )
  }
}

# The minimal angle is defined as the arccos of the largest eigenvalue of the matrix
# Z as in Krzanowski (1979).
minsub <- function (P, PP) {
  
  P <- P[,1:ncol(PP)]
  Z <- t(PP)%*%P%*%t(P)%*%PP
  
  # Eigenvalues can be different from 1 if equal due to numerical problems
  if (sum(abs(Z-diag(ncol(PP))))<10^{-14}) { 
    return(0)
  } else {
    # Take absolute value because some eigenvalues can be very close to 0
    # but negative due to numerical problems. In theory all eigenvalues are between 0 and 1
    return(  acos( sqrt(max(abs(eigen(Z)$values))) )/(pi/2)  )
  }
}


# The total similarity is defined as the sum of the eigenvalues of the matrix
# Z as in Krzanowski (1979).
totSimil <- function (P, PP) {
  
  P <- P[,1:ncol(PP)]
  Z <- t(PP)%*%P%*%t(P)%*%PP
  
  # Eigenvalues can be different from 1 if equal due to numerical problems
  if (sum(abs(Z-diag(ncol(PP))))<10^{-14}) {  
    return(0)
  } else {
    return(sum(eigen(Z)$values))
  }
}


###############################################
# Zero measure                                 #
###############################################

# Compare if two elements are both non-zero or both zero.
# prec is the used precision ("x=0" iff abs(x)<=prec)
# if a and b are vectors or matrices, the comparison is done elementwise
zeroCompare <- function (a, b, prec = 10^(-5)) {
  
  x <- abs(a)
  y <- abs(b)
  # Note that & and | are elementwise operators!
  return((x<=prec & y<=prec) | (x>prec & y>prec))
}


# Zero measure for loadings. For each loading we take the average over all indicator functions
# of the event "estimated loading is 0 <=> true loading is 0".
# prec is the used precision ("x=0" iff abs(x)<=prec)
# P is the true loadings matrix and Plist is a list of the estimated loadings matrices
zeroMeasure <- function (Plist, P, prec = 10^(-5)) {
  if (!is.list(Plist)) {
    tmp <- Plist
    Plist <- list()
    Plist[[1]] <- tmp
  }
  if (!all(dim(P)==dim(Plist[[1]]))) {
    stop('P and Plist should have the same dimensions')
  }
  
  m <- length(Plist)
  k <- ncol(Plist[[1]])
  p <- nrow(Plist[[1]])
  
  wrong <- rep(0, m)
  
  # Run over simulations
  for (t in 1:m) {  
    # Overwrite Plist by zero measure of all elements of the matrix
    Plist[[t]] <- zeroCompare(Plist[[t]], P, prec)
    # Wrong if at least one of the measures is different from 1
    wrong[t] <- any(Plist[[t]]==0)
  }
  
  # Element by element mean of a list of matrices
  s <- apply(simplify2array(Plist), 1:2, mean)
  
  st <- sapply(1:k, function (i) { paste0("PC", i)})
  dimnames(s) <- list(NULL, st) # no row names
  
  # indices of all data sets where the estimate was wrong (0 in the zero measure)
  index <- which(wrong==1)
  
  # Total zero measure
  total <- mean(s)
  
  return(list(measure=s, index=index, total=total))
}

