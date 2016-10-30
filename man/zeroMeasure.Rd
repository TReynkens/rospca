\name{zeroMeasure}
\alias{zeroMeasure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zero measure
}
\description{
 Compute the average zero measures and total zero measure for a list of matrices.
}
\usage{
zeroMeasure(Plist, P, prec = 10^(-5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Plist}{
List of estimated loadings matrices or a single estimated loadings matrix.
All these matrices should be numeric matrices of size \eqn{p} by \eqn{k}.
}
  \item{P}{
True loadings matrix, a numeric matrix of size \eqn{p} by \eqn{k}.
}
  \item{prec}{
Precision used when determining if an element is non-zero, default is \eqn{10^{-5}}.
We say that all elements with an absolute value smaller than \code{prec} are ``equal to zero''.
}
}
\details{
The \emph{zero measure} is a way to compare how correctly a PCA method estimates the sparse loadings matrix \code{P}. For each element of an estimated loadings matrix, it is equal to one if the estimated and true value are both zero or both non-zero, and 0 otherwise. We then take the average zero measure over all elements of an estimated loadings matrix and over all estimated loadings matrices which we call the \emph{total zero measure}.
}
\value{
A list with components:\cr
  \item{measure}{Numeric matrix of size \eqn{p} by \eqn{k} containg the average zero measure over all \code{length(Plist)} simulations for each element of \code{P}.}
  \item{index}{Numeric vector containing the indices of all data sets where the estimate was wrong (at least one of the zero measures for the elements of an estimated loadings matrix is equal to 0).}
  \item{total}{Total zero measure, i.e. the average zero measure over all elements of an estimated loadings matrix and over all estimated loadings matrices.}
}
\references{
Hubert, M., Reynkens, T., Schmitt, E. and Verdonck, T. (2016). ``Sparse PCA for High-Dimensional Data with Outliers,'' \emph{Technometrics}, 58, 424--434.
}
\author{
Tom Reynkens
}

\examples{
P <- cbind(c(1,1),c(0,1))
Plist <- list(matrix(1,2,2),P)

zeroMeasure(Plist,P)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ algebra }% __ONLY ONE__ keyword per line
