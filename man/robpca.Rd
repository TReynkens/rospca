\name{robpca}
\alias{robpca}

\title{
ROBust PCA algorithm
}

\description{
ROBPCA algorithm of Hubert et al. (2005) including reweighting (Engelen et al., 2005) and possible extension to skewed data (Hubert et al., 2009).
}

\usage{
robpca (x, k = 0, kmax = 10, alpha = 0.75, h = NULL, mcd = FALSE, 
        ndir = "all", skew = FALSE, ...)
}

\arguments{

\item{x}{An \eqn{n} by \eqn{p} matrix or data matrix with observations in the rows and variables in the columns.}
\item{k}{Number of principal components that will be used. When \code{k=0} (default), the number of components is selected using the criterion in Hubert et al. (2005).}
\item{kmax}{Maximal number of principal components that will be computed, default is 10.} 
\item{alpha}{Robustness parameter, default is 0.75.}
\item{h}{The number of outliers the algorithm should resist is given by \eqn{n-h}. Any value for \code{h} between \eqn{n/2} and \eqn{n} may be specified. Default is \code{NULL} which uses \code{h=ceiling(alpha*n)+1}. Do not specify \code{alpha} and \code{h} at the same time. }
\item{mcd}{Logical indicating if the MCD adaptation of ROBPCA may be applied when the number of variables is sufficiently small (see details). If \code{mcd=FALSE} (default), the full ROBPCA algorithm is always applied.}
\item{ndir}{Number of directions used when computing the outlyingness (or the adjusted outlyingness when \code{skew=TRUE}), see \code{\link[mrfDepthLight]{outlyingness}} and \code{\link[mrfDepthLight]{adjOutlyingness}} for more details.}
\item{skew}{Logical indicating if the version for skewed data (Hubert et al., 2009) is applied, default is \code{FALSE}.}
\item{...}{Other arguments to pass to methods.}
}

\details{
When \code{mcd=TRUE} and \eqn{n<5*p}, we do not apply the full ROBPCA algorithm. The loadings and eigenvalues
are then computed as the eigenvectors and eigenvalues of the MCD estimator applied to data set after the SVD step.
An extensive description of the ROBPCA algorithm can be found in Hubert et al. (2005).

This funtion is based extensively on \code{PcaHubert} from \pkg{rrcov}. Notable differences are the use of outlyingness measures from \pkg{mrfDepthLight}, the implementation of the extensions for skewed data and a different output format.
}

\value{
A list with components:\cr
  \item{loadings}{Loadings matrix containing the sparse robust loadings (eigenvectors), a numeric matrix of size \eqn{p} by \eqn{k}.}
  \item{eigenvalues}{Numeric vector of length \eqn{k} containg the robust eigenvalues.}
  \item{scores}{Scores matrix (computed as \eqn{(X-center)*loadings)}, a numeric matrix of size \eqn{n} by \eqn{k}.}
  \item{center}{Numeric vector of length \eqn{k} containg the centre of the data.}
  \item{k}{Number of (chosen) principal components.}
  \item{H0}{Logical vector of size \eqn{n} indicating if an observation is in the initial h-subset.}
  \item{H1}{Logical vector of size \eqn{n} indicating if an observation is kept in the reweighting step.}
  \item{alpha}{The robustness parameter \eqn{\alpha} used throughout the algorithm.}  
  \item{h}{The \eqn{h}-parameter used throughout the algorithm.}
  \item{sd}{Numeric vector of size \eqn{n} containing the robust score distances within the robust PCA subspace.}
  \item{od}{Numeric vector of size \eqn{n} containing the orthogonal distances to the robust PCA subspace.}
  \item{cutoff.sd}{Cut-off value for the robust score distances.}
  \item{cutoff.od}{Cut-off value for the orthogonal distances.}
  \item{flag.sd}{Numeric vector of size \eqn{n} containing the SD-flags of the observations. The observations whose score distance is larger than \code{cutoff.sd} receive an SD-flag equal to zero. The other observations receive an SD-flag equal 1.}  
  \item{flag.od}{Numeric vector of size \eqn{n} containing the OD-flags of the observations. The observations whose orthogonal distance is larger than \code{cutoff.od} receive an OD-flag equal to zero. The other observations receive an OD-flag equal 1.}
  \item{flag.all}{Numeric vector of size \eqn{n} containing the flags of the observations. The observations whose score distance is larger than  \code{cutoff.sd} or whose orthogonal distance is 
    larger than  \code{cutoff.od} can be considered as outliers and receive a flag equal to zero. 
    The regular observations receive flag 1.}
}

\references{
Hubert, M., Rousseeuw, P. J., and Vanden Branden, K. (2005), ``ROBPCA: A New Approach to Robust Principal Component Analysis,'' \emph{Technometrics}, 47, 64--79.

Engelen, S., Hubert, M. and Vanden Branden, K. (2005), ``A Comparison of Three Procedures for Robust PCA in
High Dimensions", \emph{Austrian Journal of Statistics}, 34, 117--126.

Hubert, M., Rousseeuw, P. J., and Verdonck, T. (2009), ``Robust PCA for Skewed Data and Its Outlier Map," \emph{Computational Statistics & Data Analysis}, 53, 2264--2274.

}

\author{
Tom Reynkens, based on R code from Valentin Todorov for \code{PcaHubert} in \pkg{rrcov} (released under GPL-3).
}

%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[rrcov]{PcaHubert}}, \code{\link[mrfDepthLight]{outlyingness}}, \code{\link[mrfDepthLight]{adjOutlyingness}}
}


\examples{
X <- dataGen(m=1,n=100,p=10,eps=0.2,bLength=4)$data[[1]]

resR <- robpca(X,k=2)
diagPlot(resR)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }% __ONLY ONE__ keyword per line
\keyword{ multivariate }% __ONLY ONE__ keyword per line