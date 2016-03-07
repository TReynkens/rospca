\name{rospca}
\alias{rospca}

\title{
RObust Sparse PCA algorithm
}

\description{
Sparse robust PCA algorithm based on the ROBPCA algorithm of Hubert et al. (2005).
}

\usage{
rospca(X, k, kmax = 10, alpha = 0.75, h = NULL, ndir = "all", grid = TRUE, 
       lambda = 10^(-6), sparse = "varnum", para, stand = TRUE, skew = FALSE)
}

\arguments{

\item{X}{An \eqn{n} by \eqn{p} matrix or data matrix with observations in the rows and variables in the columns.}
\item{k}{Number of principal components that will be used.}
\item{kmax}{Maximal number of principal components that will be computed, default is 10.} 
\item{alpha}{Robustness parameter, default is 0.75.}
\item{h}{The number of outliers the algorithm should resist is given by \eqn{(n-h+1)}. Any value for \code{h} between \eqn{n/2} and \eqn{n} may be specified. Default is \code{NULL} which uses \code{h=ceiling(alpha*n)+1}. Do not specify \code{alpha} and \code{h} at the same time. }
\item{ndir}{Number of directions used when computing the outlyingness (or the adjusted outlyingness when \code{skew=TRUE}), see \code{\link[mrfDepthLight]{outlyingness}} and \code{\link[mrfDepthLight]{adjOutlyingness}} for more details.}
\item{grid}{Logical indicating if the grid version of sparse PCA should be used (\code{SPcaGrid} with \code{method="sd"} from \pkg{rrcovHD}). Otherwise, the version of Zou et al. (2006) is used (\code{spca} from \pkg{elasticnet}). Default is \code{TRUE}.}
\item{lambda}{Sparsity parameter of \code{SPcaGrid} (when \code{grid=TRUE}) or ridge parameter of \code{spca} (when
      \code{grid=FALSE}), default is \eqn{10^{-6}}.}
\item{sparse}{Parameter for \code{spca} (only used when \code{grid=FALSE}), see \code{\link[elasticnet]{spca}} for more details.}
\item{para}{Parameter for \code{spca} (only used when \code{grid=FALSE}), see \code{\link[elasticnet]{spca}} for more details.}
\item{stand}{If \code{TRUE}, the data are standardised robustly in the beginning and classically before applying sparse 
            PCA. If \code{FALSE}, the data are only mean-centred before applying sparse PCA. Default is \code{TRUE}.}  
\item{skew}{Logical indicating if the version for skewed data should be applied, default is \code{FALSE}.}
}

\details{
See Hubert et al. (2016) for more details. In contrast to the paper, we allow for SPCA (Zou et al., 2006) to be used as the sparse PCA method inside ROSPCA (\code{grid=FALSE}). Moreover, we also include a skewed version of ROSPCA (\code{skew=TRUE}) similar to the skewed version of ROBPCA (Hubert et al., 2009) which is not detailed in Hubert et al. (2016).
}

\value{
A list with components:\cr
\item{loadings}{Loadings matrix containing the sparse robust loadings (eigenvectors), a numeric matrix of size \eqn{p} by \eqn{k}.}
  \item{eigenvalues}{Numeric vector of length \eqn{k} containing the robust eigenvalues.}
  \item{scores}{Scores matrix (computed as \eqn{(X-center)*loadings)}, a numeric matrix of size \eqn{n} by \eqn{k}.}
  \item{center}{Numeric vector of length \eqn{k} containing the centre of the data.}
  \item{D}{Matrix used to standardise the data before applying sparse PCA (identity matrix if \code{stand=FALSE}), a numeric matrix of size \eqn{p} by \eqn{p}.}
  \item{k}{Number of (chosen) principal components.}
  \item{H0}{Logical vector of size \eqn{n} indicating if an observation is in the initial h-subset.}
   \item{H1}{Logical vector of size \eqn{n} indicating if an observation is kept in the non-sparse reweighting step (in robust part).}
  \item{P1}{Loadings matrix before applying sparse reweighting step, a numeric matrix of size \eqn{p} by \eqn{k}.}
  \item{index}{Numeric vector containing the indices of the variables that are used in the sparse reweighting step.}
  \item{H2}{Logical vector of size \eqn{n} indicating if an observation is kept in the sparse reweighting step.}
  \item{P2}{Loadings matrix before estimating eigenvalues, a numeric matrix of size \eqn{p} by \eqn{k}.}
  \item{H3}{Logical vector of size \eqn{n} indicating if an observation is kept in the final SD reweighting step.}
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
Hubert, M., Reynkens, T., Schmitt, E. and Verdonck, T. (2016). ``Sparse PCA for High-Dimensional Data with Outliers,'' \emph{Technometrics}, \emph{accepted}.

Hubert, M., Rousseeuw, P. J., and Vanden Branden, K. (2005), ``ROBPCA: A New Approach to Robust Principal Component Analysis,'' \emph{Technometrics}, 47, 64--79.

Hubert, M., Rousseeuw, P. J., and Verdonck, T. (2009), ``Robust PCA for Skewed Data and Its Outlier Map," \emph{Computational Statistics & Data Analysis}, 53, 2264--2274.

Croux, C., Filzmoser, P., and Fritz, H. (2013), ``Robust Sparse Principal Component Analysis,'' \emph{Technometrics}, 55, 202--214.

Zou, H., Hastie, T., and Tibshirani, R. (2006), ``Sparse Principal Component Analysis,'' \emph{Journal of Computational and Graphical Statistics}, 15, 265--286.
}

\author{
Tom Reynkens, based on R code from Valentin Todorov for \code{PcaHubert} in \pkg{rrcov} (released under GPL-3).
}

%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[rrcov]{PcaHubert}}, \code{\link{robpca}}, \code{\link[mrfDepthLight]{outlyingness}}, \code{\link[mrfDepthLight]{adjOutlyingness}}, \code{\link[rrcovHD]{SPcaGrid}}, \code{\link[elasticnet]{spca}}
}

\examples{
X <- dataGen(m=1,n=100,p=10,eps=0.2,bLength=4)$data[[1]]

resRS <- rospca(X,k=2,lambda=0.4,stand=TRUE)
diagPlot(resRS)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }% __ONLY ONE__ keyword per line
\keyword{ multivariate }% __ONLY ONE__ keyword per line