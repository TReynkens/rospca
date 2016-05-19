\name{selectLambda}
\alias{selectLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of sparsity parameter using IC
}
\description{
Selection of the sparsity parameter for ROSPCA and SCoTLASS using BIC of Hubert et al. (2016), and for SRPCA using BIC of Croux et al. (2013).
}
\usage{
selectLambda(X, k, kmax = 10, method = "ROSPCA", lmin = 0, lmax = 2, lstep = 0.02,
             alpha = 0.75, stand = TRUE, skew = FALSE, multicore = FALSE, 
             mc.cores = NULL, P = NULL, ndir = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{An \eqn{n} by \eqn{p} matrix or data matrix with observations in the rows and variables in the columns.}
  \item{k}{
Number of Principal Components (PCs).
}
  \item{kmax}{
Maximal number of PCs to be computed, only used when \code{method = "ROSPCA"} or \code{method = "ROSPCAg"}. Default is 10.}
  \item{method}{
PCA method to use: ROSPCA (\code{"ROSPCA"} or \code{"ROSPCAg"}), SCoTLASS (\code{"SCoTLASS"} or \code{"SPCAg"}) or SRPCA (\code{"SRPCA"}). Default is \code{"ROSPCA"}.
}
  \item{lmin}{
Minimal value of \eqn{\lambda} to look at, default is 0.
}
  \item{lmax}{
Maximal value of \eqn{\lambda} to look at, default is 2.
}
  \item{lstep}{
Difference between two consecutive values of \eqn{\lambda}, i.e. the step size, default is 0.02.
}
  \item{alpha}{
Robustness parameter for ROSPCA, default is 0.75.
}
  \item{stand}{
Logical indicating if the data should be standardised, default is \code{TRUE}.
}
  \item{skew}{
Logical indicating if the skewed version of ROSPCA should be applied, default is \code{FALSE}.
}
  \item{multicore}{
Logical indicating if multiple cores can be used, default is \code{TRUE}. Note that this is not possible for the Windows platform, so \code{multicore} is always \code{FALSE} there.
}
  \item{mc.cores}{
Number of cores to use if \code{multicore=TRUE}, default is \code{NULL} which corresponds to the number of cores minus 1.
}
  \item{P}{
True loadings matrix, a numeric matrix of size \eqn{p} by \eqn{k}. The default is \code{NULL} which means that no true loadings matrix is specified. 
}
\item{ndir}{Number of directions used when computing the outlyingness (or the adjusted outlyingness when \code{skew=TRUE}) in \code{rospca}, see \code{\link[mrfDepthLight]{outlyingness}} and \code{\link[mrfDepthLight]{adjOutlyingness}} for more details.
}

}
\details{
We select an optimal value of \eqn{\lambda} for a certain method on a certain dataset by looking at an equidistant grid of \eqn{\lambda} values. For each value of \eqn{\lambda}, we apply the method on the dataset using this sparsity parameter, and compute an Information Criterion (IC). The optimal value of \eqn{\lambda} is then the one corresponding to the minimal IC. The ICs we consider are the BIC of for Hubert et al. (2016) for ROSPCA and SCoTLASS, and the BIC of Croux et al. (2013) for SRPCA.
}

\value{
A list with components:\cr
\item{opt.lambda}{Value of \eqn{\lambda} corresponding to minimal IC.}
\item{min.IC}{Minimal value of IC.}
\item{Lambda}{Numeric vector containg the used values of \eqn{\lambda}.}
\item{IC}{Numeric cector containing the IC values corresponding to all values of \eqn{\lambda} in \code{Lambda}.}
\item{loadings}{Loadings obtained using method with sparsity parameter \code{lambda_opt}, a numeric matrix of size \eqn{p} by \eqn{k}.}
\item{fit}{Fit obtained using method with sparsity parameter \code{lambda_opt}. This is a list containing the loadings (\code{loadings}), the eigenvalues (\code{eigenvalues}), the standardised data matrix used as input (\code{Xst}), the scores matrix (\code{scores}), the orthogonal distances (\code{od}) and the score distances (\code{sd}).}
\item{type}{Type of IC used: \code{BICod} (BIC of Hubert et al. (2016)) or \code{BIC} (BIC of Croux et al. (2013)).}
\item{measure}{A numeric vector containing the standardised angles between the true and the estimated loadings matrix for each value of \eqn{\lambda} if a loadings matrix is given. When no loadings matrix is given as input (\code{P=NULL}), \code{measure} is equal to \code{NULL}.}
}
\references{
Hubert, M., Reynkens, T., Schmitt, E. and Verdonck, T. (2016). ``Sparse PCA for High-Dimensional Data with Outliers,'' \emph{Technometrics}, \emph{accepted}.

Croux, C., Filzmoser, P., and Fritz, H. (2013), ``Robust Sparse Principal Component Analysis,'' \emph{Technometrics}, 55, 202--214.

}
\author{
Tom Reynkens
}



\seealso{
\link{selectPlot}, \link[parallel]{mclapply}, \link{angle}
}
\examples{
X <- dataGen(m=1,n=100,p=10,eps=0.2,bLength=4)$data[[1]]

sl <- selectLambda(X,k=2,method="ROSPCA",lstep=0.1)
selectPlot(sl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }

