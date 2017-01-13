\name{dataGen}
\alias{dataGen}

\title{
Generate sparse data with outliers
}

\description{
Generate sparse data with outliers using simulation scheme detailed in Hubert et al. (2016).
}

\usage{
dataGen(m = 100, n = 100, p = 10, a = c(0.9,0.5,0), bLength = 4, SD = c(10,5,2), 
        eps = 0, seed = TRUE)
}

\arguments{

\item{m}{Number of datasets to generate, default is 100.}
\item{n}{Number of observations, default is 100.}
\item{p}{Number of dimensions, default is 10.} 
\item{a}{Numeric vector containing the inner group correlations for each block. The number of useful blocks is thus given by \eqn{k=length(a)-1} which should be at least 2. By default, the correlations are equal to 0.9, 0.5 and 0, respectively.}
\item{bLength}{Length of the blocks of useful variables, default is 4.}
\item{SD}{Numeric vector containing the standard deviations of the blocks of variables, default is \code{c(10,4,2)}. Note that \code{SD} and \code{a} should have the same length.}
\item{eps}{Proportion of contamination, should be between 0 and 0.5. Default is 0 (no contamination).}
\item{seed}{Logical indicating if a seed is used when generating the datasets, default is \code{TRUE}.}
}

\details{Firstly, we generate a correlation matrix such that it has sparse eigenvectors.
We design the correlation matrix to have \eqn{length(a)=k+1} groups of variables with no correlation between variables from different groups. The first \eqn{k} groups consist of \code{bLength} variables each. The correlation between the different variables of the group is equal to \code{a[1]} for group 1, .... .  The (k+1)th group contains the remaining \eqn{p-k \times bLength} variables, which we specify to have correlation \code{a[k+1]}.  \cr
Secondly, the correlation matrix \code{R} is transformed into the covariance matrix \eqn{\Sigma= V^{0.5} \cdot R \cdot V^{0.5}}, where \eqn{V=diag(SD^2)}.\cr
Thirdly, the \code{n} observations are generated from a \eqn{p}-variate normal distribution with mean the \eqn{p}-variate zero-vector and covariance matrix \eqn{\Sigma}. Standard normally distributed noise terms are also added to each of the \code{p} variables to make the sparse structure of the data harder to detect.\cr
Finally, \eqn{(100 \times eps)\%} of the data points are randomly replaced by outliers. 
These outliers are generated from a \eqn{p}-variate normal distribution as in Croux et al. (2013). \cr
The \eqn{i}th eigenvector of \eqn{R}, for \eqn{i=1,...,k}, is given by a (sparse) vector with the \eqn{(bLength \times (i-1)+1)}th till the \eqn{(bLength \times i)}th elements equal to \eqn{1/\sqrt{bLength}} and all other elements equal to zero.\cr
See Hubert et al. (2016) for more details.
}

\value{
A list with components:\cr
  \item{data}{List of length \eqn{m} containing all data matrices.}
  \item{ind}{List of length \eqn{m} containing the numeric vectors with the indices of the contaminated observations.}
  \item{R}{Correlation matrix of the data, a numeric matrix of size \eqn{p} by \eqn{p}.}
  \item{Sigma}{Covariance matrix of the data (\eqn{\Sigma}), a numeric matrix of size \eqn{p} by \eqn{p}.}
}

\references{
Hubert, M., Reynkens, T., Schmitt, E. and Verdonck, T. (2016). ``Sparse PCA for High-Dimensional Data with Outliers,'' \emph{Technometrics}, 58, 424--434.

Croux, C., Filzmoser, P., and Fritz, H. (2013), ``Robust Sparse Principal Component Analysis,'' \emph{Technometrics}, 55, 202--214.
}

\author{
Tom Reynkens
}


\examples{
X <- dataGen(m=1, n=100, p=10, eps=0.2, bLength=4)$data[[1]]

resR <- robpca(X, k=2, skew=FALSE)
diagPlot(resR)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}