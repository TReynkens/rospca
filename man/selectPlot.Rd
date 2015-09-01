\name{selectPlot}
\alias{selectPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection plot
}
\description{
Plot Information Criterion (IC) versus values of the sparsity parameter \eqn{\lambda}.
}
\usage{
selectPlot(sl, indicate = TRUE, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sl}{
Output from \code{selectLambda} function.
}
  \item{indicate}{
Logical indicating if the value of \eqn{\lambda} corresponding to minimal IC is indicated on the plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{NULL} (no title).
}
}

\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{selectLambda}
}
\examples{
X <- dataGen(m=1,n=100,p=10,eps=0.2,bLength=4)$data[[1]]

sl <- selectLambda(X,k=2,method="ROSPCAg",lstep=0.1)
selectPlot(sl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
