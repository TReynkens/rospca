\name{diagPlot}
\alias{diagPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic plot for PCA
}
\description{
Make diagnostic plot using the output from \code{robpca} or \code{rospca}.

}
\usage{
diagPlot(res, title = "Robust PCA", col = "black", pch = 16, labelOut = TRUE, id = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{A list containing the orthogonal distances (\code{od}),
  the score distances (\code{sd}) and their respective cut-offs (\code{cutoff.od} and \code{cutoff.sd}). Output from \code{robpca} or \code{rospca} can for example be used.}
  
  \item{title}{Title of the plot, default is \code{"Robust PCA"}.}
  \item{col}{Colour of the points in the plot, this can be a single colour for all points or a vector specifying the colour for each point. The default is \code{"black"}.}
  \item{pch}{Plotting characters or symbol used in the plot, see \link[graphics]{points} for more details. The default is 16 which corresponds to filled circles.}
   \item{labelOut}{Logical indicating if outliers should be labelled on the plot, default is \code{TRUE}.} 
   \item{id}{Number of OD outliers and number of SD outliers to label on the plot, default is 3.}

}
\details{
The diagnostic plot contains the score distances on the x-axis and the orthogonal distances on the y-axis.
To detect outliers, cut-offs for both distances are added, see Hubert et al. (2005).
}

\references{
Hubert, M., Rousseeuw, P. J., and Vanden Branden, K. (2005), ``ROBPCA: A New Approach to Robust Principal Component Analysis,'' \emph{Technometrics}, 47, 64--79.
}
\author{
Tom Reynkens with code from V. Todorov for \pkg{rrcov}.
}

\examples{
X <- dataGen(m=1,n=100,p=10,eps=0.2,bLength=4)$data[[1]]

resR <- robpca(X,k=2,skew=FALSE)
diagPlot(resR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ robust }