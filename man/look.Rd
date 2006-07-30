\name{look}
\alias{look}
\title{Take a first look at a dataset}
\description{
	This takes a first look at the multivariate structure of a dataset
}
\usage{
	look(x,forceType=TRUE)
}
\arguments{
  \item{x}{ a data frame}
  \item{forceType}{ whether to force non-numeric types for better exploration}
}
\details{
	This is for initial data exploration. Currently it uses PCA, parallel coordinates, and a heatmap of the correlation matrix.
}
\value{
	a firstlook object suitable for printing or plotting
}


\author{
Micah Altman
\email{Micah\_Altman\@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah\_altman/}
}


\examples{
	print(fl<-look(swiss))
	plot(fl)
}
\keyword{manip}
\keyword{misc}
\keyword{math}
