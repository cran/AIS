\name{print.DWIMxtabs}
\alias{print.DWIMxtabs}
\title{ print DWIMxtabs}
\description{
	This is used to pretty-print DWIMxtabs
}
\usage{
	print.DWIMxtabs(x,...,
	 wantPercentages=T, wantTotals=T, wantStats=T, nameLength=15 )
}
\arguments{
  \item{x}{ object to printed}
  \item{\dots}{ additional parameters passed to print.default}
  \item{wantPercentages} { include percentages in the output }
  \item{wantTotals} { include marginals in the output }
  \item{wantStats} { include chi-squared statistics in the output }
  \item{nameLength} { maximum length of variable labels}

}

\value{
	None
}
\examples{
	# print methods called implicitly
	DWIMxtabs(Freq~.,as.data.frame(Titanic))

}
\author{
Micah Altman
\email{micah\_altman\@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah\_altman}
}


\seealso{\link{DWIMxtabs}}}}

\keyword{iplot}
\keyword{IO}
\keyword{print}
