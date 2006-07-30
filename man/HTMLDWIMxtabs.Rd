\name{HTML.DWIMxtabs}
\alias{HTML.DWIMxtabs}
\title{ pretty print DWIMxtabs}
\description{
	This is used with R2HTML to pretty-print DWIMtabs
}
\usage{
	HTML.DWIMxtabs(x,...,
	 wantPercentages=T, wantTotals=T, wantStats=T, nameLength=15)
}
\arguments{
  \item{x}{ object to converted to HTML}
  \item{\dots}{ additional parameters passed to HTML}
  \item{wantPercentages} { include percentages in the output }
  \item{wantTotals} { include marginals in the output }
  \item{wantStats} { include chi-squared statistics in the output }
  \item{nameLength} { maximum length of variable labels}

}

\value{
	None
}
\examples{
	vxt=DWIMxtabs(Freq~.,as.data.frame(Titanic))
\dontrun{
	# HTML methods called implicitly
	HTML(vxt,file="")
}

}
\author{
Micah Altman
\email{micah\_altman\@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah\_altman}
}


\seealso{\code{\link[pkg:R2HTML]{HTML}}, \link{DWIMxtabs}}}}

\keyword{iplot}
\keyword{IO}
\keyword{print}
