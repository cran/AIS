\name{DWIMxtabs}
\alias{DWIMxtabs}
\title{DWIM cross tabs, with pretty printing}
\description{
	This calculate cross tabs with marginals and proportions. 
}
\usage{
DWIMxtabs(formula,data=parent.env(),var.labels=NULL, ...) 
}

\arguments{
  \item{data}{ data to be analyzed }
  \item{formula}{ crosstab formula suitable for \code{xtabs} }
  \item{var.labels}{optional substitute varilable lables, names will be used if not supplied}
  \item{\dots}{arguments to pass to \code{xtabs}}
}
\value{
	Returns a crosstabulation, essentially a combination of \code{xtabs} with proportions and marginals. DWIMxtabs will automatcally recode continuous variables to categories prior to crosstabulation
}
\examples{
	DWIMxtabs(Freq~.,as.data.frame(Titanic))

}
\author{
Micah Altman
\email{micah\_altman\@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah\_altman}
}


\seealso{
\code{\link{xtabs}}, 
}

\keyword{IO}
\keyword{print}
