\name{Colon}
\alias{Colon}
\docType{data}
\title{
Gene expression of coon samples in samples by features format.
}
\description{
Two gene expression datasets merged together based on matched genes.
}
\usage{data("Colon")}
\format{
  The format is:
 num [1:112, 1:362] 8.68 6.14 6.85 6.47 4.52 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:112] "GSM452148" "GSM452149" "GSM452150" "GSM452151" ...
  ..$ : chr [1:362] "CCL5" "PRSS33" "FRAS1" "HS6ST2" ...
}
\details{
Two gene expression datasets merged together based on matched genes.
}
\source{
Nyamundanda, G. and Sadanandam, A., 2016. A novel and robust statistical method to diagnose and correct batch effects in large-scale genomic data.
}
\references{
Nyamundanda, G. and Sadanandam, A., 2016. A novel and robust statistical method to diagnose and correct batch effects in large-scale genomic data.
}
\examples{
data(Colon)
## maybe str(Colon) ; plot(Colon) ...
}
\keyword{datasets}
