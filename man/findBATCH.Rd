\name{findBATCH}
\alias{findBATCH}
\title{
A function to formally tests presence of batch effect.
}
\description{
It searches for directions in the principal subspace that are associated with the batch variable.
}
\usage{
findBATCH(res2, res1, nt, rerun)
}
\arguments{
  \item{res2}{
Results of running findBATCH to assess batch effect.
}
  \item{res1}{
Results of the PPCCA model with the optimal number of PCs
}
  \item{nt}{
Number of biological classes.
}
  \item{rerun}{
Results of running findBATCH to assess biological effect.
}
}
\details{
It searches for directions in the principal subspace that are associated with the batch variable.
}
\value{
  Results in findBATCH folder in the working directory.
}
\references{
Nyamundanda, G. and Sadanandam, A., 2016. A statistical method for evaluating and correcting batch effects in genomic data.
}
\author{
Gift Nyamundanda and Anguraj Sadanandam
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
