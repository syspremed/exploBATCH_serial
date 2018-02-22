\name{assessCorrectBATCH}
\alias{assessCorrectBATCH}
\title{
Assessing data corrected for batch effect using correctBATCH.
}
\description{
Assess if batch effect still exists after correctBATCH batch correction.
}
\usage{
assessCorrectBATCH(rerun2, res1, nt, rerun22)
}
\arguments{
  \item{rerun2}{
Results of running findBATCH to assess batch effect.
}
  \item{res1}{
Selection of the optimal number of PCs.
}
  \item{nt}{
Number of biological classes.
}
  \item{rerun22}{
Results of running findBATCH to assess biological effect.
}
}
\details{
Assess if batch effect still exists after correctBATCH batch correction.
}
\value{
  Results in assessCorrectBATCH folder in the working directory.
}
\references{
Nyamundanda, G. and Sadanandam, A., 2016. A novel and robust statistical method to diagnose and correct batch effects in large-scale genomic data.
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
