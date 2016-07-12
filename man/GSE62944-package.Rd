\name{GSE62944-package}

\alias{GSE62944-package}
\alias{GSE62944}

\title{GEO accession GSE62944 available as an \code{ExpressionSet} object.}

\description{
  TCGA re-processed RNA-Seq data from 9264 Tumor Samples and 741 normal 
  samples across 24 cancer types and made it available as GEO
  accession 
  [GSE62944](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944).
  These data have been parsed into a \code{SummazrixedExperiment object and
  are available in \code{ExperimentHub}.
}

\details{
  See the vignette for examples of using these data in differential
  gene expression analysis.

  browseVignettes("GSE62944")

  Details of how these data were creates are in the scripts/ directory
  of the source package.
}

\examples{
library(ExperimentHub)
hub <- ExperimentHub()
x <- query(hub, "GSE62944", "se", "tumor")
x
y <- query(hub, "GSE62944" "se", "normal")
y
\dontrun{
## download resource
se_tumor = x[[1]]
se_normal = y[[1]]
}}

\keyword{utilities}
