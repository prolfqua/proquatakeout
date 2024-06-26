% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uncorrelation_analysis.R
\name{mark_decorelated}
\alias{mark_decorelated}
\title{Marks peptides which do not correlate with other peptides of a protein}
\usage{
mark_decorelated(data, config, minCorrelation = 0.7)
}
\arguments{
\item{data}{data}

\item{config}{configuration}

\item{minCorrelation}{smallest allowed correlation}
}
\description{
It computes the pairwise correlation among all peptides and marks those
with with a corralation less then minCorrelation to all other peptides
}
\section{TODO}{
 do something with warnings of type "the standard deviation is zero".

 do investigate In max(x, na.rm = TRUE) : no non-missing arguments to max; returning -Inf
}

\examples{

bb <- prolfqua::sim_lfq_data_peptide_config()
conf <- bb$config
data <- bb$data |> dplyr::ungroup()
dataI <- mark_decorelated(data, conf)

}
\keyword{internal}
