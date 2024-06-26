% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_workflowFunctions.R
\name{workflow_correlation_preprocessing_protein_intensities}
\alias{workflow_correlation_preprocessing_protein_intensities}
\title{correlation preprocessing}
\usage{
workflow_correlation_preprocessing_protein_intensities(
  pdata,
  config,
  minCorrelation = 0.7
)
}
\arguments{
\item{pdata}{data.frame}

\item{config}{AnalysisConfiguration}

\item{minCorrelation}{correlation threshold default 0.7}
}
\description{
correlation preprocessing
}
\examples{



bb <- prolfqua::sim_lfq_data_peptide_config()
config <- bb$config
data <- bb$data

config$parameter$min_nr_of_notNA  <- 3
runLong <- TRUE
if(runLong){
  res <- workflow_correlation_preprocessing_protein_intensities(data,config)
  names(res)
}


}
\seealso{
Other workflows: 
\code{\link{workflow_DIA_NA_preprocessing}()},
\code{\link{workflow_corr_filter_impute}()}

Other deprecated: 
\code{\link{workflow_DIA_NA_preprocessing}()}
}
\concept{deprecated}
\concept{workflows}
\keyword{internal}
