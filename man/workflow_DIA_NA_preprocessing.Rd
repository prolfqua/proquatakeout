% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_workflowFunctions.R
\name{workflow_DIA_NA_preprocessing}
\alias{workflow_DIA_NA_preprocessing}
\title{filter QVAlues and NA's and factor information}
\usage{
workflow_DIA_NA_preprocessing(
  pdata,
  config,
  percent = 60,
  hierarchy_level = 2,
  factor_level = 1,
  min_peptides_protein = config$parameter$min_peptides_protein
)
}
\arguments{
\item{pdata}{data.frame}

\item{config}{AnalysisConfiguration}

\item{minCorrelation}{correlation threshold default 0.7}
}
\description{
filter QVAlues and NA's and factor information
}
\examples{


bb <- prolfqua::sim_lfq_data_peptide_config()
config <- bb$config$clone(deep=TRUE)
data <- bb$data
data$nr_peptide_Id_IN_protein_Id <- NULL

prolfqua::hierarchy_counts(data, config)
tmp <- workflow_DIA_NA_preprocessing(data, config)
prolfqua::hierarchy_counts(tmp$data, config)
tmp <- workflow_DIA_NA_preprocessing(data, config, percent=70)
prolfqua::hierarchy_counts(tmp$data, config)
stopifnot(FALSE == (dplyr::is_grouped_df(tmp$data)))

}
\seealso{
Other workflows: 
\code{\link{workflow_corr_filter_impute}()},
\code{\link{workflow_correlation_preprocessing_protein_intensities}()}

Other deprecated: 
\code{\link{workflow_correlation_preprocessing_protein_intensities}()}
}
\concept{deprecated}
\concept{workflows}
\keyword{internal}
