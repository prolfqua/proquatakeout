

#' correlation preprocessing
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param minCorrelation correlation threshold default 0.7
#' @export
#' @keywords internal
#' @family workflows
#' @family deprecated
#' @examples
#'
#'
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' config <- bb$config
#' data <- bb$data
#'
#' config$parameter$min_nr_of_notNA  <- 3
#' runLong <- TRUE
#' if(runLong){
#'   res <- workflow_correlation_preprocessing_protein_intensities(data,config)
#'   names(res)
#' }
#'
#'
workflow_correlation_preprocessing_protein_intensities <- function(
    pdata, config, minCorrelation = 0.7){
  stat_input <- prolfqua::hierarchy_counts(pdata, config)

  data_NA_QVal <- prolfqua::filter_byQValue(pdata, config)
  stat_qval <- prolfqua::hierarchy_counts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- prolfqua::rank_by_NA(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal |> dplyr::filter(.data$srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  if (nrow(data_NA_QVal) == 0) {
    warning("no rows left after filtering for min_nr_of_notNA")
  }
  stat_min_nr_of_notNA <- prolfqua::hierarchy_counts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- prolfqua::filter_proteins_by_peptide_count(data_NA_QVal, config)$data
  stat_min_peptides_protein  <- prolfqua::hierarchy_counts(data_NA_QVal, config)
  # filter decorrelated.
  data_NA_QVal <- prolfqua::transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- mark_decorelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, .data$srm_decorelated == FALSE)

  stat_correlated  <- prolfqua::hierarchy_counts(keepCorrelated, config)

  # TODO check if you are not aggregating log transformed intensities
  # rank precursors by intensity
  keepCorrelated <- prolfqua::rank_peptide_by_intensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)
  mean_na <- function(x, name=FALSE) {if (name) {return("mean_na")};mean(x, na.rm = TRUE)}
  proteinIntensities <- prolfqua::aggregate_intensity_topN(qvalFiltImputed, config, .func = mean_na, N = 3)

  # collect stats
  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- dplyr::bind_rows(stats)
  stats <- tibble::add_column(x, processing = names(stats),.before = 1)


  return(list(data = proteinIntensities$data, stats = stats, newconfig = proteinIntensities$newconfig))
}

#' Apply correlation filtering and impute missing values
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param minCorrelation correlation threshold default 0.7
#' @keywords internal
#' @family workflows
#' @export
#' @examples
#'
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' bb$config <- bb$config
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' data$nr_peptide_Id_IN_protein_Id <- NULL
#'
#' config$parameter$min_nr_of_notNA  <- 3
#' res <- workflow_corr_filter_impute(data,config)
#'
workflow_corr_filter_impute <- function(pdata, config, minCorrelation =0.6){
  stat_input <- prolfqua::hierarchy_counts(pdata, config)

  data_NA_QVal <- prolfqua::filter_byQValue(pdata, config)
  stat_qval <- prolfqua::hierarchy_counts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- prolfqua::rank_by_NA(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal |> dplyr::filter(.data$srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  stat_min_nr_of_notNA <- prolfqua::hierarchy_counts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- prolfqua::filter_proteins_by_peptide_count(data_NA_QVal, config)$data

  stat_min_peptides_protein  <- prolfqua::hierarchy_counts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- prolfqua::transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- mark_decorelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, .data$srm_decorelated == FALSE)

  stat_correlated  <- prolfqua::hierarchy_counts(keepCorrelated, config)
  keepCorrelated <- prolfqua::rank_peptide_by_intensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)


  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- dplyr::bind_rows(stats)
  stats <- tibble::add_column(x, processing = names(stats),.before = 1)
  return(qvalFiltImputed)
}

#' filter QVAlues and NA's and factor information
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param minCorrelation correlation threshold default 0.7
#' @export
#' @keywords internal
#' @family workflows
#' @family deprecated
#' @examples
#'
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' data$nr_peptide_Id_IN_protein_Id <- NULL
#'
#' prolfqua::hierarchy_counts(data, config)
#' tmp <- workflow_DIA_NA_preprocessing(data, config)
#' prolfqua::hierarchy_counts(tmp$data, config)
#' tmp <- workflow_DIA_NA_preprocessing(data, config, percent=70)
#' prolfqua::hierarchy_counts(tmp$data, config)
#' stopifnot(FALSE == (dplyr::is_grouped_df(tmp$data)))
#'
workflow_DIA_NA_preprocessing <- function(pdata,
                                          config,
                                          percent = 60,
                                          hierarchy_level = 2,
                                          factor_level = 1,
                                          min_peptides_protein = config$parameter$min_peptides_protein)
{
  stat_input <- prolfqua::hierarchy_counts(pdata, config)

  data_NA_QVal <- prolfqua::filter_byQValue(pdata, config)
  stat_qval <- prolfqua::hierarchy_counts(data_NA_QVal, config)

  resNACondition <- filter_factor_levels_by_missing(data_NA_QVal,
                                                    config,
                                                    percent = percent)

  stat_naFilter <- prolfqua::hierarchy_counts(resNACondition, config)
  protID <- prolfqua::summarize_hierarchy(resNACondition,config) |>
    dplyr::filter(!!rlang::sym(paste0(config$table$hierarchy_keys()[hierarchy_level],"_n"))
                  >= min_peptides_protein)

  data_NA_QVal_condition <- protID |>
    dplyr::select(config$table$hierarchy_keys()[1]) |>
    dplyr::inner_join(resNACondition)

  # Complete cases
  data_NA_QVal_condition <- prolfqua::complete_cases( data_NA_QVal_condition , config)
  stat_peptidFitler <- prolfqua::hierarchy_counts(data_NA_QVal_condition, config)
  stats = list(stat_input = stat_input,
               stat_qval = stat_qval,
               stat_naFilter = stat_naFilter,
               stat_peptidFitler = stat_peptidFitler
  )
  return(list(data = data_NA_QVal_condition, stats = stats))
}

