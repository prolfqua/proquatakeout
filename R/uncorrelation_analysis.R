# uncorrelation analysis ----
.findDecorrelated <- function(res, threshold = 0.65){
  if (is.null(res))
    return(NULL)
  nrtrans <- ncol(res)
  ids <- rowSums(res < threshold, na.rm = TRUE)
  names(which((nrtrans - 1) == ids))
}

.decorelatedPly <- function(pdata, corThreshold = 0.7){
  res <- prolfqua::cor_jackknife_matrix(pdata)
  decorelated <- .findDecorrelated(res,threshold = corThreshold)
  tibble( row = rownames(res), srm_decorelated = rownames(res) %in% decorelated)
}

#' Marks peptides which do not correlate with other peptides of a protein
#'
#' It computes the pairwise correlation among all peptides and marks those
#' with with a corralation less then minCorrelation to all other peptides
#'
#' @export
#' @keywords internal
#' @param data data
#' @param config configuration
#' @param minCorrelation smallest allowed correlation
#' @section TODO: do something with warnings of type "the standard deviation is zero".
#' @section TODO: do investigate In max(x, na.rm = TRUE) : no non-missing arguments to max; returning -Inf
#' @examples
#'
#' bb <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#' stopifnot( nrow(bb$data) == 25780)
#' conf <- old2new(bb$config$clone(deep=TRUE))
#' data <- bb$data |> dplyr::ungroup()
#' dim(data)
#' dataI <- mark_decorelated(data, conf)
#'
mark_decorelated <- function(data , config, minCorrelation = 0.7){
  qvalFiltX <- data |>  dplyr::group_by_at(config$table$hierarchy_keys()[1]) |> tidyr::nest()
  qvalFiltX <- qvalFiltX |>
    dplyr::mutate(spreadMatrix = purrr::map(data, response_as_matrix, config))
  #HLfigs2 <- qvalFiltX |>
  #  dplyr::mutate(srmDecor = map(.data$spreadMatrix, .decorelatedPly,  minCorrelation))
  HLfigs2 <- qvalFiltX
  HLfigs2$srmDecor <- vector(mode = "list", length = nrow(qvalFiltX))
  for (i in seq_len(nrow(qvalFiltX))) {
    HLfigs2$srmDecor[[i]] <- .decorelatedPly(qvalFiltX$spreadMatrix[[i]], minCorrelation)
  }

  unnest_res <- HLfigs2 |>
    dplyr::select(config$table$hierarchy_keys()[1], "srmDecor") |> tidyr::unnest()
  unnest_res <- unnest_res |>
    tidyr::separate(col = "row",
                    into = config$table$hierarchy_keys()[-1],
                    sep = "~lfq~")
  qvalFiltX <- dplyr::inner_join(x = data, y = unnest_res, by = c(config$table$hierarchy_keys()) )
  return(qvalFiltX)
}


# Missing Value imputation ----

simpleImpute <- function(data){
  m <- apply(data,2, mean, na.rm = TRUE )
  res <- sweep(data,2,m,"-")
  dim(data)
  dim(res)
  resMean <- apply(res, 1, mean, na.rm = TRUE)
  resid <- matrix(replicate(length(m),resMean), nrow = length(resMean))
  imp <- sweep(resid,2,m,"+")
  res <- data
  res[is.na(res)] <- imp[is.na(res)]
  return(res)
}

#' Imputate missing peptide intensities to maximize peptide correlationion
#'
#' Assumes the peptide intensities are correlation assumption
#'
#' @export
#' @param x data
#' @param config configuration
#' @keywords internal
#' @examples
#'
#'
#'
#' bb <- prolfqua::prolfqua_data('data_ionstar')$normalized()
#' bb$config <- old2new(bb$config)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' mean(is.na(data$peptide.intensity))
#' dataI <- impute_correlationBased(data, config)
#' dim(dataI)
#' stopifnot(dim(dataI) == c(dim(data)+c(0,1)))
#' stopifnot(mean(is.na(dataI$srm_ImputedIntensity)) <= mean(is.na(data$peptide.intensity)))
#'
impute_correlationBased <- function(x , config){
  x <- complete_cases(x, config)
  nestedX <- x |>  dplyr::group_by_at(config$table$hierarchy_keys_depth()) |> tidyr::nest()
  nestedX <- nestedX |> dplyr::mutate(spreadMatrix = purrr::map(data, response_as_matrix, config))

  response_matrix_as_tibble <- function(x,config){
    x <- dplyr::bind_cols(
      row = rownames(x),
      tibble::as_tibble(x)
    )
    tidyr::gather(x,key = !!config$table$sampleName, value = "srm_ImputedIntensity", 2:ncol(x))
  }

  nestedX <- nestedX |> dplyr::mutate(imputed = purrr::map(.data$spreadMatrix, simpleImpute))

  nestedX <- nestedX |> dplyr::mutate(imputed = purrr::map(.data$imputed, response_matrix_as_tibble, config))
  unnest_res <- nestedX |> dplyr::select(config$table$hierarchy_keys_depth(), "imputed") |> tidyr::unnest(cols = .data$imputed)
  unnest_res <- unnest_res |> tidyr::separate("row",config$table$hierarchy_keys()[-1], sep = "~lfq~" )

  qvalFiltX <- dplyr::inner_join(x, unnest_res,
                                 by = c(config$table$hierarchy_keys(), config$table$sampleName) )
  config$table$set_response("srm_ImputedIntensity")
  return(qvalFiltX)
}





#' Removes measurments with more than some percentage of missing values
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param percent percent of missing values
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#'
#' bb <- prolfqua::prolfqua_data('data_spectronautDIA250_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- analysis
#' data <- remove_large_QValues(data, config)
#' hc <- hierarchy_counts(data, config)
#' res <- filter_factor_levels_by_missing(data, config,percent = 80)
#' hc80 <-  hierarchy_counts(res, config)
#' res <- filter_factor_levels_by_missing(data, config,percent = 60)
#' hierarchy_counts(res, config)
#' hc60 <-  hierarchy_counts(res, config)
#' stopifnot(all(hc60 >= hc80)) # since 80% missing is more stringent than 60%
#' stopifnot(all(hc >= hc60))
#' summarize_hierarchy(res,config) |>
#'  dplyr::filter(!!rlang::sym(paste0(config$table$hierarchy_keys()[2],"_n")) > 1)
#'
filter_factor_levels_by_missing <- function(pdata,
                                            config,
                                            percent = 60){
  table <- config$table
  summaryColumn = "srm_NrNotNAs"
  column <- table$get_response()

  pdata <- complete_cases( pdata , config)
  nrNA = function(x){sum(!is.na(x))}
  summaryPerPrecursor <- pdata |>
    dplyr::group_by(!!!syms( c(table$hierarchy_keys(), table$factor_keys_depth() ))) |>
    dplyr::summarize(!!"nr" := n(), !!summaryColumn := nrNA(!!sym(column))) |>
    dplyr::mutate(fraction = !!sym(summaryColumn)/!!sym("nr") * 100 ) |>  dplyr::ungroup()

  summaryPerPrecursorFiltered <- summaryPerPrecursor |> dplyr::filter(.data$fraction > percent)
  summaryPerPrecursorFiltered <- summaryPerPrecursorFiltered |>
    dplyr::select(c(table$hierarchy_keys())) |> dplyr::distinct()
  stopifnot(all(colnames(summaryPerPrecursorFiltered) %in% table$hierarchy_keys()))
  res <- summaryPerPrecursorFiltered |> left_join(pdata)
  return( dplyr::ungroup(res))
}


