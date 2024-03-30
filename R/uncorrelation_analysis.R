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
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot( nrow(bb$data) == 25780)
#' conf <- old2new(bb$config$clone(deep=TRUE))
#' data <- bb$data |> dplyr::ungroup()
#' dim(data)
#' dataI <- mark_decorelated(data, conf)
#'
mark_decorelated <- function(data , config, minCorrelation = 0.7){
  qvalFiltX <- data |>  dplyr::group_by_at(config$table$hierarchy_keys()[1]) |> tidyr::nest()
  qvalFiltX <- qvalFiltX |>
    dplyr::mutate(spreadMatrix = map(data, response_as_matrix, config))
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
#' bb <- prolfqua_data('data_ionstar')$normalized()
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
  nestedX <- nestedX |> dplyr::mutate(spreadMatrix = map(data, response_as_matrix, config))

  response_matrix_as_tibble <- function(x,config){
    x <- dplyr::bind_cols(
      row = rownames(x),
      tibble::as_tibble(x)
    )
    tidyr::gather(x,key = !!config$table$sampleName, value = "srm_ImputedIntensity", 2:ncol(x))
  }

  nestedX <- nestedX |> dplyr::mutate(imputed = map(.data$spreadMatrix, simpleImpute))

  nestedX <- nestedX |> dplyr::mutate(imputed = map(.data$imputed, response_matrix_as_tibble, config))
  unnest_res <- nestedX |> dplyr::select(config$table$hierarchy_keys_depth(), "imputed") |> tidyr::unnest(cols = .data$imputed)
  unnest_res <- unnest_res |> tidyr::separate("row",config$table$hierarchy_keys()[-1], sep = "~lfq~" )

  qvalFiltX <- dplyr::inner_join(x, unnest_res,
                                 by = c(config$table$hierarchy_keys(), config$table$sampleName) )
  config$table$set_response("srm_ImputedIntensity")
  return(qvalFiltX)
}


.make_name_AinB <- function(levelA, levelB, prefix="nr_"){
  c_name <- paste(prefix, levelB, "_IN_", levelA, sep = "")
  return(c_name)
}

.nr_B_in_A <- function(data,
                       levelA,
                       levelB,
                       merge = TRUE){
  namA <- paste(levelA, collapse = "_")
  namB <- paste(levelB, collapse = "_")
  c_name <- .make_name_AinB(namA, namB)
  if (!c_name %in% colnames(data) ) {
    data$c_name <- NULL
  }
  tmp <- data |>
    dplyr::select_at(c(levelA, levelB)) |>
    dplyr::distinct() |>
    dplyr::group_by_at(levelA) |>
    dplyr::summarize(!!c_name := n())

  if (!merge) {
    return(tmp)
  }
  data <- dplyr::inner_join(data, tmp, by = levelA )
  message("Column added : ", c_name)
  return(list(data = data, name = c_name))
}


#' Compute nr of B per A
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data |> dplyr::select(-all_of("nr_peptide_Id_IN_protein_Id"))
#' hierarchy <- config$table$hierarchy_keys()
#' res <- nr_B_in_A(data, config)
#'
#' res$data |>
#'   dplyr::select_at(c(config$table$hierarchy_keys_depth(),  res$name)) |>
#'   dplyr::distinct() |>
#'   dplyr::pull() |> table()
#'
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' config <- old2new(bb$config_f())
#' data <- bb$data
#' data$Area[data$Area == 0] <- NA
#' analysis <- setup_analysis(data, config)
#'
#' resDataStart <- bb$analysis(bb$data, config)
#'
#'
#' nr_B_in_A(resDataStart, config)
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#' config$table$hierarchyDepth <- 2
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#'
#' bb <- prolfqua_data('data_IonstarProtein_subsetNorm')
#' bb$config <- old2new(bb$config$clone(deep=TRUE))
#' nr_B_in_A(bb$data, bb$config)
#' #undebug(nr_B_in_A)
nr_B_in_A <- function(pdata, config , merge = TRUE){
  levelA <- config$table$hkeysDepth()
  levelB <- config$table$hierarchyKeys()[length(levelA) + 1]
  if (is.na(levelB)) {
    warning("here is no B in A")
    return(NULL)
  }else{
    .nr_B_in_A(pdata, levelA, levelB , merge = merge)
  }
}



#' how many peptides per protein in each sample
#' @export
#' @keywords internal
#' @family summary
#' @examples
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' nr_B_in_A_per_sample(data, configur, nested =FALSE)
#' bb <- prolfqua_data('data_IonstarProtein_subsetNorm')
#' bb$config <- old2new(config = bb$config$clone( deep = TRUE))
#' nr_B_in_A_per_sample(bb$data, bb$config, nested=FALSE)
#'
nr_B_in_A_per_sample <- function(data, config, nested = TRUE){
  cf <- config

  levelA <- cf$table$hierarchy_keys_depth()
  levelB <- cf$table$hierarchy_keys()[length(levelA) + 1]
  if (is.na(levelB)) {
    warning("here is no B in A")
  }
  data <- prolfqua::complete_cases(data, cf)
  data <- data |>
    dplyr::mutate(presentabsent = case_when(!is.na(!!sym(cf$table$get_response())) ~ 1,
                                            TRUE ~ 0))
  pepStats <- data |> group_by_at(c(cf$table$hierarchy_keys_depth(), cf$table$sampleName)) |>
    summarize(nrPep = n(), present = sum(.data$presentabsent), .groups = "drop")

  annotColumns <- c(cf$table$fileName,
                    cf$table$sampleName,
                    cf$table$hierarchy_keys_depth(),
                    cf$table$factor_keys_depth(),
                    cf$table$isotopeLabel)
  annotation <- data |>
    dplyr::select(!!!syms(annotColumns) ) |>
    distinct()

  res <- inner_join(annotation, pepStats, by = c(cf$table$sampleName, cf$table$hierarchy_keys_depth() ))
  res <-  if (nested) {res |> group_by_at(cf$table$hierarchy_keys_depth()) |> nest()} else {res}
  return(res)
}


# Summarize Intensities by Intensity or NAs ----
.rankProteinPrecursors <- function(data,
                                   config,
                                   column = config$table$get_response(),
                                   fun = function(x){ mean(x, na.rm = TRUE)},
                                   summaryColumn = "srm_meanInt",
                                   rankColumn = "srm_meanIntRank",
                                   rankFunction = function(x){ min_rank(desc(x)) }
){
  table <- config$table

  summaryPerPrecursor <- data |>
    dplyr::group_by(!!!syms(table$hierarchy_keys())) |>
    dplyr::summarize(!!summaryColumn := fun(!!sym(column)))

  groupedByProtein <- summaryPerPrecursor |>
    dplyr::arrange(!!sym( table$hierarchy_keys()[1])) |>
    dplyr::group_by(!!sym( table$hierarchy_keys()[1]))
  rankedBySummary <- groupedByProtein |>
    dplyr::mutate(!!rankColumn := rankFunction(!!sym(summaryColumn)))

  data <- dplyr::inner_join(data, rankedBySummary)
  return(data)
}

#' ranks precursor - peptide by intensity.
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' res <- remove_large_QValues(bb$data, bb$config)
#' res <- rank_peptide_by_intensity(res,bb$config)
#' X <-res |> dplyr::select(c(bb$config$table$hierarchy_keys(),
#'  srm_meanInt, srm_meanIntRank)) |> dplyr::distinct()
#' X |> dplyr::arrange(!!!rlang::syms(c(bb$config$table$hierarchy_keys()[1], "srm_meanIntRank"  )))
rank_peptide_by_intensity <- function(pdata, config){
  summaryColumn <- "srm_meanInt"
  rankColumn <- "srm_meanIntRank"
  pdata <- .rankProteinPrecursors(pdata, config, column = config$table$get_response(),
                                  fun = function(x){ mean(x, na.rm = TRUE)},
                                  summaryColumn = summaryColumn,
                                  rankColumn = rankColumn,
                                  rankFunction = function(x){min_rank(desc(x))}
  )

  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(pdata)
}


# Summarise NAs on lowest hierarchy ----

#' Ranks peptides/precursors of a protein by NAs (adds new column .NARank)
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- prolfqua_data('data_spectronautDIA250_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' res <- remove_large_QValues(analysis, config)
#' res <- rank_by_NA(res,config)
#' colnames(res)
#' x <- res |>
#'   dplyr::select(config$table$hierarchy_keys()[1], config$table$hierarchy_keys(TRUE)[1], "srm_NrNotNAs") |>
#'   dplyr::distinct() |> dplyr::summarize(sum(srm_NrNotNAs)) |> dplyr::pull()
#' stopifnot(sum(!is.na(res[[config$table$get_response()[1]]])) == x)
#' res |> dplyr::select(c(config$table$hierarchy_keys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) |>
#'  dplyr::distinct() |>
#'  dplyr::arrange(!!!rlang::syms(c(config$table$hierarchy_keys()[1],"srm_NrNotNARank")))
rank_by_NA <- function(pdata, config){
  summaryColumn <- "srm_NrNotNAs"
  rankColumn <- "srm_NrNotNARank"
  pdata <- .rankProteinPrecursors(pdata, config,
                                  column = config$table$get_response(),
                                  fun = function(x){sum(!is.na(x))},
                                  summaryColumn = summaryColumn,
                                  rankColumn = rankColumn,
                                  rankFunction = function(x){min_rank(desc(x))}
  )
  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(pdata)
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
#' bb <- prolfqua_data('data_spectronautDIA250_A')
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


