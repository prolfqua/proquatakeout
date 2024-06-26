#!/usr/bin/Rscript

suppressMessages(library(readr))
suppressMessages(library(tidyverse))
suppressMessages(library(prolfqua))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))

doc <- "Sample Size Report from MQ file

Usage:
  lfq_MQ_SampleSizeReport.R <mqzip> [--outdir=<outdir>] [--project_Id=<project_Id>] [--order_Id=<order_Id>] [--workunit=<workunit>]

Options:
  -o --outdir=<outdir> output directory [default: samplesize_qc_report]
  -p --project_Id=<project_Id> project identifier [default: NULL]
  -q --order_Id=<order_Id> order identifier [default: NULL]
  -w --workunit=<workunit> parent workunit [default: NULL]

Arguments:
  mqzip input file
"

if (TRUE) {

  dummyargs <- c( "c:/Users/wewol/prog/prolfquaTestDataExtern/1675194.zip" )
  dummyargs <- c("D:/TEMP/checkSampleSizeEstimation/1293562.zip","-p","3000","-q","9999032")
  opt <- docopt(doc, args = dummyargs)
}else{
  opt <- docopt(doc)
}

cat("\nParameters used:\n\t",
    "     mqzip:", mqzip <- opt$mqzip, "\n\t",
    "result_dir:", outputDir <- opt[["--outdir"]], "\n\t",
    "project_Id:", project_Id <- opt[["--project_Id"]], "\n\t",
    "  order_Id:", order_Id <- opt[["--order_Id"]], "\n\t",
    "  workunit:", workunit <- opt[["--workunit"]], "\n\n\n")



if (!dir.exists(outputDir)) {
  if (!dir.create(outputDir)) {
    stop("\n could not create : ", outputDir, "\n")
  }
}

summarize_stats_raw_transformed <- function(resDataStart, config){
  config_tmp <- config$clone(deep = TRUE)
  stats_res <- summarize_stats(resDataStart, config_tmp)
  wide <- tidy_to_wide_config( resDataStart,config_tmp)
  stats_raw <- inner_join(stats_res,
                          wide$data,
                          by = config_tmp$table$hierarchy_keys())

  resDataStart <- transform_work_intensity(resDataStart, config_tmp, log2)
  data <- prolfqua::apply_to_response_matrix(resDataStart, config_tmp, .func = robust_scale)

  stats_res_transformed <- summarize_stats(data, config_tmp)
  wide <- tidy_to_wide_config( data,config_tmp)
  stats_transformed <- inner_join(stats_res_transformed,
                                  wide$data,
                                  by = config_tmp$table$hierarchy_keys())
  peptideStats <- inner_join(stats_raw, stats_transformed,
                             by = c( config_tmp$table$factor_keys(),config_tmp$table$hierarchy_keys() ),
                             suffix = c(".raw",".transformed") )
  return(peptideStats)
}


resPep <- tidyMQ_Peptides(mqzip)

config <- prolfqua::create_config_MQ_peptide()
config$table$factors[["FACTOR"]] = "QC"
config$table$factorDepth <- 1

resPep$QC <- "QC"
resPep$isotope <- "light"

project_conf <- list()
project_conf$project_Id <- project_Id
project_conf$order_Id <- order_Id
project_conf$workunit_Id <- workunit


resPep <- setup_analysis(resPep, config)
resPep <- remove_small_intensities(resPep, config) |>
  complete_cases(config)
resPepL <- list()
resPepL$data <- resPep
resPepL$config <- config
filename <- basename(tools::file_path_sans_ext(mqzip))

# debug(render_MQSummary_rmd)
prolfqua::render_MQSummary_rmd(
  resPep$data,
  resPep$config$clone(deep = TRUE),
  project_conf,
  pep = TRUE,
  dest_path = outputDir,
  dest_file_name = paste0("r_", filename, "_Peptide"),
  workdir = ".",
  format = "html"
)

peptideStats <- summarize_stats_raw_transformed(resPep$data, resPep$config)
lfq_write_table(separate_hierarchy(peptideStats, resPep$config),
                path = outputDir,
                name = paste0("r_", filename, "_PeptideStats"), format = "xlsx")


# Perform filtering and Protein Aggregation...


path <- ""
res <- filter_proteins_by_peptide_count(resPep$data, resPep$config)
results <- normalize_log2_robscale(res$data, resPep$config)


#rmarkdown::render("MQSummary2.Rmd", params=list(configuration = config$clone(deep=TRUE), data = resDataStart), output_format = bookdown::pdf_document2())
protintensity <- medpolish_protein_estimates( results$data, results$config )
xx <- protintensity("unnest")
data <- prolfqua::apply_to_response_matrix(xx$data, xx$config, .func = robust_scale)


stats_res <- summarize_stats(data, xx$config, all = FALSE)


data_wide <-
  inner_join(
    stats_res,
    tidy_to_wide_config(data, xx$config)$data,
    by = xx$config$table$hierarchy_keys_depth(),
    suffix = c(".factor", "")
  )


lfq_write_table(
  separate_hierarchy(data_wide, xx$config),
  path = outputDir,
  name = paste0("r_", filename, "_Protein"),
  format = "xlsx"
)

lfq_write_table(
  tidy_to_wide_config(data, xx$config)$annotation,
  path = outputDir,
  name = paste0("r_", filename, "_annotation"),
  format = "xlsx"
)

##```{r corrplot, fig.cap="Correlation plot."}
##prolfqua::plot_sample_correlation(dataTransformed, config)
##```

prolfqua::render_MQSummary_rmd(
  data,
  xx$config$clone(deep = TRUE),
  project_conf = project_conf,
  pep = FALSE,
  dest_path = outputDir,
  dest_file_name = paste0("r_", filename, "_Protein"),
  workdir = ".",
  format = "html"
)

#rmarkdown::render("MQSummary2.Rmd", params=list(data =xx$data, configuration = xx$config), output_file = paste0("Protein_",mqzip,".pdf"))#,envir =new.env())
#break()
#}

