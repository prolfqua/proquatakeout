#' copy all files need to run mixed model analysis.
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#' @examples
#' copy_mixed_model_analysis_script(workdir = tempdir())
copy_mixed_model_analysis_script <- function(workdir = getwd()){
  runscripts <- c("fgcz_formatting/fgcz_header.html",
                  "fgcz_formatting/fgcz_footer.html",
                  "fgcz_formatting/fgcz.css",
                  "fgcz_formatting/fgcz_banner.png",
                  "rmarkdown/mixed_model_analysis_script_Report.Rmd",
                  "rmarkdown/bibliography.bib"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir)
}


#' render Filtering Summary.
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @family vignetteHelpers
#' @keywords internal
#'
#' @export
#'
render_SummarizeFiltering_rmd <-
  function(results,
           dest_path = ".",
           dest_file_name = "Summarize_Filtering.pdf",
           workdir = tempdir())
  {
    dist_file_path <- .run_markdown_with_params(
      results,
      markdown_path = "rmarkdown/Summarize_Filtering.Rmd",
      dest_path = dest_path,
      dest_file_name = dest_file_name,
      workdir = workdir,
      packagename = "prolfqua"
    )
  }
