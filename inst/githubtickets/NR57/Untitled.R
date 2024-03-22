library(tibble)
library(dplyr)
library(prolfqua)

#read input annotation
#read input protein

source("../../../R/tidyMSFragger.R")

tmp <- read.csv("combined_protein.tsv",
         header = TRUE,
         sep = "\t",
         stringsAsFactors = FALSE,
         check.names = FALSE)

debug(tidy_FragPipe_combined_protein)
res <- tidy_FragPipe_combined_protein("combined_protein.tsv",
                               intnames = "Razor Intensity",
                               maxlfqnames = "MaxLFQ Razor Intensity")
