# this R script should jumpstart CP for prolfqua analysis ;)

# libraries
library(prolfqua)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(prolfquapp)
library(seqinr)

# read in DIANN LONGformat
# prlfqapDiaNN <- prolfquapp::read_DIANN_output(diann.path = "../2289908_diannOutput/WU288688/out-2023-05-05/diann-output.tsv", isUniprot = "FALSE", fasta.file = "../prozor_proteinInference_fromDIANNoutput/Homo_sapiens.GRCh38.pep.all_noAsteriks.fasta")
# Error: 'diann_read_output' is not an exported object from 'namespace:prolfqua'

diannPrec <- read_tsv("../2289908_diannOutput/WU288688/out-2023-05-05/diann-output.tsv")
colnames(diannPrec)
keepCols <-  c("File.Name", "Run", "Protein.Group", "Protein.Ids", "PG.Quantity", "Stripped.Sequence", "Modified.Sequence", "Protein.Q.Value")

diannSlim <- diannPrec[,keepCols]

# read in mapping tables for prozorRazorProtein and ENSP2GeneNames
myProzorRazors <- read_tsv("../prozor_proteinInference_fromDIANNoutput/p2699_ENSP2PeptideMapping_prozor.txt")
myENSP2GeneSymbol <- read_tsv("../prozor_proteinInference_fromDIANNoutput/p2699_allFasta_ENSP2GeneNameMapping.txt")


# join in razor protein from prozor
colnames(myProzorRazors) <- c("RazorProtein", "Stripped.Sequence")
colnames(diannSlim)
diannSlimRazor <- left_join(x = diannSlim, y = myProzorRazors)

# join in GeneName and ENSG
colnames(diannSlimRazor)
colnames(myENSP2GeneSymbol)
diannSlim_ready <- left_join(x = diannSlimRazor, y = myENSP2GeneSymbol, join_by(RazorProtein == allENSP))


# annotation
dataset <- read_csv("BFdataset.txt")
dataset$Run <- gsub("^x|.d.zip$|.d$|.raw$|.mzML$", "", basename(gsub("\\\\", "/", dataset$`Relative Path`)))
# are we having the same? -> we only have a subset in the dataset! let's only use these
unique(diannSlim_ready$Run) %in% unique(dataset$Run)

longPepDatawAnno <- left_join(x = diannSlim_ready, y = dataset)

# get rid of all where we do not have annotation
longPepAllAnno <- longPepDatawAnno |>  filter(!is.na(`Grouping Var`))
longPepAllAnno$Condition <- longPepAllAnno$`Grouping Var`

colnames(longPepAllAnno)

longPepAllAnno <- longPepAllAnno |> select("Run", "Protein.Group","PG.Quantity","Stripped.Sequence","Modified.Sequence","Protein.Q.Value","RazorProtein","allENSG","Genesymbols","Name","Condition")
# parse subject
longPepAllAnno$Subject <- gsub(x = sapply(strsplit((longPepAllAnno$Run), split = "_"), function(x)x[4]),pattern = "-.*",  replacement = "")


# define analysis table for prolfqua
tablePep <- prolfqua::AnalysisTableAnnotation$new()
tablePep$set_response("PG.Quantity")
# define here the roll up!
tablePep$hierarchy[["protein_Id"]] = c("Genesymbols")
tablePep$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")

# define factors that are important
tablePep$fileName = "Run"
tablePep$factors[["Condition_"]] <- "Condition"
#tablePep$factors[["Subject_"]] <- "Subject"
tablePep$factorDepth <- 1
#tablePep$factorDepth <- 2 # if subject should be included!
 

# FILTERING
dim(longPepAllAnno)
myHardIntensityThreshold <- 0
data_r.filtered_pep <- longPepAllAnno  %>%  filter(PG.Quantity > myHardIntensityThreshold)
dim(data_r.filtered_pep)

# prolfqua internals

# create configuration prolfqua
Myconfig_pep <- prolfqua::AnalysisConfiguration$new(tablePep)

# witold helps troubleshooting
save(data_r.filtered_pep, Myconfig_pep, file = "RobjectsForTroubleShooting.RData")
load("RobjectsForTroubleShooting.RData")


# setup analysis with config and data
data_r.filtered_pep$PG.Quantity # here everything is still ok
prlfquaData_pep <- setup_analysis(data_r.filtered_pep , Myconfig_pep)
# wrong sampleNames
prlfquaData_pep$sampleName <- prlfquaData_pep$Run
prlfquaData_pep$PG.Quantity # suddenly here we get very often same values!! AND many NAs??






# prolfqua object
lfqdata_pep <- LFQData$new(prlfquaData_pep, Myconfig_pep)

# roll it up
transformed <- lfqdata_pep$get_Transformer()$intensity_array(log2)$lfq
Myaggregator <- transformed$get_Aggregator()
Myaggregator$medpolish() # why is this not working witold???
lfqdata_pep <- Myaggregator$lfq_agg

# generate plotter and summariser objects
datSmrz <- lfqdata_pep$get_Summariser()
lfqpl <- lfqdata_pep$get_Plotter()

# number cruncing
datSmrz$plot_hierarchy_counts_sample()

# see stats (e.g. violings for CVs)
st_pep <- lfqdata_pep$get_Stats()
st_pep$violin()
st_pep$stats()

# make a copy for different normalizations
lfqcopy_pep_rob <- lfqdata_pep$get_copy()


# rob scale
tr_pep_rob <- lfqcopy_pep_rob$get_Transformer()
transformed_Data_pep_lg2 <- tr_pep_rob$log2()
transformed_Data_pep_rob <- transformed_Data_pep_lg2$robscale()



