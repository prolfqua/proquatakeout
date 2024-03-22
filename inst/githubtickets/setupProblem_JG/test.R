library(prolfqua)
tmp <- load("TroubleshootThanksToWitold/RobjectsForTroubleShooting.RData")

data <- data_r.filtered_pep
hist(data$PG.Quantity)
class(Myconfig_pep)
xx <- setup_analysis( data, Myconfig_pep)
xx |> dplyr::filter(n > 1)
Myconfig_pep$table$workIntensity
