args <- commandArgs(trailingOnly=TRUE)
if (args[1] %in% installed.packages(repos = BiocManager::repositories()[["Bioconductor"]])[, "Package"]) {
  print("TRUE")
} else {
  print("FALSE")
}