args <- commandArgs(trailingOnly=TRUE)
bmpkgs <- installed.packages(repos = BiocManager::repositories()[["Bioconductor"]])
print(bmpkgs[bmpkgs[,"Package"]==args[1], "Version"][1])
