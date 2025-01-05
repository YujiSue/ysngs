args <- commandArgs(trailingOnly=TRUE)
bmpkgs <- installed.packages(repos = BiocManager::repositories()[["Bioconductor"]])
print(bmpkgs[bmpkgs[,"Package"]=="DESeq2", "Version"])
