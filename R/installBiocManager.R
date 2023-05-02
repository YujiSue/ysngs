args <- commandArgs(trailingOnly = T)
repos_url=""
if (1 < length(args)) { repos_url=args[1] } else { repos_url="http://cran.us.r-project.org" }
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager", repos = repos_url)
