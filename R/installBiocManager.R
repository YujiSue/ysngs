if (!requireNamespace("BiocManager"))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
print(BiocManager::version())
BiocManager::install(c("RSQLite", "cummeRbund", "edgeR"))
