# CummeRbund




# DESeq2
deseqAnalyze <- function(data, group, formula, method = 'auto', reduced = as.formula('~1')) {
    dds <- DESeqDataSetFromMatrix(countData = data, colData = group, design = formula)
    if (method == 'auto') {
        dds <- DESeq(dds)
    }
    else if (method == 'lrt') {
        dds <- DESeq(dds, test = "LRT", reduced = reduced)
        #dds <- estimateSizeFactors(dds)
        #dds <- estimateDispersions(dds)
        #dds <- nbinomLRT(dds, full = formula, reduced = reduced)
    }
    return(dds)                 
}

deseqResult <- function(data, target, export_csv=T, suffix='') {
    # Get result
    res <- results(data, contrast = target)
    # Export data as CSV
    if (export_csv) {
        path = "DEG_result.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(res, file = path)
    }
    return(res)
}

# EdgeR



# Common
geneConversion <- function(genes, db, src, dest, select = 'first') {
    return(mapIds(db, keys = genes, keytype = src, column = dest, multiVals = select))
}

countData <- function(data, restype = 'edgeR', export_csv=T, suffix='') {
    # Get count
    ncounts <- NA
    if (restype == 'edgeR') {}
    else if (restype == 'DESeq2') {
        ncounts <- counts(data, normalized = TRUE)
    }
    # Export data as CSV
    if (export_csv) {
        # Make dir.
        path = "counts.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(ncounts, file = path)
    }
    return(ncounts)
}

zScoreData <- function(counts, group, factors,
                        collabel = NA,
                        rowlabel = NA,
                        cluster = T, 
                        export_csv=T, suffix='') {
    means <- list()
    for (f in factors) {
        means[[f]] <- rowMeans(counts[, group == f])
    }
    zscores <- do.call(cbind, means)
    m <- apply(zscores, 1, mean)
    s <- apply(zscores, 1, sd)
    zscores <- (zscores - m) / s
    #
    if (collabel) colnames(zscores) <- collabel
    else colnames(zscores) <- factors
    if (rowlabel) rownames(zscores) <- rowlabel
    else rownames(zscores) <- rownames(counts)
    #
    if (cluster) {
        d <- dist(scores)
        hc <- hclust(d)
        zscores <- zscores[hc$order,]
    }
    if (export_csv) {
        path = "zscores.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(ncounts, file = path)
    }
    return(zscores)
}

expressionChanged <- function(result, 
                              restype = 'edgeR',
                              fc_limit = c(1, Inf ),
                              p_threshold = 0.05, 
                              q_threshold = 0.05, 
                              export_csv=T, suffix='') {
    # Preset container
    up <- NA
    down <- NA
    if (restype == 'edgeR') {}
    else if (restype == 'DESeq2') {
        ## Up-regulated
        upfilter <- !is.na(result$padj) & 
                result$pvalue < p_threshold & 
                result$padj < q_threshold & 
                result$log2FoldChange >= fc_limit[1] &
                result$log2FoldChange <= fc_limit[2]
        up <- result[upfilter,]
        up <- up[order(up$log2FoldChange, decreasing = T),]
        
        ## Down-regulated
        downfilter <- !is.na(result$padj) & 
                result$pvalue < p_threshold & 
                result$padj < q_threshold & 
                result$log2FoldChange <= -fc_limit[1] &
                result$log2FoldChange >= -fc_limit[2]
        down <- res[downfilter,]
        down <- down[order(down$log2FoldChange, decreasing = F),]
    }
    #
    if (export_csv) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'DEG')
        dir.create(dir, showWarnings = FALSE)
        setwd(dir)
        ## 
        path <- "up.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(up, file = path)
        path <- "down.csv" 
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(down, file = path)
        ##
        setwd(parent)
    }
    return(c(up=up, down=down))
}

enrichmentGO <- function(data,
                         org,
                         genetype = 'ENTREZID',
                         gotype = c('CC', 'MF', 'BP'),
                         q_threshold = 0.05, 
                         export_csv=T,
                         export_img=T,
                         img_width = 8,
                         img_height = 6,
                         dot_plot=T,
                         ont_network=T,
                         suffix='') {
    
    # Preset gene list
    upgenes <- rownames(data$up)
    downgenes <- rownames(data$down)
    ## Gene ID conversion
    if (genetype != "ENTREZID") {
        upgenes <- geneConversion(gene = upgenes, db = org, src = genetype, dest = "ENTREZID")
        downgenes <- geneConversion(gene = downgenes, db = org, src = genetype, dest = "ENTREZID")
    }
    # GO enrichment
    ## CC:Cellular Component, MF:Molecular Function, BP:Biological Process   
    for (t in gotype) {
        go_up <- enrichGO(gene = upgenes,
                 OrgDb = org,
                 ont = t,
                 pAdjustMethod = "BH",
                 qvalueCutoff = q_threshold)
        
        go_down <- enrichGO(gene = downgenes,
                 OrgDb = org,
                 ont = t,
                 pAdjustMethod = "BH",
                 qvalueCutoff = q_threshold)
        
        #
        if (export_csv) {
            # Make dir.
            parent <- getwd()
            dir <- file.path(parent, 'GO')
            dir.create(dir, showWarnings = FALSE)
            setwd(dir)
            ## 
            path <- paste("GO", t, "up.csv", sep="_")
            if (suffix != '') path = paste(suffix, path, sep="_")
            write.csv(go_up, file = path)
            path <- paste("GO", t, "down.csv", sep="_")
            if (suffix != '') path = paste(suffix, path, sep="_")
            write.csv(go_down, file = path)
            ##
            setwd(parent)
        }
        #
        if (export_img) {
            # Make dir.
            parent <- getwd()
            dir <- file.path(parent, 'Chart')
            dir.create(dir, showWarnings = FALSE)
            setwd(dir)
            ##
            if (dot_plot) {
                ##                
                up_dot <- dotplot(go_up)
                path <- paste("GO", t, "up", "dot", sep="_")
                if (suffix != '') path = paste(suffix, path, sep="_")
                exportImage(up_dot, path, "png", width = img_width, height = img_height, resolution=150)
                ##
                down_dot <- dotplot(go_down)
                path <- paste("GO", t, "down", "dot", sep="_")
                if (suffix != '') path = paste(suffix, path, sep="_")
                exportImage(down_dot, path, "png", width = img_width, height = img_height, resolution=150)
            }
            ##
            if (ont_network) {
                ##
                up_net <- goplot(go_up)
                path <- paste("GO", t, "up", "net", sep="_")
                if (suffix != '') path = paste(suffix, path, sep="_")
                exportImage(up_net, path, "png", width = img_width, height = img_height, resolution=150)
                ##
                down_net <- goplot(go_down)
                path <- paste("GO", t, "down", "net", sep="_")
                if (suffix != '') path = paste(suffix, path, sep="_")
                exportImage(down_net, path, "png", width = img_width, height = img_height, resolution=150)
            }
            setwd(parent)
        }
    }   
}

enrichmentPathway <- function(data,
                              org,
                              genetype = 'ENTREZID',
                              target = 'KEGG',
                              gs_size = c(10, 500),
                              q_threshold = 0.05, 
                              export_csv=T,
                              export_img=T,
                              img_width = 8,
                              img_height = 6,
                              suffix='') {
                             
    # Preset gene list
    upgenes <- rownames(data$up)
    downgenes <- rownames(data$down)
    ## Gene ID conversion
    if (genetype != "ENTREZID") {
        upgenes <- geneConversion(gene = upgenes, db = org, src = genetype, dest = "ENTREZID")
        downgenes <- geneConversion(gene = downgenes, db = org, src = genetype, dest = "ENTREZID")
    }
    
    # KEGG enrichment
    if (target == 'KEGG') {
        kegg_up <- enrichKEGG(
            gene = upgenes,
            minGSSize = gs_size[1],
            maxGSSize = gs_size[2],
            pAdjustMethod = "BH",
            qvalueCutoff = q_threshold)
        
        kegg_down <- enrichKEGG(
            gene = downgenes,
            minGSSize = gs_size[1],
            maxGSSize = gs_size[2],
            pAdjustMethod = "BH",
            qvalueCutoff = q_threshold)
        
        #
        if (export_csv) {
            # Make dir.
            parent <- getwd()
            dir <- file.path(parent, 'KEGG')
            dir.create(dir, showWarnings = FALSE)
            setwd(dir)
            ## 
            path <- paste("KEGG", "up.csv", sep="_")
            if (suffix != '') path = paste(suffix, path, sep="_")
            write.csv(kegg_up, file = path)
            path <- paste("KEGG", "down.csv", sep="_")
            if (suffix != '') path = paste(suffix, path, sep="_")
            write.csv(kegg_down, file = path)
            ##
            setwd(parent)
        }
        #
        if (export_img) {
            # Make dir.
            parent <- getwd()
            dir <- file.path(parent, 'Chart')
            dir.create(dir, showWarnings = FALSE)
            setwd(dir)
            ##
            up_dot <- dotplot(kegg_up)
            path <- "KEGG_up_net"
            if (suffix != '') path = paste(suffix, path, sep="_")
            exportImage(up_dot, path, "png", width = img_width, height = img_height, resolution=150)
            ##
            down_dot <- dotplot(kegg_down)
            path <- "KEGG_down_net"
            if (suffix != '') path = paste(suffix, path, sep="_")
            exportImage(down_dot, path, "png", width = img_width, height = img_height, resolution=150)
            ##
            setwd(parent)
        }
    }
}

scoreHeatMap <- function(scores, 
                         export_img = T, 
                         img_width = 8,
                         img_height = 6,
                         suffix = '') {
    genes <- rep(rownames(scores), ncol(scores))
    conditions <- c()
    values <- c()
    for (c in c(1:ncol(scores))) {
        conditions <- c(conditions, rep(colnames(scores)[c], nrow(scores)))
        values <- c(values, scores[, c])
    }
    data <- data.frame(Condition = conditions, Genes = genes, Value = values)
    plt = ggplot(data, aes(Condition, Genes, fill = Value)) + 
            geom_tile()                    
    if (export_img) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart')
        dir.create(dir, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "hmap"
        if (suffix != '') path = paste(suffix, path, sep="_")
        exportImage(plt, path, "png", width = img_width, height = img_height, resolution=150)
        setwd(parent)
    }
    return(plt)
}

vennPlot <- function(data, 
                     color = c("white", "turquoise"), 
                     export_img = T, 
                     img_width = 8,
                     img_height = 6,
                     suffix = '') {
    plt = ggVennDiagram(data) + 
            scale_fill_gradient(low = color[1], high = color[2])
    #
    if (export_img) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart')
        dir.create(dir, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "venn"
        if (suffix != '') path = paste(suffix, path, sep="_")
        exportImage(plt, path, "png", width = img_width, height = img_height, resolution=150)
        setwd(parent)
    }
    return(plt)
}

maPlot <- function(result,
                   counts,
                   restype = 'edgeR',
                   fc_limit = c(1, Inf ),
                   p_threshold = 0.05, 
                   q_threshold = 0.05, 
                   export_img=T, 
                   img_width = 8,
                   img_height = 6,
                   suffix='') {
    
                             
    
}

volcanoPlot <- function(result,
                        restype = 'edgeR',
                        fc_limit = c(1, Inf ),
                        p_threshold = 0.05, 
                        q_threshold = 0.05, 
                        export_img=T, 
                        img_width = 8,
                        img_height = 6,
                        suffix='') {
    #
    plt <- NA
    #
    if (restype == 'edgeR') {}
    else if (restype == 'DESeq2') {
        # Set color flag
        upfilter <- !is.na(result$padj) & 
                result$pvalue < p_threshold & 
                result$padj < q_threshold & 
                result$log2FoldChange >= fc_limit[1] &
                result$log2FoldChange <= fc_limit[2]
        downfilter <- !is.na(result$padj) & 
                result$pvalue < p_threshold & 
                result$padj < q_threshold & 
                result$log2FoldChange <= -fc_limit[1] &
                result$log2FoldChange >= -fc_limit[2]
        result$DEG <- ifelse(upfilter, "Up", ifelse(downfilter, "Down", "No"))
        # Make plot
        plt = ggplot(data=result, aes(x=log2FoldChange, y=-log10(padj), col=DEG)) +
              theme_minimal() +
              geom_point() +
              scale_color_manual(values=c("blue", "black", "red"))
    }
    #
    if (export_img) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart')
        dir.create(dir, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "volcano"
        if (suffix != '') path = paste(suffix, path, sep="_")
        exportImage(plt, path, "png", width = img_width, height = img_height, resolution=150)
        setwd(parent)
    }
    return(plt)    
}

pathwayView <- function(genes, target, species, suffix = '') {
    pathview(gene.data = genes, 
               pathway.id = target,
               species = species,
               out.suffix = suffix)                           
}
