# CummeRbund




# DESeq2
deseqAnalyze <- function(data, group, formula, method = 'auto', reduced = as.formula('~1')) {
    dds <- DESeqDataSetFromMatrix(countData = data, colData = group, design = formula)
    if (method == 'auto') {
        dds <- DESeq(dds)
    }
    else if (method == 'lrt') {
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomLRT(dds, full = formula, reduced = reduced)
    }
    return(dds)                 
}

deseqResult <- function(data, target, export_csv=T, export_obj=T, suffix='') {
    # Get result
    res <- results(data, contrast = target)
    # Make dir.
    parent = getwd()
    dir.create(file.path(parent, 'Result'), showWarnings = FALSE)
    setwd(file.path(parent, 'Result'))
    # Save result object
    if (export_obj) saveRDS(res, file = paste(suffix, "_DEG.obj", sep=""))
    # Export data as CSV
    if (export_csv) write.csv(res, file = paste(suffix, "_DEG.csv", sep=""))
    setwd(parent)
    return(res)
}

# EdgeR



# Common
countData <- function(data, restype = 'edgeR', normalized=T, export_csv=T, suffix='') {
    # Get count
    counts <- NA
    if (restype == 'edgeR') {}
    else if (restype == 'DESeq2') {
        counts <- counts(data, normalized=normalized)
    }
    # Export data as CSV
    if (export_csv) {
        # Make dir.
        parent = getwd()
        dir.create(file.path(parent, 'Result', 'Count'), showWarnings = FALSE)
        setwd(file.path(parent, 'Result', 'Count'))
        write.csv(counts, file = paste(suffix, "_counts.csv", sep=""))
        setwd(parent)
    }
    return(counts)
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
        parent = getwd()
        dir.create(file.path(parent, 'Result'), showWarnings = FALSE)
        setwd(file.path(parent, 'Result'))
        ## 
        write.csv(up, file = paste(suffix, "_DEG_up.csv", sep=""))
        write.csv(down, file = paste(suffix, "_DEG_down.csv", sep=""))
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
                         dot_plot=T,
                         ont_network=T,
                         suffix='') {
    
    # Preset gene list
    upgenes <- rownames(data$up)
    downgenes <- rownames(data$down)
    ## Gene ID conversion
    if (genetype != "ENTREZID") {
        upgenes <- bitr(upgenes, fromType=genetype, toType="ENTREZID", OrgDb=org)
        upgenes <- upgenes[, 2]
        downgenes <- bitr(downgenes, fromType=genetype, toType="ENTREZID", OrgDb=org)
        downgenes <- downgenes[, 2]
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
            parent = getwd()
            dir.create(file.path(parent, 'Result', 'GO'), showWarnings = FALSE)
            setwd(file.path(parent, 'Result', 'GO'))
            ## 
            write.csv(go_up, file=paste(suffix, "_GO_", t, "_up.csv", sep=""))
            write.csv(go_down, file = paste(suffix, "_GO_", t, "_down.csv", sep=""))
            setwd(parent)
        }
        #
        if (export_img) {
            # Make dir.
            parent = getwd()
            dir.create(file.path(parent, 'Result', 'Chart'), showWarnings = FALSE)
            setwd(file.path(parent, 'Result', 'Chart'))
            ##
            if (dot_plot) {
                up_dot <- dotplot(go_up)
                down_dot <- dotplot(go_down)
                exportImage(up_dot, paste(suffix, "_GO_", t, "_up_dot", sep=""), "png", width = 3, height = 3, resolution=150)
                exportImage(down_dot, paste(suffix, "_GO_", t, "_down_dot", sep=""), "png", width = 3, height = 3, resolution=150)
            }
            ##
            if (ont_network) {
                up_net <- goplot(go_up)
                down_net <- goplot(go_down)
                exportImage(up_net, paste(suffix, "_GO_", t, "_up_net", sep=""), "png", width = 3, height = 3, resolution=150)
                exportImage(down_net, paste(suffix, "_GO_", t, "_down_net", sep=""), "png", width = 3, height = 3, resolution=150)               
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
                              suffix='') {
                             
    # Preset gene list
    upgenes <- rownames(data$up)
    downgenes <- rownames(data$down)
    ## Gene ID conversion
    if (genetype != "ENTREZID") {
        upgenes <- bitr(upgenes, fromType=genetype, toType="ENTREZID", OrgDb=org)
        print(upgenes[is.na(upgenes[,2]), 1])
        upgenes <- upgenes[, 2]
        downgenes <- bitr(downgenes, fromType=genetype, toType="ENTREZID", OrgDb=org)
        downgenes <- downgenes[, 2]
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
            parent = getwd()
            dir.create(file.path(parent, 'Result', 'KEGG'), showWarnings = FALSE)
            setwd(file.path(parent, 'Result', 'KEGG'))
            ## 
            write.csv(kegg_up, file=paste(suffix, "_KEGG_up.csv", sep=""))
            write.csv(kegg_down, file = paste(suffix, "_KEGG_down.csv", sep=""))
            setwd(parent)
        }
        #
        if (export_img) {
            # Make dir.
            parent = getwd()
            dir.create(file.path(parent, 'Result', 'Chart'), showWarnings = FALSE)
            setwd(file.path(parent, 'Result', 'Chart'))
            ##
            up_dot <- dotplot(kegg_up)
            down_dot <- dotplot(kegg_down)
            exportImage(up_dot, paste(suffix, "_KEGG_up_dot", sep=""), "png", width = 6, height = 6, resolution=150)
            exportImage(down_dot, paste(suffix, "_KEGG_down_dot", sep=""), "png", width = 6, height = 6, resolution=150)
            ##
            setwd(parent)
        }
    }
}

heatMap <- function(counts, export_img = T, suffix = '') {
    plt = ggplot(counts, aes(Condition, Genes, fill = Count)) + 
            geom_tile()                    
    if (export_img) {
        # Make dir.
        parent = getwd()
        dir.create(file.path(parent, 'Result', 'Chart'), showWarnings = FALSE)
        setwd(file.path(parent, 'Result', 'Chart'))
        ##
        exportImage(plt, paste(suffix, "_hmap", sep=""), "png", width = 8, height = 6, resolution=150)
        setwd(parent)
    }
    return(plt)
}

vennPlot <- function(data, target, color = c("white", "turquoise"), export_img = T, suffix = '') {
    plt = ggVennDiagram(x[1:3],label_alpha=0) + 
            scale_fill_gradient(low = color[1], high = color[2])
    #
    if (export_img) {
        # Make dir.
        parent = getwd()
        dir.create(file.path(parent, 'Result', 'Chart'), showWarnings = FALSE)
        setwd(file.path(parent, 'Result', 'Chart'))
        ##
        exportImage(plt, paste(suffix, "_venn", sep=""), "png", width = 8, height = 6, resolution=150)
        setwd(parent)
    }
    return(plt)
}

maPlot <- function(counts) {
    
                             
    
}

volcanoPlot <- function(result,
                       restype = 'edgeR',
                       fc_limit = c(1, Inf ),
                       p_threshold = 0.05, 
                       q_threshold = 0.05, 
                       export_img=T, suffix='') {
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
        parent = getwd()
        dir.create(file.path(parent, 'Result', 'Chart'), showWarnings = FALSE)
        setwd(file.path(parent, 'Result', 'Chart'))
        ##
        exportImage(plt, paste(suffix, "_volcano", sep=""), "png", width = 6, height = 6, resolution=150)
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
