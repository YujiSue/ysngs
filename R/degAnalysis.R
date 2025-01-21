# DESeq2

deseqData <- function(data, group, formula) {
    dds <- DESeqDataSetFromMatrix(countData = data, colData = group, design = formula)
    return(dds)
}

deseqAnalyze <- function(data, group, formula, method = 'auto', reduced = as.formula('~1'), level = NA, ref = NA) {
    dds <- deseqData(data, group, formula)
    if (!is.na(level)) dds[[level]] <- relevel(dds[[level]], ref)
    if (method == 'auto') {
        dds <- DESeq(dds)
    }
    else if (method == 'lrt') {
        dds <- DESeq(dds, test = "LRT", reduced = reduced)
    }
    return(dds)
}

deseqReAnalyze <- function(dds, formula = c(), method = 'auto', reduced = as.formula('~1'), level = NA, ref = NA) {
    if (!is.na(level)) dds[[level]] <- relevel(dds[[level]], ref)
    if (0 < length(formula)) design(dds) <- formula
    if (method == 'auto') {
        dds <- DESeq(dds)
    }
    else if (method == 'lrt') {
        dds <- DESeq(dds, test = "LRT", reduced = reduced)
    }
    return(dds)
}

deseqCount <- function(data, export_csv=T, suffix='') {
    ncounts <- counts(data, normalized = TRUE)
    # Export data as CSV
    if (export_csv) {
        # Make dir.
        path = "counts.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(ncounts, file = path)
    }
    return(ncounts)
}

deseqResult <- function(data, target, export_csv=T, suffix='') {
    res <- NA
    if (is.list(target)) res <- results(data, target)
    else res <- results(data, contrast = target)
    # Export data as CSV
    if (export_csv) {
        path = "DEG_result.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(res, file = path)
    }
    return(res)
}

# EdgeR
edgerAnalyze <- function(data, group, formula) {
    dge <- DGEList(counts = data, group = group)
    dge <- calcNormFactors(dge)
    
    dge <- estimateGLMCommonDisp(dge, formula)
    dge <- estimateGLMTrendedDisp(dge, formula)
    dge <- estimateGLMTagwiseDisp(dge, formula)


    fit <- glmFit(dge, formula)
    
    
    return(list(result = dge, test = lrt))
}

edgerCount <- function(data, export_csv=T, suffix='') {
    ncounts <- cpm(data$result, normalized.lib.sizes = TRUE)
    # Export data as CSV
    if (export_csv) {
        # Make dir.
        path = "counts.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(ncounts, file = path)
    }
    return(ncounts)
}

egderResult <- function(data, target, export_csv=T, suffix='') {
    lrt <- glmLRT(data, coef = target)
    res <- as.data.frame(topTags(lrt, n = nrow(data$result)))
    # Export data as CSV
    if (export_csv) {
        path = "DEG_result.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(res, file = path)
    }
    return(res)
}

# Common
geneConversion <- function(genes, db, src, dest, select = 'first') {
    conv <- mapIds(db, keys = genes, keytype = src, column = dest, multiVals = select)
    conv[is.na(conv)] <- genes[is.na(conv)]
    return(conv)
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
    if (restype == 'edgeR') {
        ## Up-regulated
        upfilter <- !is.na(result$FDR) &
                result$PValue < p_threshold &
                result$FDR < q_threshold &
                result$logFC >= fc_limit[1] &
                result$logFC <= fc_limit[2]
        up <- result[upfilter,]
        up <- up[order(up$logFC, decreasing = T),]

        ## Down-regulated
        downfilter <- !is.na(result$FDR) &
                result$PValue < p_threshold &
                result$FDR < q_threshold &
                result$logFC <= -fc_limit[1] &
                result$logFC >= -fc_limit[2]
        down <- result[downfilter,]
        down <- down[order(down$logFC, decreasing = F),]
    }
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
        down <- result[downfilter,]
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

zScoreData <- function(counts, group, factors,
                        collabel = c(),
                        rowlabel = c(),
                        cluster = T,
                        export_csv=F,
                        export_html=T,
                        image_width=640,
                        image_height=NA,
                        suffix='') {
    means <- list()
    for (f in factors) {
        means[[f]] <- rowMeans(counts[, group == f])
    }
    zscores <- do.call(cbind, means)
    m <- apply(zscores, 1, mean)
    s <- apply(zscores, 1, sd)
    zscores <- (zscores - m) / s
    #
    if (length(collabel) == 0) colnames(zscores) <- factors
    else colnames(zscores) <- collabel
    if (length(rowlabel) == 0) rownames(zscores) <- rownames(counts)
    else rownames(zscores) <- rowlabel
    #
    if (cluster) {
        d <- dist(zscores)
        hc <- hclust(d)
        zscores <- zscores[hc$order,]
    }
    if (export_csv) {
        path = "zscores.csv"
        if (suffix != '') path = paste(suffix, path, sep="_")
        write.csv(ncounts, file = path)
    }
    if (export_html) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart', 'HTML')
        dir.create(dir, recursive = T, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "heatmap.html"
        if (suffix != '') path = paste(suffix, path, sep="_")
        hmap <- plot_ly(x = colnames(zscores),
                        y = rownames(zscores),
                        z = as.matrix(zscores), colors = "OrRd",
                        type = "heatmap") %>%
                layout(xaxis = list(side ="top" ) )

        if (is.na(image_height)) image_height <- nrow(zscores) * 2
        hmap <- hmap %>% layout(autosize = F, width = image_width, height = image_height) %>% hide_colorbar()
        style <- paste("function(el,x){el.style.height='", image_width, "px';el.style.width='", image_height, "px';el.style.overflow='auto';}", sep="")
        exportHTML(hmap, path, style)
        setwd(parent)
    }
    return(zscores)
}

volcanoPlot <- function(result,
                        org,
                        genetype = 'ENTREZID',
                        restype = 'edgeR',
                        fc_limit = c(1, Inf ),
                        p_threshold = 0.05,
                        q_threshold = 0.05,
                        export_img=T,
                        img_width = 8,
                        img_height = 6,
                        export_html = T,
                        suffix='') {
    #
    plt <- NA
    plt2 <- NA
    #
    if (restype == 'edgeR') {
    }
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

        if (export_html) {
            plt2 <- plot_ly(type = 'scatter') %>%
            add_trace(
                name = "Upregulated",
                x = result[upfilter,"log2FoldChange"],
                y = -log10(result[upfilter,"padj"]),
                text = paste("Gene :", geneConversion(rownames(result[upfilter,]), db = org, src = genetype, dest = "SYMBOL"),
                            "<br>Fold change :", 2^result[upfilter, "log2FoldChange"],
                            "<br>Adjust. p-value :", result[upfilter, "padj"]),
                hoverinfo = 'text',
                xaxis = list(title = 'log2FoldChange'),
                yaxis = list(title = '-log10(padj)'),
                marker = list(color = 'red')
            ) %>%
            add_trace(
                name = "Downregulated",
                x = result[downfilter,"log2FoldChange"],
                y = -log10(result[downfilter,"padj"]),
                text = paste("Gene :", geneConversion(rownames(result[downfilter,]), db = org, src = genetype, dest = "SYMBOL"),
                            "<br>Fold change :", 2^result[downfilter, "log2FoldChange"],
                            "<br>Adjust. p-value :", result[downfilter, "padj"]),
                hoverinfo = 'text',
                xaxis = list(title = 'log2FoldChange'),
                yaxis = list(title = '-log10(padj)'),
                marker = list(color = 'blue')
            ) %>%
            add_trace(
                name = "Not regulated",
                x = result[!(upfilter | downfilter),"log2FoldChange"],
                y = -log10(result[!(upfilter | downfilter),"padj"]),
                text = paste("Gene :", geneConversion(rownames(result[!(upfilter | downfilter),]), db = org, src = genetype, dest = "SYMBOL"),
                            "<br>Fold change :", 2^result[!(upfilter | downfilter), "log2FoldChange"],
                            "<br>Adjust. p-value :", result[!(upfilter | downfilter), "padj"]),
                hoverinfo = 'text',
                xaxis = list(title = 'log2FoldChange'),
                yaxis = list(title = '-log10(padj)'),
                marker = list(color = 'black')
            )
        }
    }
    #
    if (export_img) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart')
        dir.create(dir, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "Volcano"
        if (suffix != '') path = paste(suffix, path, sep="_")
        exportImage(plt, path, "png", width = img_width, height = img_height, resolution=150)
        setwd(parent)
    }
    if (export_html) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart', 'HTML')
        dir.create(dir, recursive = T, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "Volcano.html"
        if (suffix != '') path = paste(suffix, path, sep="_")
        exportHTML(plt2, path)
        setwd(parent)
    }
    return(plt)
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
        upgenes <- geneConversion(genes = upgenes, db = org, src = genetype, dest = "ENTREZID")
        downgenes <- geneConversion(genes = downgenes, db = org, src = genetype, dest = "ENTREZID")
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
                if (0 < nrow(go_up)) {
                    up_dot <- dotplot(go_up)
                    path <- paste("GO", t, "up", "dot", sep="_")
                    if (suffix != '') path = paste(suffix, path, sep="_")
                    exportImage(up_dot, path, "png", width = img_width, height = img_height, resolution=150)
                }
                ##
                if (0 < nrow(go_down)) {
                    down_dot <- dotplot(go_down)
                    path <- paste("GO", t, "down", "dot", sep="_")
                    if (suffix != '') path = paste(suffix, path, sep="_")
                    exportImage(down_dot, path, "png", width = img_width, height = img_height, resolution=150)
                }
            }
            ##
            if (ont_network) {
                ##
                if (0 < nrow(go_up)) {
                    up_net <- goplot(go_up)
                    path <- paste("GO", t, "up", "net", sep="_")
                    if (suffix != '') path = paste(suffix, path, sep="_")
                    exportImage(up_net, path, "png", width = img_width, height = img_height, resolution=150)
                }
                ##
                if (0 < nrow(go_down)) {
                    down_net <- goplot(go_down)
                    path <- paste("GO", t, "down", "net", sep="_")
                    if (suffix != '') path = paste(suffix, path, sep="_")
                    exportImage(down_net, path, "png", width = img_width, height = img_height, resolution=150)
                }
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
        upgenes <- geneConversion(genes = upgenes, db = org, src = genetype, dest = "ENTREZID")
        downgenes <- geneConversion(genes = downgenes, db = org, src = genetype, dest = "ENTREZID")
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
            if (0 < nrow(kegg_up)) {
                up_dot <- dotplot(kegg_up)
                path <- "KEGG_up_dot"
                if (suffix != '') path = paste(suffix, path, sep="_")
                exportImage(up_dot, path, "png", width = img_width, height = img_height, resolution=150)
            }
            ##
            if (0 < nrow(kegg_down)) {
                down_dot <- dotplot(kegg_down)
                path <- "KEGG_down_dot"
                if (suffix != '') path = paste(suffix, path, sep="_")
                exportImage(down_dot, path, "png", width = img_width, height = img_height, resolution=150)
            }
            ##
            setwd(parent)
        }
    }
}

vennPlot <- function(data,
                     color = c("white", "turquoise"),
                     export_img = T,
                     img_width = 8,
                     img_height = 6,
                     suffix = '') {
    plt <- ggVennDiagram(data) +
            scale_fill_gradient(low = color[1], high = color[2]) +
            scale_x_continuous(expand = expansion(mult = .2))
    #
    if (export_img) {
        # Make dir.
        parent <- getwd()
        dir <- file.path(parent, 'Chart')
        dir.create(dir, showWarnings = FALSE)
        setwd(dir)
        ##
        path = "Venn"
        if (suffix != '') path = paste(suffix, path, sep="_")
        exportImage(plt, path, "png", width = img_width, height = img_height, resolution=150)
        setwd(parent)
    }
    return(plt)
}

pathwayView <- function(genes, target, species, suffix = '') {
    pview <- pathview(gene.data = genes,
                      pathway.id = target,
                      species = species,
                      out.suffix = suffix)
    return(pview)
}
