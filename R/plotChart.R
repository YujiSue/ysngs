



# Simple heatmap
clusteredMap <- function(data, #matrix
                        color = heat.colors(256)) {
    dist <- listdist(data)
    cluster <- hclust(dist)
    heatmap(data, Rowv=as.dendrogram(cluster), col = color)
}
# Heatmap

#' Make heatmap object
#'
#' @param data Dataframe
#' @param labels List of axis labels [x,y,z]
#' @param color List of colors (2 or 3 colors), or colorspace name (for only plot_ly)
#' @param zrange Limit and midpoint value
#' @param showscale boolean to show color scale
#' @param lib 'ggplot2' or 'plotly'
#' @return ggplot2 or plot_ly object
#' @export
hmapPlot <- function(data, #data.frame
                    labels,
                    color,
                    zrange = list(min=NA,max=NA,mid=NA),
                    showscale = TRUE,
                    lib = "ggplot") {
    #
    range <- zrange
    if (is.na(range[["max"]])) range[["max"]] <- max(unlist(data), na.rm = TRUE)
    if (is.na(range[["min"]])) range[["min"]] <- min(unlist(data), na.rm = TRUE)
    if (is.na(range[["mid"]])) range[["mid"]] <- (range[["max"]] + range[["min"]]) / 2
    #
    if (lib == "ggplot") {
        df <- melt(as.matrix(data))
        colnames(df) <- c(labels[["y"]], labels[["x"]], labels[["z"]])
        plt <- ggplot(df, aes(x = !!sym(labels[["x"]]), y = !!sym(labels[["y"]]), fill = !!sym(labels[["z"]]))) +
                geom_tile()
        #
        if (is.vector(color)) {
            if (length(color) == 2) {
                plt <- plt + scale_fill_gradient(low = color[1], high = color[2], limits = c(range[["min"]], range[["max"]]))
            }
            else if (length(color) == 3) {
                plt <- plt + scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = range[["mid"]], limits = c(range[["min"]], range[["max"]]))
            }
            else stop("Number of colors should be 2 or 3.")
        }
        else stop("Set list of colors for displaying heatmap using ggplot2.")
        plt <- plt + theme_minimal()
        return(plt)
    }
    else if (lib == "plotly") {
        cs <- NA        
        if (is.character(color) && length(color) == 1) cs = color
        else if (is.vector(color)) {
            if (length(color) == 2) {
                cs = list(
                    list(0, color[1]),
                    list(1, color[2])
                )
            }
            else if (length(color) == 3) {
                cs = list(
                    list(0, color[1]),
                    list(0.5, color[2]),
                    list(1, color[3])
                )
            }
            else stop("Number of colors should be 2 or 3.")
        }
        else stop("Set colorspace name or list of colors for displaying heatmap using plotly.")
        plt <- plot_ly(x = colnames(data),
                        y = rownames(data),
                        z = as.matrix(data),
                        type = "heatmap",
                        colorscale = cs,
                        showscale = showscale)
        
        plt <- plt %>% layout(
            coloraxis = list(
                zmin = range[["min"]],
                zmax = range[["max"]],
                zmid = range[["mid"]]
            ),
            xaxis = list(side ="top", title = labels[["x"]]),
            yaxis = list(title = labels[["y"]], autorange = "reversed"),
            scene = list(zaxis = list( title = labels[["z"]]))
        )
        return(plt)
    }
    else stop("'lib' argment should be 'ggplot' or 'plotly' under the current version.")
}


# Venn chart
vennPlot <- function(data,
                     color = c("white", "turquoise")) {
    plt <- ggVennDiagram(data) +
            scale_fill_gradient(low = color[1], high = color[2]) +
            scale_x_continuous(expand = expansion(mult = .2))
    return(plt)
}