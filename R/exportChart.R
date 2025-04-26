#
exportImage <- function(plot,output,format="png",width=3,height=3,resolution=150) {
    ggsave(paste(output, format, sep = "."), plot, width=width,height=height,dpi=resolution)
}
#
exportImageColab <- function(plot,output,format="png",width=3,height=3,resolution=150) {
    if (format == "svg" || format == "wmf") {
        temp = paste(output, "eps", sep = ".")
        ggsave(temp, plot, width=width, height=height, dpi=resolution)
        opt = ""
        if (format == "svg") opt = paste("--export-plain-svg=", output, ".svg", sep = "")
        else if (format == "wmf") opt = paste("--export-wmf=", output, ".wmf", sep = "")
        system(paste("dbus-run-session",
              "inkscape",
              opt,
              temp,
              sep=" "))
        system(paste("rm", temp, sep = " "))
    }
    else exportImage(plot,output,format=format,width=width,height=height,resolution=resolution)
}
#
exportPlotly2JSON <- function(plot, output) {
    plotjs <- plotly_json(plot, jsonedit = FALSE)
    write(plotjs, output)
}
#
exportPlotly2HTML <- function(plot, output, style = NA) {
    if (is.na(style)) widget <- plot
    else widget <- htmlwidgets::onRender(plot, style)
    htmlwidgets::saveWidget(
        widget,
        file = output,
        selfcontained = TRUE
    )
}
#
exportHTML <- function(plot, output, style = NA) {
    if (is.na(style)) widget <- plot
    else widget <- htmlwidgets::onRender(plot, style)
    htmlwidgets::saveWidget(
        widget,
        file = output,
        selfcontained = TRUE
    )
}
