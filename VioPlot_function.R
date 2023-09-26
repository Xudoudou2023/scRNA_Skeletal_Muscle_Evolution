library(reshape2)
VioPlot <- function (object, features, ident = Idents(object), assay = "RNA", 
                     facet.toward = "row", do.log = F, strip.text.size = 8, strip.text.angle = NULL, 
                     strip.text.hjust = NULL, strip.text.vjust = NULL, color.use = NULL, 
                     jitter = F, cell.order = NULL, add.ave.point = T, add.line = T, 
                     line.size = 0.5, add.box = F, fill.expression = F, fill.col = NA, 
                     mid.color.point = "auto", legend.position = "top", 
                     legend.name = "", jitter.size = 0.1, text.x.size = 8, text.y.size = 8, 
                     hjust = 0, vjust = 0, legend.scale = NULL, legend.title.size = 8, 
                     axis.text.x.angle = NULL, ...) 
{
  data <- as.data.frame(t(as.matrix(object@assays[[assay]]@data[features, 
                                                                , drop = FALSE])))
  p <- vioplot(data, features, ident = ident, facet.toward = facet.toward, 
               do.log = do.log, strip.text.size = strip.text.size, strip.text.angle = strip.text.angle, 
               strip.text.hjust = strip.text.hjust, strip.text.vjust = strip.text.vjust, 
               color.use = color.use, jitter = jitter, cell.order = cell.order, 
               add.ave.point = add.ave.point, add.line = add.line, line.size = line.size, 
               add.box = add.box, fill.expression = fill.expression, 
               fill.col = fill.col, mid.color.point = mid.color.point, 
               legend.position = legend.position, legend.name = legend.name, 
               jitter.size = jitter.size, text.x.size = text.x.size, 
               text.y.size = text.y.size, hjust = hjust, vjust = vjust, 
               legend.scale = legend.scale, legend.title.size = legend.title.size, 
               axis.text.x.angle = axis.text.x.angle, ...)
  return(p)
}

vioplot <- function (data, paths, ident, facet.toward = "row", do.log = F, 
                     strip.text.size = 8, strip.text.angle = NULL, strip.text.hjust = NULL, 
                     strip.text.vjust = NULL, color.use = NULL, jitter = F, cell.order = NULL, 
                     add.ave.point = T, add.line = T, line.size = 0.5, add.box = F, fill.expression = F,
                     fill.col = NA, mid.color.point = "auto", legend.position = "top", 
                     legend.name = "", jitter.size = 0.1, text.x.size = 8, text.y.size = 8, 
                     hjust = 0, vjust = 0, legend.scale = NULL, legend.title.size = 8, 
                     axis.text.x.angle = NULL) 
{
  sig_path <- as.data.frame(data[, paths, drop = FALSE])
  if (do.log) {
    sig_path <- log2(sig_path + 1)
  }
  sig_path$groups <- ident
  matm <- melt(sig_path, id.vars = c("groups"), measure.vars = paths, 
               variable.name = "pathway", value.name = "TPM")
  colnames(matm) <- c("groups", "pathway", "TPM")
  ave <- tapply(matm$TPM, list(matm$groups, matm$pathway), 
                mean)
  ave <- as.data.frame(ave)
  ave$groups <- rownames(ave)
  ave <- reshape2::melt(ave, id.vars = c("groups"), measure.vars = paths, 
                        variable.name = "pathway", value.name = "ave", factorsAsStrings = T)
  matm <- suppressMessages(inner_join(matm, ave))
  if (!is.null(cell.order)) {
    matm$groups <- factor(matm$groups, levels = cell.order)
  }
  color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", 
                    "#3C5488FF", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", 
                    "#848482", "#008856", "#E68FAC", "#0067A5", "#604E97", 
                    "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", 
                    "#654522", "#E25822", "#2B3D26", "#848482", "#008856", 
                    "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", 
                    "#DCD300", "#882D17")
  if (fill.expression) {
    if (mid.color.point == "auto") {
      mid.color <- (max(matm$ave) - min(matm$ave))/2
    }
    else {
      mid.color <- mid.color.point
    }
    p <- ggplot(matm, aes(x = groups, y = TPM, fill = ave))
    p <- p + geom_violin(trim = T, scale = "width") + theme_bw()
    if (!is.null(legend.scale)) {
      p <- p + scale_fill_gradientn(colours = fill.col,
                                    values = scales::rescale(legend.scale), guide = "colorbar",
                                    limits = legend.scale, name = legend.name)
    }
    else {
      p <- p + scale_fill_gradient2(name = legend.name,
                                    low = fill.col[1], mid = fill.col[2], high = fill.col[3],
                                    midpoint = mid.color)
    }
  }
  else {
    p <- ggplot(matm, aes(x = groups, y = TPM, color = groups,fill = groups)) + scale_fill_manual(values = fill.col)
    p <- p + geom_violin(trim = TRUE,width = 0.6,scale = 'width') + theme_bw()
    if (is.null(color.use) & length(unique(matm$groups)) <
        31) {
      p <- p + scale_color_manual(values = color_palette)
    }
    else {
      p <- p + scale_color_manual(values = color.use)
    }
    legend.position <- "none"
  }
  if (jitter) {
    p <- p + geom_jitter(size = jitter.size)
  }
  if (add.ave.point) {
    p <- p + geom_point(data = matm, aes(x = groups, y = ave, 
                                         group = 1), color = "Black", shape = 15, size = 1.5)
  }
  if (add.line) {
    p <- p + geom_line(data = matm, aes(x = groups, y = ave, 
                                        group = 1), color = "Black", size = line.size, linetype = "dashed")
  }
  if (add.box) {
    p <- p + geom_boxplot(width = 0.1, fill = "white", color = "Black", 
                          outlier.shape = NA)
  }
  if (facet.toward == "row") {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
                                90, axis.text.x.angle)
    p <- p + facet_grid(rows = vars(pathway), scales = "free") + 
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA), 
            panel.grid.minor = element_line(colour = NA)) + 
      theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                                       size = text.x.size, hjust = hjust, vjust = vjust), 
            axis.text.y = element_text(size = text.y.size))
  }
  else if (facet.toward == "col") {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
                                270, axis.text.x.angle)
    p <- p + coord_flip()
    p <- p + facet_grid(cols = vars(pathway), scales = "free") + 
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA), 
            panel.grid.minor = element_line(colour = NA)) + 
      theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                                       size = text.x.size, hjust = hjust, vjust = vjust), 
            axis.text.y = element_text(size = text.y.size))
  }
  p <- p + xlab("") + ylab("")
  p <- p + theme(legend.title = element_text(size = legend.title.size), 
                 strip.text = element_text(size = strip.text.size, angle = strip.text.angle, 
                                           hjust = strip.text.hjust, vjust = strip.text.vjust), 
                 strip.background = element_blank(), strip.placement = "outside")
  p <- p + guides(color = FALSE)
  return(p)
}
