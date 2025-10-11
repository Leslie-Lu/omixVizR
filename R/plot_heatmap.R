#' @title plot_heatmap
#' @description Generate a heatmap for a given correlation matrix.
#' @author Zhen Lu <luzh29@mail2.sysu.edu.cn>
#' @param cor_mat A numeric matrix representing the correlation values.
#' @param scale Character string indicating how to scale the data. Default is 'none'.
#' @param show_rownames Logical indicating whether to show row names. Default is FALSE.
#' @param show_colnames Logical indicating whether to show column names. Default is FALSE.
#' @param labels_row Character vector of labels for the rows. Default is NULL.
#' @param labels_col Character vector of labels for the columns. Default is NULL.
#' @param breaks Numeric vector of breakpoints for the heatmap colors. Default is seq(-1, 1, length.out = 101).
#' @param legend_breaks Numeric vector of breakpoints for the legend. Default is c(-1, -0.5, 0, 0.5, 1).
#' @param legend_labels Character vector of labels for the legend. Default is c("-1.0", "-0.5", "0", "0.5", "1.0").
#' @param legend_name Character string for the legend title. Default is 'Correlation'.
#' @param output_path Character string for the output file path. Default is NULL.
#' @param display Logical indicating whether to display the plot. Default is TRUE.
#' @param render_as Character string indicating the file format to save the plot. Default is 'png'.
#' @param width Numeric indicating the width of the plot. Default is 18.
#' @param height Numeric indicating the height of the plot. Default is 16.
#' @param dpi Numeric indicating the resolution of the plot. Default is 600.
#' @return A heatmap plot.
#' @details This function generates a heatmap for a given correlation matrix.
#' @rdname plot_heatmap
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom gtable gtable_filter
#' @importFrom grid textGrob gpar unit gTree gList viewport
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
plot_heatmap = function(cor_mat,
                        scale = "none",
                        show_rownames = FALSE,
                        show_colnames = FALSE,
                        labels_row = NULL,
                        labels_col = NULL,
                        breaks = seq(-1, 1, length.out = 101),
                        legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                        legend_labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
                        legend_name = "Correlation",
                        output_path = NULL,
                        display = TRUE,
                        render_as = "png",
                        width = 18,
                        height = 16,
                        dpi = 600) {
  if(show_rownames & is.null(labels_row)) stop("Please provide labels_row when show_rownames is TRUE!")
  if(show_colnames & is.null(labels_col)) stop("Please provide labels_col when show_colnames is TRUE!")
  cor_plot = pheatmap::pheatmap(cor_mat,
                                scale = scale,
                                show_rownames = show_rownames,
                                show_colnames = show_colnames,
                                labels_row = labels_row,
                                labels_col = labels_col,
                                breaks = breaks,
                                legend_breaks = legend_breaks,
                                legend_labels = legend_labels
                                )
  gt <- cor_plot$gtable
  gt <- gtable::gtable_filter(gt, "row_tree", invert = TRUE)
  gt <- gtable::gtable_filter(gt, "dend", invert = TRUE)
  legend_pos <- which(grepl("legend|guide", gt$layout$name))
  if(length(legend_pos) > 0) {
    legend_layout <- gt$layout[legend_pos[1], ]
    legend_grob <- gt$grobs[[legend_pos[1]]]
    
    title_grob <- grid::textGrob(legend_name,
                                gp = grid::gpar(fontsize = 12, fontface = "bold"),
                                just = "left",
                                x = grid::unit(-0.02, "npc"))
    
    new_legend_grob <- gridExtra::arrangeGrob(
      title_grob,
      legend_grob,
      heights = grid::unit(c(2.8, 4), c("lines", "null")),
      ncol = 1
    )
    
    gt$grobs[[legend_pos[1]]] <- new_legend_grob
  }
  final_plot <- grid::gTree(
    children = grid::gList(gt),
    vp = grid::viewport(x = 0.49, y = 0.5, width = 1, height = 1)
  )
  if (!(render_as == "rmarkdown")) {
    if(is.null(output_path)){
        output_path <- tempdir()
      }
    filename = file.path(output_path, paste0("heatmap_plot.", render_as))
    message(paste0("Saving image to ", filename))
    ggplot2::ggsave(
      filename = filename,
      plot = final_plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = render_as
    )
    if (display == TRUE) {
      invisible(system(paste0('open "', filename, '"')))
    }
  }else{
    final_plot
  }
}
