#' @title plot_forest
#' @description Creates a forest plot overlaid on a table.
#' @author Zhen Lu <luzh29@mail2.sysu.edu.cn>
#' @param p_left_data Data frame (required). The information to be displayed to the left of the forest plot.
#' @param point_estimate Vector. The point estimates to be displayed in the forest plot.
#' @param ci_lower_bound Vector. The lower confidence bounds.
#' @param ci_upper_bound Vector. The upper confidence bounds.
#' @param ci_sep String. What should separate the low and high confidence bounds? Default " to ".
#' @param p_right_data Data frame (optional). Information to be displayed on the right side of the table. If not supplied, an Estimate column is generated automatically.
#' @param precision_digits Integer. The number of decimal places on the point_estimate (default 3)
#' @param p_mid_width Integer. The width of forest plot in characters (default 30)
#' @param null_line_at Numeric. Default 1. Change to 0 if using absolute measures.
#' @param output_path String. Where to save the image, default tempdir().
#' @param dpi Numeric. The image resolution in dpi, default 600
#' @param display Logical. Should the file be opened? Default TRUE.
#' @param font_family String or character vector. The font to use for the ggplot and table. Default c("MetroSans", "mono"). The first available font is used.
#' @param p_left_data_name String or String Vector. The name vector for the left side data. Default NULL, which uses the column names of p_left_data.
#' @param p_right_data_name String. The name vector for the right side data. Default "Effect size (95% CI)"
#' @param stripe_colour Hex String. Colour to use for the table stripes, default "#eff3f2".
#' @param background_colour Hex String or Colour Name. The colour of the background, default "white".
#' @param x_scale_linear Logical. Default TRUE, change to FALSE for a log scale.
#' @param xlim Vector. Manually specify limits for the x axis as a vector length 2, i.e. c(low, high)
#' @param xbreaks Vector. X axis breaks to label. Specify limits in xlim if using this option.
#' @param nudge_x Numeric. Nudge the alignment horizontally. Default 1. Higher values make the entire plot wider and consequently space out the elements of the figure.
#' @param nudge_y Numeric. Allows small changes to the vertical alignment of the forest plot points. 1 unit is approximately the height of 1 row.
#' @param nudge_height Numeric. Adjust the overall height of the plot output. Default is 0.
#' @param nudge_width Numeric. Adjust the overall width of the plot output. Default is 0.
#' @param justify Numeric Vector. This should be a numeric vector either of length 1 (in which case it will apply to all columns) or of length equal to the number of columns in p_left_data + 1 (for the point_estimate column). Each number in the vector dictates the column justification, with 0 being left, 0.5 being center, and 1 being right.
#' @param arrows Logical. Should there be arrows displayed below the ggplot? Default FALSE. Specify xlim if using arrows.
#' @param arrow_labels String Vector, length 2. Labels for the arrows. Set arrows to TRUE or this will have no effect.
#' @param risk_colors Vector. Length 2. Colors for the effect measure. Default is ggsci::pal_npg("nrc")(2).
#' @param arrow_nudge_y Numeric. Nudge the vertical position of the arrows. Default 0.
#' @param add_plot A ggplot object to add to the right side of the table. To align correctly with rows, 1 unit is the height of a row and y = 0 for the center of the bottom row.
#' @param add_plot_width Numeric. Width to display add_plot. Relative to the width of the forest plot, where 1 (the default) is the same width.
#' @param add_plot_gap Logical. Should there be space added between the plot and the main figure? Default FALSE.
#' @param point_sizes Vector. Length should be equal to 1 or nrow(p_left_data). The sizes of the points in the center plot, where 3.25 is the default.
#' @param point_shapes Vector. Length should be equal to 1 or nrow(p_left_data). The shapes of the points in the center plot, where 16 (a filled circle) is the default.
#' @param p_mid_forest A ggplot object to use instead of the central plot.
#' @param lower_header_row Logical. If TRUE, drops the header down one row (In the table rather than above it, like the default value (FALSE))
#' @param render_as String or Function. What output format should be used? Option is passed to ggplot2::ggsave() as the argument "device". Either pass a device function (e.g. png) or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#' @param table_theme A gridExtra table theme. If specified, overwrites all table theme customization in other options. The default is a modified version of ttheme_minimal.
#'
#' @return image
#' @details This function creates a forest plot overlaid on a table. It is highly customizable, allowing for adjustments to fonts, colors, and layout.
#' @rdname plot_forest
#' @importFrom rlang .data
#' @importFrom graphics text
#' @export
#'
plot_forest <- function(p_left_data,
                        point_estimate,
                        ci_lower_bound,
                        ci_upper_bound,
                        ci_sep = ", ",
                        p_right_data = NULL,
                        precision_digits = 3,
                        p_mid_width = 30,
                        null_line_at = 1,
                        output_path = NULL,
                        dpi = 600,
                        display = TRUE,
                        font_family = c("MetroSans", "mono"),
                        p_left_data_name = NULL,
                        p_right_data_name = "Effect size (95% CI)",
                        stripe_colour = "#eff3f2",
                        background_colour = "white",
                        x_scale_linear = TRUE,
                        xlim = NULL,
                        xbreaks = NULL,
                        nudge_y = 0,
                        nudge_x = 1,
                        nudge_height = 0,
                        nudge_width = 0,
                        justify = 0,
                        arrows = FALSE,
                        arrow_labels = c("Lower", "Higher"),
                        risk_colors = ggsci::pal_npg("nrc")(2),
                        arrow_nudge_y = 0,
                        add_plot = NULL,
                        add_plot_width = 1, #TODO
                        add_plot_gap = FALSE,
                        point_sizes = 3,
                        point_shapes = 16,
                        p_mid_forest = NULL,
                        lower_header_row = FALSE,
                        render_as = "png",
                        table_theme = NULL){
  resolve_font_family <- function(candidates) {
    if (length(candidates) == 0) {
      candidates <- "MetroSans"
    }

    candidates <- as.character(candidates)
    candidates <- unique(candidates[!is.na(candidates)])
    candidates <- candidates[nzchar(candidates)]

    if (!length(candidates)) {
      candidates <- "MetroSans"
    }

    for (candidate in candidates) {
      available_fonts <- sysfonts::font_families()
      if (candidate %in% available_fonts) {
        return(candidate)
      }
    }

    return("MetroSans")
  }
  font_family <- resolve_font_family(font_family)
  font_family_group <- dplyr::case_when(
    identical(font_family, "mono") ~ "mono",
    identical(font_family, "serif") ~ "serif",
    grepl("sans", font_family, ignore.case = TRUE) ~ "sans",
    TRUE ~ "other"
  )
  font_dir <- system.file("extdata", package = "omixVizR")
  if (!"MetroSans" %in% sysfonts::font_families()) {
    sysfonts::font_add(
      family = "MetroSans",
      regular = file.path(font_dir, "MetroSans-Regular.ttf"),
      bold = file.path(font_dir, "MetroSans-Bold.ttf"),
      bolditalic = file.path(font_dir, "MetroSans-BoldItalic.ttf")
    )
  }
  showtext::showtext_auto(enable = TRUE)
  showtext::showtext_opts(dpi = 600)

  if(!length(justify) == 1){
    justify <- c(justify[1:(ncol(p_left_data))], 0, justify[(ncol(p_left_data)+1):length(justify)])
    justify <- matrix(justify, ncol=length(justify), nrow = nrow(p_left_data) + 3, byrow=TRUE)
    justify <- as.vector(justify)
  }

  if (!is.null(table_theme)) {
    theme <- table_theme
  } else if (lower_header_row == FALSE) {
    theme <- gridExtra::ttheme_minimal(core = list(
      fg_params = list(hjust = justify, x = (0.05 + (0.45/0.5) * justify), fontfamily = font_family),
      bg_params = list(fill = c(rep(c(stripe_colour, background_colour), length.out = nrow(p_left_data)), background_colour, background_colour, background_colour))
    ),
      colhead = list(fg_params = list(hjust = c(0,0, rep(0.5, ncol(p_left_data)+ncol(p_right_data))), x = c(0.05,0.05, rep(0.5, ncol(p_left_data)+ncol(p_right_data))),
                                      fontfamily = font_family,
                                      # parse = TRUE,
                                    fontface=c(rep("bold", ncol(p_left_data)+ncol(p_right_data)-1), "bold.italic", "bold")),
                   bg_params = list(fill = background_colour))
    )
  }else{
    theme <- gridExtra::ttheme_minimal(core = list(
      fg_params = list(hjust = justify, x = (0.05 + (0.45/0.5) * justify), fontfamily = font_family),
      bg_params = list(fill = c(rep(c(background_colour, stripe_colour), length.out=nrow(p_left_data)), background_colour, background_colour, background_colour))
    ),
    colhead = list(fg_params = list(hjust = 0, x = 0.05,
                                    fontfamily = font_family),
                   bg_params = list(fill = background_colour))
    )
  }

  gdata <- data.frame(point_estimate = point_estimate,
                      ci_lower_bound = ci_lower_bound,
                      ci_upper_bound = ci_upper_bound,
                      effect_color = ifelse(point_estimate >= null_line_at, risk_colors[1], risk_colors[2])
                      )
  if (lower_header_row){
    gdata <- tibble::add_row(gdata, .before = 1)
  }

  if (is.null(p_right_data)){
    tdata <- gdata

    tdata <- dplyr::mutate_all(tdata, ~sprintf(.,
        fmt = paste0('%#.', precision_digits,'f')
    ))

    tdata[tdata == "NA"] <- " "
    # pretty formatting for confidence intervals
    p_right_data <- data.frame(`Effect size (95% CI)` = ifelse(tdata$point_estimate == " ",
                                  " ", paste0(tdata$point_estimate, " (", tdata$ci_lower_bound,
                                      ci_sep, tdata$ci_upper_bound, ")")))
  }
  names(p_right_data) <- p_right_data_name

  # finds width in number of characters for monospaced font
  find_width_mono <- function(data){
    num_of_rows <- nrow(data)
    num_of_cols <- ncol(data)

    print_data <- dplyr::mutate_all(data, as.character)

    num_char_across <- 0
    width <- 0

    for (i in 1:num_of_cols) {
      for (j in 1:num_of_rows) {
        num_char_across[j] <- nchar(print_data[j, i])
      }
      width[i] <- max(max(num_char_across, na.rm = TRUE),
                      nchar(colnames(print_data)[i]), na.rm = TRUE)
    }
    return(sum(width, na.rm = TRUE))
  }
  # finds width using shape_string from the systemfonts package
  # if not using monospaced font
  find_width <- function(data){
    num_of_rows <- nrow(data)
    num_of_cols <- ncol(data)

    print_data <- dplyr::mutate_all(data, as.character)

    width <- 0

    names <- colnames(print_data)

    for (i in 1:num_of_cols) {
      temp <- systemfonts::shape_string(print_data[[names[i]]], family = font_family)
      temp_col <- systemfonts::shape_string(names[i], family = font_family)
      width[i] <- max(max(temp$metrics$width, na.rm = TRUE),
                      temp_col$metrics$width, na.rm = TRUE)
    }
    return(sum(width, na.rm = TRUE)/7.2)
  }
  # calculate widths for each side with the appropriate function
  if (font_family == "mono") {
    left_width <- find_width_mono(p_left_data)
    right_width <- find_width_mono(p_right_data)
  }else{
    left_width <- find_width(p_left_data)
    right_width <- find_width(p_right_data)
  }

  # insert a blank column so we can put the ggplot object on top
  # and correctly order columns
  total_width <- left_width + right_width + p_mid_width

  names(p_left_data) = p_left_data_name
  tdata_print <- p_left_data

  if (lower_header_row) {
    rbind.data.frame(colnames(tdata_print), tdata_print)
  }

  tdata_print$`  ` <- paste(rep(" ", times = round(p_mid_width, 0)),
                           collapse = '')
  tdata_print <- cbind(tdata_print, p_right_data)
  tdata_print <- tibble::add_row(tdata_print)
  tdata_print <- tibble::add_row(tdata_print)
  tdata_print <- tibble::add_row(tdata_print)
  tdata_print <- dplyr::mutate_all(tdata_print, as.character)
  tdata_print[is.na(tdata_print)] <- " "

  ## formatting functions
  mono_column <- function(table, col, font_family){
    col_indexes <- function(table, col, name="core-fg"){
      l <- table$layout
      which(l$l == col & l$name == name)
    }

    ind <- col_indexes(table, col, "core-fg")

    for (i in ind) {
      table$grobs[i][[1]][["gp"]] <- grid::gpar(fontfamily = font_family)
    }
    return(table)
  }
  white_column <- function(table, col, font_family){
    col_indexes <- function(table, col, name="core-bg"){
      l <- table$layout
      which(l$l == col & l$name == name)
    }

    ind <- col_indexes(table, col, "core-bg")
    ind_fg <- col_indexes(table, col, "core-fg")

    for (i in ind) {
      table$grobs[i][[1]][["gp"]] <- grid::gpar(fill = background_colour, col = background_colour)
    }

    for (i in ind_fg) {
      table$grobs[i][[1]][["gp"]] <- grid::gpar(fontfamily = font_family)
    }
    return(table)
  }

  ######## calculations for the top and bottom of the plot
  gdata$row_num <- (nrow(gdata) - 1):0
  h_adj <- dplyr::case_when(
    font_family_group == "mono" ~ 0.2,
    font_family_group == "serif" ~ .43,
    font_family_group == "sans" ~ .37,
    TRUE ~ 0
  )
  h_adj <- nudge_y + h_adj
  slope_adj <- dplyr::case_when(
    font_family_group == "mono" ~ -0.175,
    font_family_group == "serif" ~ -.19,
    font_family_group == "sans" ~ -.16,
    TRUE ~ 0
  )
  font_adj <- 0.3 + h_adj + log(nrow(gdata)) * slope_adj
  y_low <- -.5 + font_adj + -.1381 * log(nrow(gdata))
  y_high <- 1.017 * nrow(gdata) - 0.6

  #### add shapes and sizes to gdata ########
  gdata$shape <- point_shapes
  gdata$sizes <- point_sizes

  #### if a ci will be out of bounds, add arrow on the oob side  ###############
  g_oob <- tibble::tibble()
  if (!is.null(xlim)) {
    oob_arrows <- gdata

    oob_arrows$x_low <- xlim[1]
    oob_arrows$x_high <- xlim[2]

    ra <- sum(oob_arrows$ci_upper_bound > oob_arrows$x_high, na.rm = T) > 0
    la <- sum(oob_arrows$ci_lower_bound < oob_arrows$x_low, na.rm = T) > 0

    if (ra) {
      right_arrows <- dplyr::select(dplyr::filter(oob_arrows, ci_upper_bound > .data$x_high), start = point_estimate, end = .data$x_high, y = .data$row_num)
    }
    if (la) {
      left_arrows <- dplyr::select(dplyr::filter(oob_arrows, ci_lower_bound < .data$x_low), start = point_estimate, end = .data$x_low, y = .data$row_num)
    }

    if (ra && !la) {
      g_oob <- right_arrows
    }else if (!ra && la) {
      g_oob <- left_arrows
    }else if (ra && la) {
      g_oob <- rbind.data.frame(right_arrows, left_arrows)
    }
  }

  ########## the main figure - this will be overlaid on the table ##############
  p_mid <- ggplot2::ggplot() +
    ggplot2::geom_point(data = gdata, ggplot2::aes(y = row_num, x = point_estimate, size = sizes, shape = shape), na.rm = TRUE, color = gdata$effect_color) +
    ggplot2::geom_errorbarh(data = gdata, ggplot2::aes(y = row_num,
          xmin = ci_lower_bound,
          xmax = ci_upper_bound),
          height = .25,
          color = gdata$effect_color,
          na.rm = TRUE) +
    ggplot2::theme_classic() + # base theme
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), # remove axis, make bg transparent
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          axis.ticks.length.x = grid::unit(.07, "in"),
          text = ggplot2::element_text(family = font_family, size = 12),
          panel.background = ggplot2::element_rect(fill = "transparent"),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "transparent"),
          legend.box.background = ggplot2::element_rect(fill = "transparent")) +
    ggplot2::geom_vline(xintercept = null_line_at, linetype = "dashed") +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_size_identity() +
      ggplot2::xlab("")

  ### add oob arrows if required ###
  if (nrow(g_oob) > 0) {
    p_mid <- p_mid +
      ggplot2::geom_segment(data = g_oob,
                            ggplot2::aes(x = start,
                                             xend = end,
                                             y = y,
                                             yend = y),
                            arrow = ggplot2::arrow(angle = 15,
                                             type = "closed",
                                             length = grid::unit(0.1, "in")))
  }

  ####### fix plot zoom ######
  if (is.null(xlim)) {
    p_mid <- p_mid + ggplot2::coord_cartesian(ylim = c(y_low, y_high))
  }else{
    p_mid <- p_mid + ggplot2::coord_cartesian(ylim = c(y_low, y_high), xlim = xlim)
  }

  ######## handle breaks, log vs linear scales ########
  if (x_scale_linear) {
    if (is.null(xbreaks)) {
      p_mid <- p_mid + ggplot2::scale_x_continuous(labels = scales::number_format(accuracy = 0.1),
                                            expand = c(0,0))
    }else{
      p_mid <- p_mid + ggplot2::scale_x_continuous(labels = scales::number_format(accuracy = 0.1),
                                            breaks = xbreaks,
                                            expand = c(0,0))
    }
  }else{
    if (is.null(xbreaks)) {
      p_mid <- p_mid + ggplot2::scale_x_log10(labels = scales::number_format(accuracy = 0.1),
                                            expand = c(0,0))
    }else{
      p_mid <- p_mid + ggplot2::scale_x_log10(labels = scales::number_format(accuracy = 0.1),
                                            breaks = xbreaks,
                                            expand = c(0,0))
    }
  }

  #### allow overwrite of central plot #######################
  if(!is.null(p_mid_forest)) {p_mid <- p_mid_forest}

  ######################## Arrows ##############################
  if (arrows == TRUE) {
    # this df has the text labels
    xlab_df <- data.frame(text = arrow_labels,
                          x = xlim,
                          y = c(0, 0) + arrow_nudge_y,
                          hjust = c(0, 1))
    a_small_amount <- abs(xlim[1] - xlim[2])/35

    # this df has the arrows
    if (x_scale_linear == TRUE) {
      arrow_df <- data.frame(id = c(1,2),
                           xstart = c(null_line_at - a_small_amount, null_line_at + a_small_amount),
                           xend = c(xlim[1] + a_small_amount, xlim[2] - a_small_amount),
                           y = c(1, 1) + arrow_nudge_y)
    }else{
      arrow_df <- data.frame(id = c(1,2),
                             xstart = c(null_line_at - a_small_amount, null_line_at + a_small_amount),
                             xend = c(xlim[1], xlim[2]),
                             y = c(1, 1) + arrow_nudge_y)
    }
    # create the arrow/label ggplot object
    arrows_plot <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = arrow_df, ggplot2::aes(x = .data$xstart, xend = .data$xend, y = .data$y, yend = .data$y),
                 arrow = ggplot2::arrow(angle = 15, type = "closed", length = grid::unit(0.1, "in"))) +
      ggplot2::geom_text(data = xlab_df, ggplot2::aes(x = .data$x, y = .data$y, label = .data$text, hjust = .data$hjust),
              family = font_family, size = 3) +
      ggplot2::scale_y_continuous(expand = c(0,0), limits = c(-0.5, 1.75)) +
      ggplot2::scale_x_continuous(expand = c(0,0), limits = xlim) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent"),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "transparent"),
          legend.box.background = ggplot2::element_rect(fill = "transparent"),
          panel.border = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank())

    if (x_scale_linear == FALSE) {
      arrows_plot <- arrows_plot + ggplot2::scale_x_log10(expand = c(0,0), limits = xlim)
    }

  }


  ######### using patchwork, overlay the ggplot on the table ###################
  png_width <- total_width/10 + nudge_x
  png_height <- (nrow(gdata) + 3)/3.8
  if (is.null(add_plot)) {
    table_final <- mono_column(gridExtra::tableGrob(tdata_print, theme = theme, rows = NULL), ncol(p_left_data) + 1, font_family)

    table_final$widths[ncol(p_left_data) + 1] <- grid::unit(p_mid_width/10, "in")
    table_final$heights <- grid::unit(rep(0.255, times = length(table_final$heights)), "in")

    final <- patchwork::wrap_elements(table_final) +
                    patchwork::inset_element(p_mid,
                             align_to = "full",
                             left = (left_width/total_width),
                             right = ((p_mid_width + left_width)/total_width),
                             top = 1,
                             bottom = 0.35/nrow(gdata))

    if (arrows == TRUE) {
      final <- final + patchwork::inset_element(arrows_plot,
                                              align_to = "full",
                                              left = (left_width/total_width),
                                              right = ((p_mid_width + left_width)/total_width),
                                              top = 1.5/nrow(gdata),
                                              bottom = 0)
    }

  }else{
    tdata_print$`  ` <- paste(rep(" ", times = round(p_mid_width, 0)),
                              collapse = '')

    table_final <- mono_column(gridExtra::tableGrob(tdata_print, theme = theme, rows = NULL), ncol(p_left_data) + 1, font_family)

    table_final <- white_column(table_final, ncol(table_final), font_family)

    table_final$widths[ncol(p_left_data) + 1] <- grid::unit(p_mid_width/10, "in")
    table_final$widths[ncol(table_final)] <- grid::unit(p_mid_width/10, "in")

    table_final$heights <- grid::unit(rep(0.255, times = length(table_final$heights)), "in")

    new_full_width <- total_width + p_mid_width

    png_width <- new_full_width/10 + nudge_x

    if (add_plot_gap){
      add_plot <- add_plot + ggplot2::scale_y_continuous(limits = c(y_low, y_high), expand = c(0,0)) +
        ggplot2::theme_classic() + # base theme
        ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "transparent"), # make axis transparent rather than removing it
                     axis.text.y = ggplot2::element_text(colour = "transparent"), # this makes alignment much easier
                     axis.ticks.y = ggplot2::element_line(colour = "transparent"),
                     axis.line.y = ggplot2::element_line(colour = "transparent"),
                     axis.title.x = ggplot2::element_text(colour = "transparent"),
                     axis.text.x = ggplot2::element_text(colour = "transparent"),
                     axis.ticks.x = ggplot2::element_line(colour = "transparent"),
                     axis.line.x = ggplot2::element_line(colour = "transparent"),
                     panel.background = ggplot2::element_rect(fill = "transparent"),
                     plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"),
                     legend.box.background = ggplot2::element_rect(fill = "transparent"),
                     legend.position = "none")
    }else{
      add_plot <- add_plot + ggplot2::scale_y_continuous(limits = c(y_low, y_high), expand = c(0,0)) +
        ggplot2::theme_classic() + # base theme
        ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "transparent"), # make x axis (only) transparent
                       axis.text.x = ggplot2::element_text(colour = "transparent"),
                       axis.ticks.x = ggplot2::element_line(colour = "transparent"),
                       axis.line.x = ggplot2::element_line(colour = "transparent"),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "transparent"),
                       plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(fill = "transparent"),
                       legend.box.background = ggplot2::element_rect(fill = "transparent"),
                       legend.position = "none")
    }

    final <- patchwork::wrap_elements(table_final) +
      patchwork::inset_element(p_mid,
                               align_to = "full",
                               left = (left_width/new_full_width),
                               right = ((p_mid_width + left_width)/new_full_width),
                               top = 1,
                               bottom = 0.35/nrow(gdata)) +
      patchwork::inset_element(add_plot,
                               align_to = "full",
                               left = total_width/new_full_width,
                               right = 1,
                               top = 1,
                               bottom = 0.35/nrow(gdata))

    if (arrows == TRUE) {
      final <- final + patchwork::inset_element(arrows_plot,
                                                align_to = "full",
                                                left = (left_width/new_full_width),
                                                right = ((p_mid_width + left_width)/new_full_width),
                                                top = 1.5/nrow(gdata),
                                                bottom = 0)
    }
  }

  ######### save the plot as a png, then display it with magick ################
  if (!(render_as == "rmarkdown")) {
    if(is.null(output_path)){
      output_path <- tempdir()
    }
    filename = file.path(output_path, paste0("forest_plot.", render_as))
    message(paste0("Saving image to ", filename))
    ggplot2::ggsave(
         dpi = dpi,
         height = png_height + nudge_height,
         width = png_width + nudge_width,
         units = "in",
         filename = filename,
         device = render_as
    )

    if (display == TRUE) {
      invisible(system(paste0('open "', filename, '"')))
    }
  }else{
    final
  }
}
