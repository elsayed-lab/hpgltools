## plot_misc.r:  Silly plots and plotting helpers.

#' Plot a picture, with hopefully useful options for most(any) format.
#'
#' This calls svg/png/postscript/etc according to the filename provided.
#'
#' @param file Filename to write
#' @param image Optionally, add the image you wish to plot and this will both
#'  print it to file and screen.
#' @param width How wide?
#' @param height How high?
#' @param res The chosen resolution.
#' @param ... Arguments passed to the image plotters.
#' @return a png/svg/eps/ps/pdf with height = width=9 inches and a high resolution
#' @seealso [png()] [svg()] [postscript()] [cairo_ps()] [cairo_pdf()] [tiff()] [devEMF::emf()]
#'  [jpg()] [bmp()]
#' @export
pp <- function(file, image = NULL, width = 9, height = 9, res = 180, ...) {
  ext <- tolower(tools::file_ext(file))
  file_dir <- dirname(file)
  if (!file.exists(file_dir)) {
    warning("The directory: ", file_dir, " does not exist, will attempt to create it.")
    dir.create(file_dir, recursive = TRUE)
  }

  start_dev <- dev.list()
  result <- NULL
  switchret <- switch(
    ext,
    "png" = {
      result <- png(filename = file, width = width, height = height,
                    units = "in", res = res,
                    ...)
    },
    "bmp" = {
      result <- bmp(filename = file, ...)
    },
    "jpg" = {
      result <- jpeg(filename = file, ...)
    },
    "webp" = {
      result <- webp::write_webp(target = file, ...)
    },
    "svg" = {
      result <- svg(filename = file, ...)
    },
    "ps" = {
      result <- postscript(file = file, width = width, height = height, ...)
    },
    "eps" = {
      result <- cairo_ps(filename = file, width = width, height = height, ...)
    },
    "pdf" = {
      result <- cairo_pdf(filename = file, ...)
    },
    "tif" = {
      result <- tiff(filename = file, width = width, height = height,
                     units = "in", res = res, ...)
    },
    "tiff" = {
      result <- tiff(filename = file, width = width, height = height,
                     units = "in", res = res, ...)
    },
    "emf" = {
      result <- devEMF::emf(file = file, width = width, height = height, ...)
    },
    {
      mesg("Defaulting to tiff.")
      result <- tiff(filename = file, width = width, height = height,
                     units = "in", res = res, ...)
    }) ## End of the switch
  ## Find the new device for closing later.
  now_dev <- dev.list()
  new_dev_idx <- ! names(now_dev) %in% names(start_dev)
  new_dev <- now_dev[new_dev_idx]

  ## Check and make sure I am not looking at something containing a plot, as a bunch of
  ## my functions are lists with a plot slot.
  if (class(image)[[1]] == "list") {
    if (!is.null(image[["plot"]])) {
      image <- image[["plot"]]
    }
  }

  if (is.null(image)) {
    mesg("Going to write the image to: ", file, " when dev.off() is called.")
    mesg("Do not forget to close the device when you are done.")
    return(invisible(result))
  }


  mesg("Writing the image to: ", file, " and calling dev.off().")
  if (class(image)[[1]] == "recordedplot") {
    print(image)
  } else {
    plot(image)
  }

  if (length(new_dev) > 0) {
    dev.off(which = new_dev)
  } else {
    warning("There is no device to shut down.")
  }
  return(image)
}

#' Plot metadata factors as a sankey diagram.
#'
#' This provides two implementations of a sankey plot, one interactive and one using ggplot2.
#'
#' @param design Metadata from which to extract the categories/numbers.
#' @param factors Factors/columns in the metadata to count and plot.
#' @param fill Use either the current or next node for coloring the transitions.
#' @param color_choices Either a named vector of states and colors, or NULL
#'  (in which case it will use viridis.)
#' @param html Write the interactive plot to this file.
#' @param drill_down When true, this will end in the product of the
#'  factor levels number of final states. (e.g. if there are 2 sexes,
#'  3 visits, and 4 genotypes, there will be 2, 6, 24 states going
#' from right to left).  If FALSE, there will be 2,3,4 states going
#' from right to left.
#' @return List containing a couple of plots, one interactive, one gg.
#' @export
plot_meta_sankey <- function(design, factors = c("condition", "batch"), fill = "node",
                             font_size = 18, node_width = 30,
                             color_choices = NULL,
                             drill_down = TRUE) {
  warning("FIXME: I separated the interactive and ggplot functions, but haven't figured out what need to be kept.")
  found <- factors %in% colnames(design)
  if (sum(found) < length(factors)) {
    missing <- factors[!found]
    message("These columns are not in the metadata: ", toString(missing))
    factors <- factors[found]
  }

  permutations <- c()
  states <- list()
  for (f in factors) {
    state_levels <- levels(as.factor(design[[f]]))
    states[[f]] <- state_levels
    new_permutations <- tidyr::expand_grid(!!!states) |> purrr::pmap_chr(paste)
    permutations <- c(permutations, new_permutations)
  }

  start_node <- c("all", "0")
  names(start_node) <- c("name", "node")
  my_nodes <- data.frame(name = permutations)
  my_nodes[["node"]] <- rownames(my_nodes)
  my_nodes <- rbind(start_node, my_nodes)

  my_links <- data.frame(row.names = 0)
  my_links[["source"]] <- 0
  my_links[["target"]] <- 0
  my_links[["value"]] <- nrow(design)

  result <- list("all" = nrow(design))
  for (p in seq_len(length(permutations))) {
    element <- permutations[p]
    pieces <- strsplit(x = element, split = " ")[[1]]

    working_meta <- design
    for (cat_num in seq_len(length(pieces))) {
      category <- pieces[cat_num]
      factor <- factors[cat_num]
      idx <- working_meta[[factor]] == category
      working_meta <- working_meta[idx, ]
    }
    if (nrow(working_meta) > 0) {
      result[[element]] <- nrow(working_meta)

      target_node_idx <- my_nodes[["name"]] == element
      target_node <- my_nodes[target_node_idx, "node"]
      if (length(pieces) > 1) {
        source_node_pieces <- pieces[1:length(pieces) - 1]
        source_node_name <- stringi::stri_paste(source_node_pieces, collapse = " ")
        source_node_idx <- my_nodes[["name"]] == source_node_name
        source_node <- my_nodes[source_node_idx, "node"]
      } else {
        source_node <- "0"
      }
      link <- c(source_node, target_node, nrow(working_meta))
      my_links <- rbind(my_links, link)
    }
  }

  my_links <- my_links[-1, ]
  my_links[["value"]] <- as.numeric(my_links[["value"]])
  my_links[["source"]] <- as.numeric(my_links[["source"]])
  my_links[["target"]] <- as.numeric(my_links[["target"]])

  links_to_nodes <- merge(my_links, my_nodes, by.x = "target",
                          by.y = "node")
  links_to_nodes <- links_to_nodes[, c("name", "value")]

  sub_design <- design[, c(factors)]
  if (isTRUE(drill_down)) {
    for (col in seq_len(ncol(sub_design))) {
      if (col == 1) {
        next
      }
      col_name <- colnames(sub_design)[col]
      previous_col <- colnames(sub_design)[col - 1]
      new_values <- paste0(sub_design[[previous_col]], " ",
                           sub_design[[col_name]])
      sub_design[[col_name]] <- new_values
    }
  }

  plot_df <- sub_design %>%
    ggsankey::make_long(factors)
  plot_df[["name"]] <- plot_df[["node"]]
  plot_df <- merge(plot_df, links_to_nodes, by.x = "node", by.y = "name")
  plot_df[["name"]] <- paste0(plot_df[["name"]], ":",
                              plot_df[["value"]])

  color_fact <- NULL
  if (!is.null(color_choices)) {
    color_levels <- levels(as.factor(plot_df[["node"]]))
    all_colors <- unlist(color_choices)
    names(all_colors) <- gsub(x = names(all_colors),
                              pattern = ".*\\.", replacement = "")
    color_fact_idx <- names(all_colors) %in% color_levels
    color_fact <- all_colors[color_fact_idx]
    if (isTRUE(drill_down)) {
      ## Set the plot name and color names:
      plot_df[["name"]] <- gsub(x = plot_df[["node"]],
        pattern = "^.* (\\w+$)",
        replacement = "\\1")
      plot_df[["name"]] <- paste0(plot_df[["name"]], ":",
                                  plot_df[["value"]])
      color_level_suffixes <- gsub(x = color_levels,
        pattern = "^.* (\\w+$)",
        replacement = "\\1")
      color_fact_idx <- names(all_colors) %in% color_level_suffixes
      color_suffix_fact <- all_colors[color_fact_idx]
      new_color_fact <- rep("#000000", length(color_level_suffixes))
      names(new_color_fact) <- color_levels
      for (col in seq_len(length(new_color_fact))) {
        color_name <- color_level_suffixes[col]
        color_suffix <- as.character(color_suffix_fact[color_name])
        new_color_fact[col] <- color_suffix
      }
      color_fact <- new_color_fact
    }
  }
  retlist <- list(
    "design" = design,
    "factor" = factors,
    "observed_nodes" = unique(plot_df[["node"]]))

  if (fill == "node") {
    ggplt <- ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                                 next_node = next_node, fill = factor(node), label = name))
  } else if (fill == "next") {
    message("Filling to next node?")
    ggplt <- ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                                 next_node = next_node, fill = factor(next_node), label = name))
  }
  ggplt <- ggplt +
    ggsankey::geom_sankey(flow.alpha = 0.6,
                          node.color = "gray30") +
    ggsankey::geom_sankey_label()

  ## I want to figure out how to set up my own colors...
  if (is.null(color_choices)) {
    ggplt <- ggplt +
      ggplot2::scale_fill_viridis_d() +
      ggsankey::theme_sankey(base_size = 18) +
      ggplot2::labs(x = NULL) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
  } else {
    ggplt <- ggplt +
      ggplot2::scale_fill_manual(values = color_fact) +
      ggsankey::theme_sankey(base_size = font_size) +
      ggplot2::labs(x = NULL) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
  }
  retlist[["ggplot"]] <- ggplt
  class(retlist) <- "meta_sankey"
  return(retlist)
}

## I might want to just delete this, I think the ggplot version is better.
plot_meta_interactive_sankey <- function(design, factors = c("condition", "batch"), fill = "node",
                                         font_size = 18, node_width = 30,
                                         color_choices = NULL, html = NULL,
                                         drill_down = TRUE) {

  found <- factors %in% colnames(design)
  if (sum(found) < length(factors)) {
    missing <- factors[!found]
    message("These columns are not in the metadata: ", toString(missing))
    factors <- factors[found]
  }

  permutations <- c()
  states <- list()
  for (f in factors) {
    state_levels <- levels(as.factor(design[[f]]))
    states[[f]] <- state_levels
    new_permutations <- tidyr::expand_grid(!!!states) |> purrr::pmap_chr(paste)
    permutations <- c(permutations, new_permutations)
  }

  start_node <- c("all", "0")
  names(start_node) <- c("name", "node")
  my_nodes <- data.frame(name = permutations)
  my_nodes[["node"]] <- rownames(my_nodes)
  my_nodes <- rbind(start_node, my_nodes)

  my_links <- data.frame(row.names = 0)
  my_links[["source"]] <- 0
  my_links[["target"]] <- 0
  my_links[["value"]] <- nrow(design)

  result <- list("all" = nrow(design))
  for (p in seq_len(length(permutations))) {
    element <- permutations[p]
    pieces <- strsplit(x = element, split = " ")[[1]]

    working_meta <- design
    for (cat_num in seq_len(length(pieces))) {
      category <- pieces[cat_num]
      factor <- factors[cat_num]
      idx <- working_meta[[factor]] == category
      working_meta <- working_meta[idx, ]
    }
    if (nrow(working_meta) > 0) {
      result[[element]] <- nrow(working_meta)

      target_node_idx <- my_nodes[["name"]] == element
      target_node <- my_nodes[target_node_idx, "node"]
      if (length(pieces) > 1) {
        source_node_pieces <- pieces[1:length(pieces) - 1]
        source_node_name <- stringi::stri_paste(source_node_pieces, collapse = " ")
        source_node_idx <- my_nodes[["name"]] == source_node_name
        source_node <- my_nodes[source_node_idx, "node"]
      } else {
        source_node <- "0"
      }
      link <- c(source_node, target_node, nrow(working_meta))
      my_links <- rbind(my_links, link)
    }
  }

  my_links <- my_links[-1, ]
  my_links[["value"]] <- as.numeric(my_links[["value"]])
  my_links[["source"]] <- as.numeric(my_links[["source"]])
  my_links[["target"]] <- as.numeric(my_links[["target"]])

  links_to_nodes <- merge(my_links, my_nodes, by.x = "target",
                          by.y = "node")
  links_to_nodes <- links_to_nodes[, c("name", "value")]


  plt <- NULL
  if (!is.null(html)) {
    plt = networkD3::sankeyNetwork(Links = my_links, Nodes = my_nodes,
                                   Source = "source", Target = "target",
                                   Value = "value", fontSize = font_size, nodeWidth = node_width)
  }
  retlist <- list(
    "permutations" = permutations,
    "plot" = plt)
  if (!is.null(html)) {
    retlist[["html"]] <- htmlwidgets::saveWidget(plt, file = html, selfcontained = TRUE)
  }
  return(retlist)
}

#' Make spirographs!
#'
#' Taken (with modifications) from:
#' http://menugget.blogspot.com/2012/12/spirograph-with-r.html#more
#' A positive value for 'B' will result in a epitrochoid, while a negative value
#' will result in a hypotrochoid.
#'
#' @param radius_a The radius of the primary circle.
#' @param radius_b The radius of the circle travelling around a.
#' @param dist_bc A point relative to the center of 'b' which rotates with the turning of 'b'.
#' @param revolutions How many revolutions to perform in the plot
#' @param increments The number of radial increments to be calculated per revolution
#' @param center_a The position of the center of 'a'.
#' @return something which I don't yet know.
#' @export
plot_spirograph <- function(radius_a = 1, radius_b=-4, dist_bc=-2,
                            revolutions = 158, increments = 3160, center_a = list(x = 0, y = 0)) {
  center_b_start <- list(x = 0, y = center_a$y + radius_a + radius_b)
  angle_a <- seq(0, 2 * pi * revolutions, revolutions * increments)
  circum_a <- 2 * pi * radius_a
  circum_b <- 2 * pi * radius_b
  center_b <- c()
  hypotenuse <- radius_a + radius_b
  adjacent <- sin(angle_a) * hypotenuse
  opposite <- cos(angle_a) * hypotenuse
  center_b[["x"]] <- center_a[["x"]] + adjacent
  center_b[["y"]] <- center_a[["y"]] + opposite
  point_c <- c()
  circle_a_dist <- (circum_a * angle_a) / (2 * pi)
  angle_b_point <- circle_a_dist / (circum_b * (2 * pi))
  hypotenuse <- dist_bc
  adjacent <- sin(angle_b_point) * hypotenuse
  opposite <- cos(angle_b_point) * hypotenuse
  point_c[["x"]] <- center_b[["x"]] + adjacent
  point_c[["y"]] <- center_b[["y"]] + opposite
  points <- data.frame(point_c)
  points[["counter"]] <- seq(1, nrow(points))
  spiro <- ggplot(data = points, aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_point(ggplot2::aes(colour = .data[["counter"]]), size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_gradientn(colours = grDevices::rainbow(4)) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank())
  return(spiro)
}

#' Make hypotrochoid plots!
#'
#' 3,7,1 should give the classic 7 leaf clover
#'
#' @param radius_a Radius of the major circle
#' @param radius_b And the smaller circle.
#' @param dist_b between b and the drawing point.
#' @param revolutions How many times to revolve through the spirograph.
#' @param increments How many dots to lay down while writing.
#' @export
plot_hypotrochoid <- function(radius_a = 3, radius_b = 7, dist_b = 1,
                              revolutions = 7, increments = 6480) {
  points <- seq(0, revolutions * increments)
  radians <- points / (2 * pi)
  getx <- function(t) {
    x <- ((radius_a - radius_b) * cos(t)) +
      (dist_b * cos((t * ((radius_a - radius_b) / radius_b))))
    return(x)
  }
  gety <- function(t) {
    y <- ((radius_a - radius_b) * sin(t)) +
      (dist_b * sin((t * ((radius_a - radius_b) / radius_b))))
    return(y)
  }
  x_points <- as.numeric(lapply(radians, getx))
  y_points <- as.numeric(lapply(radians, gety))
  positions <- cbind(points, x_points)
  positions <- cbind(positions, y_points)
  positions <- as.data.frame(positions)
  petals <- dist_b / radius_b
  message("The spirograph will have ", petals, " petals.")
  image <- ggplot(
    data = positions, aes(x = .data[["x_points"]], y = .data[["y_points"]])) +
    ggplot2::geom_point(size = 1) + ggplot2::theme_bw() +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank())
  return(image)
}

#' Make epitrochoid plots!
#'
#' 7, 2, 6, 7 should give a pretty result.
#'
#' @param radius_a Radius of the major circle
#' @param radius_b And the smaller circle.
#' @param dist_b between b and the drawing point.
#' @param revolutions How many times to revolve through the spirograph.
#' @param increments How many dots to lay down while writing.
#' @export
plot_epitrochoid <- function(radius_a = 7, radius_b = 2, dist_b = 6,
                             revolutions = 7, increments = 6480) {
  points <- seq(0, revolutions * increments)
  radians <- points / (2 * pi)
  getx <- function(t) {
    x <- ((radius_a + radius_b) * cos(t)) -
      (dist_b * cos((t * ((radius_a + radius_b) / radius_b))))
    return(x)
  }
  gety <- function(t) {
    y <- ((radius_a + radius_b) * sin(t)) -
      (dist_b * sin((t * ((radius_a + radius_b) / radius_b))))
    return(y)
  }
  x_points <- as.numeric(lapply(radians, getx))
  y_points <- as.numeric(lapply(radians, gety))
  positions <- cbind(points, x_points)
  positions <- cbind(positions, y_points)
  positions <- as.data.frame(positions)
  image <- ggplot(data = positions,
                  aes(x = .data[["x_points"]], y = .data[["y_points"]])) +
    ggplot2::geom_point(size = 1) + ggplot2::theme_bw() +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank())
  return(image)
}

#' Add a little logic to ggplotly to simplify adding clicky link.
#'
#' There are some other ease of life improvements I have in a few of my plotly
#' invocations which I should add here.
#'
#' @param plot Plot generated via ggplot2.
#' @param filename filename to save the output html plot.
#' @param id_column Column containing the gene IDs.
#' @param plot_title Provide a title for the generated html file.
#' @param url_info Either a glue() string or column of urls.
#' @param tooltip Passed to ggplotly().
#' @param url_column Column in the url_info containing URLs.
#' @return plotly with clicky links.
#' @export
ggplotly_url <- function(plot, filename = "ggplotly_url.html", id_column = "id", plot_title = NULL,
                         url_info = NULL, tooltip = "all", url_column = "url") {
  first_tooltip_column <- "label"
  if (is.null(tooltip) || tooltip == "all") {
    tooltip_columns <- "label"
  } else {
    first_tooltip_column <- tooltip[1]
  }

  if (is.null(plot[["data"]][[id_column]])) {
    plot[["data"]][[id_column]] <- rownames(plot[["data"]])
  }
  if (is.null(url_info) && is.null(url_column)) {
    warning("No url information was provided.")
  } else if ("character" %in% class(url_info) && length(url_info) > 1) {
    message("url_info has multiple entries, assuming it is a character vector with 1 url/entry.")
    ## Then this should contain all the URLs
    plot[["data"]][[url_column]] <- url_info
  } else if ("glue" %in% class(url_info) && length(url_info) == 1) {
    message("url_info has length 1, assuming it is a glue specification including {ids}.")
    ## Assuming url_data looks like: 'http://useast.ensembl.org/Mus_musculus/Gene/Summary?q={ids}'
    ids <- plot[["data"]][[id_column]]
    plot[["data"]][[url_column]] <- glue(url_info)
  } else if ("data.frame" %in% class(url_info)) {
    ## This assumes url data has a column named whatever is url_column
    message("Merging the url data with the plot data.")
    plot[["data"]] <- merge(plot[["data"]], url_info, by = "row.names", all.x = TRUE)
    rownames(plot[["data"]]) <- plot[["Row.names"]]
    plot[["data"]][["Row.names"]] <- NULL
  }

  if (is.null(plot[["data"]][[first_tooltip_column]])) {
    plot[["data"]][[first_tooltip_column]] <- rownames(plot[["data"]])
  }
  plotly <- plotly::ggplotly(plot, tooltip = tooltip)
  for (i in seq_along(plotly[["x"]][["data"]])) {
    plotly[["x"]][["data"]][[i]][["customdata"]] <- plot[["data"]][[url_column]]
  }
  plotly <- htmlwidgets::onRender(plotly, "
function(el, x) {
  el.on('plotly_click', function(d) {
    var url = d.points[0].customdata;
    console.log(url);
    window.open(url);
  });
}")
  out <- htmlwidgets::saveWidget(plotly, filename, title = plot_title)
  retlist <- list(
    "out" = out,
    "plotly" = plotly,
    "modified_df" = plot[["data"]])
  return(retlist)
}

## EOF
