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
  } else {
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
  }

  return(image)
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
  spiro <- ggplot2::ggplot(data = points, ggplot2::aes_string(x = "x", y = "y")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "counter"), size = 0.5) +
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
  image <- ggplot2::ggplot(data = positions, ggplot2::aes_string(x = "x_points", y = "y_points")) +
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
  image <- ggplot2::ggplot(data = positions, ggplot2::aes_string(x = "x_points", y = "y_points")) +
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
  if (is.null(tooltip) | tooltip == "all") {
    tooltip_columns <- "label"
  } else {
    first_tooltip_column <- tooltip[1]
  }

  if (is.null(plot[["data"]][[id_column]])) {
    plot[["data"]][[id_column]] <- rownames(plot[["data"]])
  }
  if (is.null(url_info) & is.null(url_column)) {
    warning("No url information was provided.")
  } else if ("character" %in% class(url_info) & length(url_info) > 1) {
    message("url_info has multiple entries, assuming it is a character vector with 1 url/entry.")
    ## Then this should contain all the URLs
    plot[["data"]][[url_column]] <- url_info
  } else if ("glue" %in% class(url_info) & length(url_info) == 1) {
    message("url_info has length 1, assuming it is a glue specification including {ids}.")
    ## Assuming url_data looks like: 'http://useast.ensembl.org/Mus_musculus/Gene/Summary?q={ids}'
    ids <- plot[["data"]][[id_column]]
    plot[["data"]][[url_column]] <- glue::glue(url_info)
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
  for (i in 1:length(plotly[["x"]][["data"]])) {
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
