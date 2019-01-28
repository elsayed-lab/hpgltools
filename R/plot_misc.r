## plot_misc.r:  Silly plots and plotting helpers.

#' Plot a picture, with hopefully useful options for most(any) format.
#'
#' This calls svg/png/postscript/etc according to the filename provided.
#'
#' @param file Filename to write
#' @param image Optionally, add the image you wish to plot and this will both
#'   print it to file and screen.
#' @param width  How wide?
#' @param height  How high?
#' @param res  The chosen resolution.
#' @param ...  Arguments passed to the image plotters.
#' @return a png/svg/eps/ps/pdf with height=width=9 inches and a high resolution
#' @export
pp <- function(file, image=NULL, width=9, height=9, res=180, ...) {
  ext <- tools::file_ext(file)
  start_dev <- dev.list()
  result <- NULL
  if (ext == "png") {
    result <- png(filename=file, width=width, height=height, units="in", res=res, ...)
  } else if (ext == "svg") {
    result <- svg(filename=file, ...)
  } else if (ext == "ps") {
    result <- postscript(file=file, width=width, height=height, ...)
  } else if (ext == "eps") {
    result <- cairo_ps(filename=file, width=width, height=height, ...)
  } else if (ext == "pdf") {
    result <- cairo_pdf(filename=file, ...)
  } else {
    message("Defaulting to tiff.")
    result <- tiff(filename=file, width=width, height=height, units="in", res=res, ...)
  }
  now_dev <- dev.list()
  new_dev <- now_dev[length(now_dev)]

  ## Check and make sure I am not looking at something containing a plot, as a bunch of
  ## my functions are lists with a plot slot.
  if (class(image)[[1]] == "list") {
    if (!is.null(image[["plot"]])) {
      image <- image[["plot"]]
    }
  }

  if (is.null(image)) {
    message("Going to write the image to: ", file, " when dev.off() is called.")
    return(invisible(image))
  } else {
    message("Writing the image to: ", file, " and calling dev.off().")
    if (class(image)[[1]] == "recordedplot") {
      print(image)
    } else {
      plot(image)
    }
    dev.off(which=new_dev)
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
#' @param radius_a   The radius of the primary circle.
#' @param radius_b   The radius of the circle travelling around a.
#' @param dist_bc   A point relative to the center of 'b' which rotates with the turning of 'b'.
#' @param revolutions   How many revolutions to perform in the plot
#' @param increments   The number of radial increments to be calculated per revolution
#' @param center_a   The position of the center of 'a'.
#' @return something which I don't yet know.
#' @export
plot_spirograph <- function(radius_a=1, radius_b=-4, dist_bc=-2,
                            revolutions=158, increments=3160, center_a=list(x=0, y=0)) {
  center_b_start <- list(x=0, y=center_a$y + radius_a + radius_b)
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
  spiro <- ggplot2::ggplot(data=points, ggplot2::aes_string(x="x", y="y")) +
    ggplot2::geom_point(ggplot2::aes_string(colour="counter"), size=0.5) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_gradientn(colours=grDevices::rainbow(4)) +
    ggplot2::theme(axis.line=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   legend.position="none",
                   panel.background=ggplot2::element_blank(),
                   panel.border=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   plot.background=ggplot2::element_blank())
  return(spiro)
}

#' Make hypotrochoid plots!
#'
#' 3,7,1 should give the classic 7 leaf clover
#'
#' @param radius_a  Radius of the major circle
#' @param radius_b  And the smaller circle.
#' @param dist_b between b and the drawing point.
#' @param revolutions  How many times to revolve through the spirograph.
#' @param increments  How many dots to lay down while writing.
plot_hypotrochoid <- function(radius_a=3, radius_b=7, dist_b=1,
                              revolutions=7, increments=6480) {
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
  image <- ggplot2::ggplot(data=positions, ggplot2::aes_string(x="x_points", y="y_points")) +
    ggplot2::geom_point(size=1) + ggplot2::theme_bw() +
    ggplot2::theme(axis.line=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   legend.position="none",
                   panel.background=ggplot2::element_blank(),
                   panel.border=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   plot.background=ggplot2::element_blank())
  return(image)
}

#' Make epitrochoid plots!
#'
#' 7, 2, 6, 7 should give a pretty result.
#'
#' @param radius_a  Radius of the major circle
#' @param radius_b  And the smaller circle.
#' @param dist_b between b and the drawing point.
#' @param revolutions  How many times to revolve through the spirograph.
#' @param increments  How many dots to lay down while writing.
plot_epitrochoid <- function(radius_a=7, radius_b=2, dist_b=6,
                             revolutions=7, increments=6480) {
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
  image <- ggplot2::ggplot(data=positions, ggplot2::aes_string(x="x_points", y="y_points")) +
    ggplot2::geom_point(size=1) + ggplot2::theme_bw() +
    ggplot2::theme(axis.line=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   legend.position="none",
                   panel.background=ggplot2::element_blank(),
                   panel.border=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   plot.background=ggplot2::element_blank())
  return(image)
}

## EOF
