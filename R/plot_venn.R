#' A quick wrapper around venneuler to help label stuff
#'
#' venneuler makes pretty venn diagrams, but no labels!
#'
#' @param ones Character list of singletone categories
#' @param twos Character list of doubletone categories
#' @param threes Character list of tripletone categories
#' @param fours Character list of quad categories
#' @param fives Character list of quint categories
#' @param factor Currently unused, but intended to change the radial distance to
#'  the label from the center of each circle.
#' @return Two element list containing the venneuler data and the plot.
#' @seealso [venneuler]
#' @export
plot_fun_venn <- function(ones = c(), twos = c(), threes = c(),
                          fours = c(), fives = c(), factor = 0.9) {
  venn_sets <- ones
  venn_intersect_label <- ""
  do_doubles <- FALSE
  do_triples <- FALSE
  do_quads <- FALSE

  if (is.null(venn_sets)) {
    stop("I received nothing to play with.")
  }

  if (!is.null(twos)) {
    venn_sets <- append(venn_sets, twos)
  } else {
    if (!is.null(ones)) {
      venn_intersect_label <- as.character(ones[[length(ones)]])
    }
  }

  if (!is.null(threes)) {
    venn_sets <- append(venn_sets, threes)
    do_doubles <- TRUE
  } else {
    if (!is.null(twos)) {
      venn_intersect_label <- as.character(twos[[length(twos)]])
    }
  }

  if (!is.null(fours)) {
    venn_sets <- append(venn_sets, fours)
    do_triples <- TRUE
  } else {
    if (!is.null(threes)) {
      venn_intersect_label <- as.character(threes[[length(threes)]])
    }
  }

  if (!is.null(fives)) {
    venn_sets <- append(venn_sets, fives)
    do_quads <- TRUE
  } else {
    if (!is.null(fours)) {
      venn_intersect_label <- as.character(fours[[length(fours)]])
    }
  }

  all_venn <- venneuler::venneuler(venn_sets)
  plot(all_venn)

  center_x <- mean(all_venn[["centers"]][, 1])
  center_y <- mean(all_venn[["centers"]][, 2])

  text(center_x, center_y, glue("all:{venn_intersect_label}"))
  all_centers <- all_venn[["centers"]]

  ## To get a number placed at the edge of each region, I must
  ## find where on the unit circle the lm_center is with respect to the actual
  ## center in radians. If that number is calculated as deg_lm,
  ## then I can take the ~0.9 * lm_diameter * sin(deg_lm) and 0.9 *
  ## lm_diameter * cos(deg_lm) and add it to center_lm to get reasonable
  ## coordinates for putting the lm-only number once I have these coordinates
  ## for each lm/tc/tb, I can average them to get lm/tc and lm/tb
  get_single_edge <- function(name) {
    message("hmm")
  }

  angles <- list()
  radii <- list()
  edges_x <- list()
  edges_y <- list()
  for (i in seq_along(ones)) {
    factor <- 0.98
    single <- ones[i]
    single_name <- names(single)
    single_value <- ones[[i]]
    single_center_x <- all_centers[single_name, "x"]
    single_center_y <- all_centers[single_name, "y"]
    single_rise <- single_center_y - center_y
    single_run <- single_center_x - center_x
    single_angle <- atan2(single_rise, single_run)
    angles[[single_name]] <- single_angle
    single_radius <- all_venn[["diameters"]][[single_name]] / 2.0
    radii[[single_name]] <- single_radius
    single_x_add <- factor * single_radius * cos(single_angle)
    single_y_add <- factor * single_radius * sin(single_angle)
    single_x_edge <- all_centers[single_name, "x"] + single_x_add
    single_y_edge <- all_centers[single_name, "y"] + single_y_add
    edges_x[[single_name]] <- single_x_edge
    edges_y[[single_name]] <- single_y_edge
    text(single_x_edge, single_y_edge,
         glue("{single_name}:{as.character(single_value)}"))
  }

  if (isTRUE(do_doubles)) {
    for (i in seq_along(twos)) {
      factor <- 0.97
      double <- twos[i]
      double_name <- names(double)
      double_value <- twos[[i]]
      name_pair <- strsplit(x = double_name, split = "&")[[1]]
      first_name <- name_pair[[1]]
      second_name <- name_pair[[2]]
      middle_x <- (edges_x[[first_name]] + edges_x[[second_name]]) / 2.0
      middle_y <- (edges_y[[first_name]] + edges_y[[second_name]]) / 2.0
      middle_rise <- middle_y - center_y
      middle_run <- middle_x - center_x
      middle_angle <- atan2(middle_rise, middle_run)
      middle_radius <- (radii[[first_name]] + radii[[second_name]]) / 2.0
      middle_x_add <- factor * middle_radius * cos(middle_angle)
      middle_y_add <- factor * middle_radius * sin(middle_angle)
      middle_x_edge <- center_x + middle_x_add
      middle_y_edge <- center_y + middle_y_add
      text(middle_x_edge, middle_y_edge,
           glue("{double_name}:{as.character(double_value)}"))
    }
  }

  if (isTRUE(do_triples)) {
    for (i in seq_along(threes)) {
      factor <- 0.8
      triple <- threes[i]
      triple_name <- names(triple)
      triple_value <- threes[[i]]
      name_pair <- strsplit(x = triple_name, split = "&")[[1]]
      first_name <- name_pair[[1]]
      third_name <- name_pair[[3]]
      ## this assumes they are given as 1,2,3 where 2 is between 1 and 3 on the circle
      middle_x <- (edges_x[[first_name]] + edges_x[[third_name]]) / 2.0
      middle_y <- (edges_y[[first_name]] + edges_y[[third_name]]) / 2.0
      middle_rise <- middle_y - center_y
      middle_run <- middle_x - center_x
      middle_angle <- atan2(middle_rise, middle_run)
      middle_radius <- (radii[[first_name]] + radii[[third_name]]) / 2.0
      middle_x_add <- factor * middle_radius * cos(middle_angle)
      middle_y_add <- factor * middle_radius * sin(middle_angle)
      middle_x_edge <- center_x + middle_x_add
      middle_y_edge <- center_y + middle_y_add
      text(middle_x_edge, middle_y_edge,
           glue("{triple_name}:{as.character(triple_value)}"))
    }
  }

  if (isTRUE(do_quads)) {
    for (i in seq_along(fours)) {
      factor <- 0.7
      quad <- fours[i]
      quad_name <- names(quad)
      quad_value <- fours[[i]]
      name_pair <- strsplit(x = triple_name, split = "&")[[1]]
      first_name <- name_pair[[1]]
      fourth_name <- name_pair[[4]]
      ## this assumes they are given as 1,2,3,4 where 2,3 is between 1 and 4 on the circle
      middle_x <- (edges_x[[first_name]] + edges_x[[fourth_name]]) / 2.0
      middle_y <- (edges_y[[first_name]] + edges_y[[fourth_name]]) / 2.0
      middle_rise <- middle_y - center_y
      middle_run <- middle_x - center_x
      middle_angle <- atan2(middle_rise, middle_run)
      middle_radius <- (radii[[first_name]] + radii[[fourth_name]]) / 2.0
      middle_x_add <- factor * middle_radius * cos(middle_angle)
      middle_y_add <- factor * middle_radius * sin(middle_angle)
      middle_x_edge <- center_x + middle_x_add
      middle_y_edge <- center_y + middle_y_add
      text(middle_x_edge, middle_y_edge,
           glue("{quad_name}:{as.character(triple_value)}"))
    }
  }
  retlist <- list(
      "venn_data" = all_venn,
      "plot" = grDevices::recordPlot())
  return(retlist)
}

rename_vennerable_intersections <- function(venn, lst) {
  intersects <- venn@IntersectionSets
  list_names <- names(lst)
  for (i in seq(from = 2, to = length(intersects) - 1)) {
    characters <- names(intersects)[i]
    characters <- strsplit(x = characters, split = "")[[1]]
    new_name <- ""
    for (c in seq_along(characters)) {
      char <- characters[c]
      if (char == "1") {
        list_name <- list_names[c]
        new_name <- glue("{new_name}{list_name}, ")
      }
    }
    new_name <- gsub(pattern = ", $", replacement = "", x = new_name)
    names(intersects)[i] <- new_name
  } ## Iterating through every intersection
  names(intersects)[1] <- "none"
  names(intersects)[length(intersects)] <- "all"
  return(intersects)
}

get_vennerable_rows <- function(tables, intersections) {
  ## Skip the 'none' table.
  int_tables <- intersections
  table_names <- names(tables)
  int_names <- names(intersections)
  ## Skip 'none'
  for (t in seq(from = 2, to = length(intersections))) {
    int_name <- int_names[t]
    chosen_table_name <- strsplit(x = int_name, split = ", ")[[1]][1]
    chosen_table <- tables[[chosen_table_name]]
    chosen_rows <- intersections[[t]]
    rows <- chosen_table[chosen_rows, ]
    int_tables[[t]] <- rows
  }
  return(int_tables)
}
