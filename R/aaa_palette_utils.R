# Helper utilities for building distinct annotation palettes ----

.annotation_custom_palettes <- function() {
  list(
    "Okabe Ito" = c(
      "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#999999"
    )
  )
}

.available_hcl_palettes <- function() {
  colorspace::hcl_palettes(type = "qualitative")$palette
}

annotation_palette_values <- function(n, palette_name = NULL) {
  if (n <= 0) {
    return(character(0))
  }

  custom <- .annotation_custom_palettes()
  if (!is.null(palette_name) && palette_name %in% names(custom)) {
    pal <- custom[[palette_name]]
    return(rep_len(pal, n))
  }

  palette_to_use <- palette_name
  if (is.null(palette_to_use) || !(palette_to_use %in% .available_hcl_palettes())) {
    palette_to_use <- "Dark 3"
  }

  colorspace::qualitative_hcl(
    n = n,
    palette = palette_to_use,
    c = 90,
    l = 60
  )
}

build_annotation_colors <- function(levels, palette_name = NULL) {
  levels <- unique(as.character(levels))
  if (length(levels) == 0) {
    return(setNames(character(0), character(0)))
  }

  stats::setNames(
    annotation_palette_values(length(levels), palette_name),
    levels
  )
}
