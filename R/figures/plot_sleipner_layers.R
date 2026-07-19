# Generate a static Sleipner plume-layer plot from parsed CSV exports.
# Source data are not redistributed here; see data/README.md.

source(file.path("R", "project_setup.R"))

required_packages <- c("sf", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

csv_folder <- Sys.getenv(
  "SLEIPNER_CSV_DIR",
  unset = file.path(
    PROJECT_ROOT,
    "data",
    "sleipner_2019",
    "CSV_exports"
  )
)

file_pattern <- "_parsed\\.csv$"
parsed_files <- list.files(
  csv_folder,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(parsed_files) == 0L) {
  stop(
    "No parsed Sleipner CSV files found in: ", csv_folder,
    "\nSee data/README.md."
  )
}

csv_to_sf <- function(fpath) {
  df <- read.csv(fpath, stringsAsFactors = FALSE)
  required <- c("SEG_ID", "X_EASTING", "Y_NORTHING", "Layer")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(basename(fpath), " is missing: ", paste(missing, collapse = ", "))
  }

  rings <- lapply(split(df, df$SEG_ID), function(d) {
    coords <- as.matrix(d[, c("X_EASTING", "Y_NORTHING")])
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }
    coords
  })

  sf::st_sf(
    Layer = unique(df$Layer)[1],
    geometry = sf::st_sfc(sf::st_multipolygon(list(rings))),
    crs = 25832
  )
}

plume_sf <- do.call(rbind, lapply(parsed_files, csv_to_sf))

p <- ggplot2::ggplot(plume_sf) +
  ggplot2::geom_sf(
    ggplot2::aes(fill = Layer),
    alpha = 0.55,
    color = "black",
    linewidth = 0.25
  ) +
  ggplot2::coord_sf() +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Sleipner CO2 plume boundaries by numerical layer",
    fill = "Layer"
  )

output_file <- file.path(
  RESULT_ROOT,
  "figures",
  "sleipner_plume_layers.png"
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggplot2::ggsave(
  output_file,
  p,
  width = 7,
  height = 6,
  dpi = 350
)

message("Saved: ", output_file)
