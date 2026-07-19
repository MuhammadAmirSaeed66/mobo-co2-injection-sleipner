# Install the R packages used by the main analysis and repository examples.

required_packages <- c(
  "dplyr",
  "tidyr",
  "lhs",
  "nnet",
  "JuliaCall",
  "ggplot2",
  "ranger",
  "sf"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) == 0L) {
  message("All required R packages are already installed.")
} else {
  message("Installing: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}
