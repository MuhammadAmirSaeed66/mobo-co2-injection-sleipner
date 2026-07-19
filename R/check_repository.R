# Lightweight repository checks that do not run the full experiment.

required_files <- c(
  "R/project_setup.R",
  "R/01_run_mobo_final.R",
  "R/examples/01_one_well.R",
  "R/examples/03_three_wells.R",
  "results/reported/tables/manuscript_multiseed_summary_table.csv",
  "results/reported/sessionInfo.txt"
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0L) {
  stop("Missing repository files:\n", paste("-", missing, collapse = "\n"))
}

scripts <- list.files("R", pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
for (script in scripts) {
  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  if (grepl("/Users/dr.amirsaeed/", text, fixed = TRUE)) {
    warning("Personal path remains in: ", script)
  }
}

message("Basic repository file checks passed.")
