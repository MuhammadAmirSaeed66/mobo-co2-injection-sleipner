# ============================================================
# Download the exact Julia simulator used by this repository
# ============================================================

options(stringsAsFactors = FALSE)

SIMULATOR_REPOSITORY <- paste0(
  "https://github.com/",
  "MuhammadAmirSaeed66/CO2InjectionModeling.jl.git"
)

find_repository_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(current, "mobo-co2-injection-sleipner.Rproj"),
      file.path(current, "README.md")
    )

    if (any(file.exists(markers))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Repository root not found. Open the .Rproj file or run R ",
        "from the repository root."
      )
    }

    current <- parent
  }
}

PROJECT_ROOT <- find_repository_root()

TARGET_DIR <- file.path(
  PROJECT_ROOT,
  "julia",
  "CO2InjectionModeling.jl"
)

COMMIT_FILE <- file.path(
  PROJECT_ROOT,
  "julia",
  "SIMULATOR_COMMIT.txt"
)

if (dir.exists(TARGET_DIR) &&
    length(list.files(TARGET_DIR, all.files = TRUE, no.. = TRUE)) > 0L) {
  message("Julia simulator already exists at:")
  message("  ", TARGET_DIR)
  message("Nothing was downloaded.")
} else {
  git_path <- Sys.which("git")

  if (!nzchar(git_path)) {
    stop(
      "Git was not found on this computer.\n",
      "Install Git or GitHub Desktop, restart RStudio, and run this ",
      "script again."
    )
  }

  dir.create(dirname(TARGET_DIR), recursive = TRUE, showWarnings = FALSE)

  if (dir.exists(TARGET_DIR)) {
    unlink(TARGET_DIR, recursive = TRUE, force = TRUE)
  }

  message("Downloading the Julia simulator from:")
  message("  ", SIMULATOR_REPOSITORY)

  clone_status <- system2(
    git_path,
    c(
      "clone",
      "--depth", "1",
      SIMULATOR_REPOSITORY,
      TARGET_DIR
    )
  )

  if (!identical(clone_status, 0L)) {
    stop(
      "Git could not clone the Julia simulator repository.\n",
      "Check the internet connection and confirm that the GitHub ",
      "repository is public."
    )
  }

  simulator_commit <- system2(
    git_path,
    c(
      "-C", TARGET_DIR,
      "rev-parse", "HEAD"
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  if (length(simulator_commit) < 1L ||
      !grepl("^[0-9a-fA-F]{40}$", simulator_commit[1])) {
    simulator_commit <- "Commit hash could not be determined."
  } else {
    simulator_commit <- simulator_commit[1]
  }

  # Remove the nested Git metadata so the simulator source is committed
  # directly inside the paper-reproducibility repository. The exact source
  # commit is preserved in SIMULATOR_COMMIT.txt.
  unlink(
    file.path(TARGET_DIR, ".git"),
    recursive = TRUE,
    force = TRUE
  )

  writeLines(
    c(
      paste0("Repository: ", SIMULATOR_REPOSITORY),
      paste0("Downloaded commit: ", simulator_commit),
      paste0("Downloaded on: ", format(Sys.time(), tz = "UTC")),
      "",
      "The simulator source is vendored into julia/CO2InjectionModeling.jl.",
      "Its original license and attribution files must be retained."
    ),
    COMMIT_FILE
  )

  message("Julia simulator downloaded successfully.")
}

required_paths <- c(
  file.path(TARGET_DIR, "Project.toml"),
  file.path(TARGET_DIR, "src")
)

missing_paths <- required_paths[
  !file.exists(required_paths) &
    !dir.exists(required_paths)
]

if (length(missing_paths) > 0L) {
  stop(
    "The downloaded repository is missing required Julia project items:\n  ",
    paste(missing_paths, collapse = "\n  ")
  )
}

source_files <- list.files(
  file.path(TARGET_DIR, "src"),
  pattern = "\\.jl$",
  full.names = TRUE
)

if (length(source_files) < 1L) {
  stop(
    "No Julia source file was found under:\n  ",
    file.path(TARGET_DIR, "src")
  )
}

if (!file.exists(file.path(TARGET_DIR, "Manifest.toml"))) {
  warning(
    "Manifest.toml was not found. Project.toml is sufficient for ",
    "Pkg.instantiate(), but committing a Manifest.toml provides stronger ",
    "version reproducibility."
  )
}

message("")
message("Simulator installation check passed.")
message("Project.toml: ", file.path(TARGET_DIR, "Project.toml"))
message("Source directory: ", file.path(TARGET_DIR, "src"))
message("")
message("Next run:")
message('  source("R/install_packages.R")')
message('  source("R/examples/01_one_well.R")')
