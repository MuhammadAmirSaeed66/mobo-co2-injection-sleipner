# ============================================================
# Portable repository and Julia setup
# ============================================================

options(stringsAsFactors = FALSE)

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
        "Repository root not found. Open the .Rproj file or run R from ",
        "the repository root directory."
      )
    }

    current <- parent
  }
}

find_julia_home <- function() {
  env_home <- Sys.getenv("JULIA_HOME", unset = "")
  path_julia <- Sys.which("julia")

  candidates <- c(
    env_home,
    path.expand("~/.juliaup/bin"),
    if (nzchar(path_julia)) dirname(path_julia) else "",
    "/opt/homebrew/bin",
    "/usr/local/bin",
    Sys.glob(
      "/Applications/Julia-*.app/Contents/Resources/julia/bin"
    )
  )

  candidates <- unique(candidates[nzchar(candidates)])

  for (candidate in candidates) {
    julia_executable <- file.path(candidate, "julia")

    if (file.exists(julia_executable)) {
      return(
        normalizePath(
          candidate,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  stop(
    "Julia could not be found.\n\n",
    "On macOS with juliaup, Julia is normally located at:\n",
    "  ~/.juliaup/bin/julia\n\n",
    "Confirm the file exists by running in R:\n",
    '  file.exists(path.expand("~/.juliaup/bin/julia"))\n\n',
    "If it returns TRUE, restart R and run:\n",
    '  Sys.setenv(JULIA_HOME = path.expand("~/.juliaup/bin"))\n',
    '  source("R/examples/01_one_well.R")\n\n',
    "Otherwise install Julia or set JULIA_HOME to the directory that ",
    "contains the julia executable."
  )
}

PROJECT_ROOT <- find_repository_root()

# The simulator can be supplied in either of two ways:
# 1. Download it into julia/CO2InjectionModeling.jl by running
#      source("R/install_julia_simulator.R")
# 2. Set CO2_JULIA_PROJECT to another complete local Julia project.
JULIA_PROJECT <- Sys.getenv(
  "CO2_JULIA_PROJECT",
  unset = file.path(
    PROJECT_ROOT,
    "julia",
    "CO2InjectionModeling.jl"
  )
)

required_julia_items <- c(
  file.path(JULIA_PROJECT, "Project.toml"),
  file.path(JULIA_PROJECT, "src")
)

missing_julia_items <- required_julia_items[
  !file.exists(required_julia_items) &
    !dir.exists(required_julia_items)
]

if (length(missing_julia_items) > 0L) {
  stop(
    "The Julia simulator is not installed correctly.\n\n",
    "From the repository root, run:\n",
    '  source("R/install_julia_simulator.R")\n\n',
    "Expected Julia project:\n  ",
    JULIA_PROJECT,
    "\n\nMissing:\n  ",
    paste(missing_julia_items, collapse = "\n  ")
  )
}

if (!requireNamespace("JuliaCall", quietly = TRUE)) {
  stop(
    "R package 'JuliaCall' is not installed. Run:\n",
    "  source('R/install_packages.R')"
  )
}

JULIA_HOME <- find_julia_home()
Sys.setenv(JULIA_HOME = JULIA_HOME)

julia_executable <- file.path(JULIA_HOME, "julia")
julia_version <- tryCatch(
  system2(
    julia_executable,
    "--version",
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) conditionMessage(e)
)

message("Julia executable: ", julia_executable)
message("Julia version: ", paste(julia_version, collapse = " "))

JuliaCall::julia_setup(
  JULIA_HOME = JULIA_HOME
)

julia_project_normalized <- normalizePath(
  JULIA_PROJECT,
  winslash = "/",
  mustWork = TRUE
)

JuliaCall::julia_command(
  sprintf(
    'using Pkg; Pkg.activate(raw"%s"); Pkg.instantiate()',
    julia_project_normalized
  )
)
JuliaCall::julia_command("using CO2InjectionModeling")

RESULT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "reported"
)
dir.create(
  RESULT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Repository root: ", PROJECT_ROOT)
message("Julia home: ", JULIA_HOME)
message("Julia project: ", JULIA_PROJECT)
message("Results directory: ", RESULT_ROOT)
