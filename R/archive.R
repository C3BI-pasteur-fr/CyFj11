# Copyright (c) 2026 Institut Pasteur
# Author: Bernd Jagla
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' @title Extract FlowJo v11 Workspace
#' @name archive
#' @keywords internal
NULL

#' Process FlowJo v11 ZIP Archive
#'
#' Extracts and parses the contents of a .flowjo file:
#' - manifest.txt: File listing
#' - analysis-*.json: Workspace data (gates, populations, samples)
#'
#' @param zip_path Path to the .flowjo file
#' @return List containing manifest and parsed JSON data
#' @keywords internal
process_zip_archive <- function(zip_path) {
  # Validate input
  if (!file.exists(zip_path)) {
    stop("Workspace file does not exist: ", zip_path)
  }
  # Create unique temporary directory for extraction
  work_dir <- file.path(tempdir(), paste0("processing_", Sys.getpid(), "_", round(runif(1, 10000, 99999))))
  dir.create(work_dir, recursive = TRUE)
  
  # Ensure cleanup happens even if function errors
  on.exit({
    if(dir.exists(work_dir)) {
      unlink(work_dir, recursive = TRUE)
      if (.pkgenv$verbose) message("Cleaned up temporary directory:", work_dir, "\n") # nocov
    }
  })
  
  if (.pkgenv$verbose) message("Created temporary directory:", work_dir, "\n") # nocov
  
  # Extract all files from the ZIP archive
  unzip(zip_path, exdir = work_dir)
  zip_info <- list.files(work_dir, recursive = TRUE, full.names = TRUE)
  if (.pkgenv$verbose) { # nocov
    message("Archive contains ", length(zip_info), " files:")
    if (length(zip_info)) message(paste(zip_info, collapse = "\n"))
  }
  # Find target files (manifest and JSON)
  manifest_files <- grep("manifest\\.txt$", zip_info, value = TRUE)
  json_files <- grep("\\.json$", zip_info, value = TRUE)
  
  if (.pkgenv$verbose) { # nocov
    message("\nFound", length(manifest_files), "manifest file(s)\n")
  message("Found", length(json_files), "JSON file(s)\n")
  }
  # Initialize results structure
  results <- list(
    manifests = list(),
    json = list()
  )
  
  # Read manifest files (plain text)
  for(manifest_file in manifest_files) {
    if(file.exists(manifest_file)) {
      results$manifests[[basename(manifest_file)]] <- readLines(manifest_file, warn = FALSE)
    }
  }
  
  # Read and parse JSON files
  for(json_file in json_files) {
    if(file.exists(json_file)) {
      results$json[[basename(json_file)]] <- tryCatch({
        jsonlite::fromJSON(json_file, simplifyVector = FALSE, simplifyMatrix = FALSE)  # Parse JSON
      }, error = function(e) {
        # If parsing fails, return error and raw content for debugging
        list(error = e$message, raw_content = readLines(json_file, warn = FALSE))
      })
    }
  }
  
  return(results)
}

#' Read FlowJo v11 Workspace
#'
#' Main wrapper function to read and parse FlowJo v11 workspace files.
#'
#' @param workspace_path Path to the FlowJo workspace file (.flowjo)
#' @return Parsed workspace object containing manifest and JSON data
#' @examples
#' ws_path <- system.file("extdata", "min_test.flowjo", package = "CyFj11")
#' ws <- read_flowjo11_workspace(ws_path)
#' length(ws$samples)  # Number of samples in workspace
#' names(ws$groups)    # Group names
#' @export
read_flowjo11_workspace <- function(workspace_path) {
  # Validate input
  if (!file.exists(workspace_path)) {
    stop("Workspace file does not exist: ", workspace_path)
  }
  
  # Check file extension
  if (!grepl("\\.(fjw|flowjo)$", workspace_path)) {
    warning("Workspace file does not have expected .flowjo or .flowjo extension")
  }
  
  # Process the ZIP archive
  if (.pkgenv$verbose) message("Reading FlowJo v11 workspace:", workspace_path, "\n") # nocov
  results <- process_zip_archive(workspace_path)
  
  # Get the main analysis JSON (find the first analysis JSON file)
  main_json_name <- grep("^analysis-.*\\.json$", names(results$json), value = TRUE)
  if (length(main_json_name) == 0) {
    stop("No analysis JSON file found in workspace")
  }
  
  main_json <- results$json[[main_json_name[1]]]
  
  # Create structured workspace object with properly organized components
  workspace <- list(
    path = normalizePath(workspace_path),
    manifest = results$manifests,
    json = results$json,
    # Extract the key components that the rest of the codebase expects
    groups = main_json$groups,
    dataSources = main_json$dataSources,
    populationDefinitions = main_json$populationDefinitions,
    populations = main_json$populations,
    # Include all top-level JSON fields for completeness
    schemaVersion = main_json$schemaVersion,
    analysisUUID = main_json$analysisUUID,
    uri = main_json$uri,
    reports = main_json$reports,
    compoundParameterSets = main_json$compoundParameterSets,
    compoundPopulations = main_json$compoundPopulations,
    paramsetDefinitions = main_json$paramsetDefinitions,
    platforms = main_json$platforms,
    cytometers = main_json$cytometers,
    analysisRoot = main_json$analysisRoot,
    timestamp = Sys.time()
  )
  
  # Add class attribute for S3 methods
  class(workspace) <- "flowjo11_workspace"
  if (.pkgenv$verbose) { # nocov
    message("Successfully parsed FlowJo v11 workspace\n")
    message("  - Manifest files:", length(workspace$manifest), "\n")
    message("  - JSON files:", length(workspace$json), "\n")
    message("  - Groups:", length(workspace$groups), "\n")
    message("  - DataSources:", length(workspace$dataSources), "\n")
    message("  - PopulationDefinitions:", length(workspace$populationDefinitions), "\n")
    message("  - Populations:", length(workspace$populations), "\n")
    message("  - Reports:", length(workspace$reports), "\n")
    message("  - Platforms:", length(workspace$platforms), "\n")
    message("  - Cytometers:", length(workspace$cytometers), "\n")
  }
  return(workspace)
}
