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

#' @title Convert FlowJo v11 Workspace to GatingSet
#' @name conversion
#' @keywords internal
#' @importFrom flowWorkspace GatingSet cytoset load_cytoset_from_fcs gs_pop_add gs_get_pop_paths recompute sampleNames gs_pop_get_gate
#' @importFrom dplyr filter enquo
NULL

#' Convert FlowJo v11 Workspace to GatingSet
#'
#' Parses a FlowJo v11 workspace and creates a GatingSet object with gating hierarchy,
#' compensation, and transformations. Supports biexponential (Biex) transformations
#' which are mapped to logicle transforms.
#'
#' @param fj11_workspace List returned by \code{read_flowjo11_workspace()} containing
#'   the complete workspace structure
#' @param group_name Character or numeric. Name or index of the sample group to import.
#'   If NULL, will list available groups interactively.
#' @param subset Numeric vector, character vector of FCS filenames, or expression
#'   to filter samples. See Details.
#' @param execute Logical. Should gates be executed immediately? Default TRUE.
#'   If FALSE, creates GatingSet structure without computing cell counts.
#' @param path Character. Root directory to search for FCS files. Required if
#'   FCS files are not in workspace location.
#' @param cytoset A cytoset object providing alternative data source (instead of FCS files).
#'   Useful for preprocessed data.
#' @param backend_dir Directory for storing h5 backend data. Default is tempdir().
#' @param backend Character. Storage backend, either "h5" or "tile". Default "h5".
#' @param include_gates Logical. Should gates be imported? Default TRUE.
#' @param compensation Compensation object, matrix, data.frame, or named list of these.
#'   Overrides workspace compensation. Names should match sample GUIDs.
#' @param additional.keys Character vector of FCS keywords to include in sample GUID
#'   (concatenated with filename). Default is "$TOT".
#' @param additional.sampleID Logical. Include FlowJo sample ID in GUID? Default FALSE.
#' @param keywords Character vector of keywords to extract as pData.
#' @param keyword.ignore.case Logical. Case-insensitive keyword matching? Default FALSE.
#' @param channel.ignore.case Logical. Case-insensitive channel matching? Default FALSE.
#' @param extend_val Numeric. Threshold for extending gate coordinates. Default 0.
#' @param extend_to Numeric. Value to extend gates to. Default -4000.
#' @param leaf.bool Logical. Compute leaf boolean gates? Default TRUE.
#' @param include_empty_tree Logical. Include samples without gates? Default FALSE.
#' @param correct_faulty_gate Numerical, zero = don't correct, otherwise use this value as transformation max value and Try to correct faulty gates instead of erroring? Default TRUE.
#' @param transform Logical. Apply transformations? Default TRUE.
#' @param use_transformed_coords Logical. Keep gate coordinates in FlowJo display
#'   (transformed) space rather than raw space? Default is the value of \code{transform}.
#'   Only relevant when \code{transform = TRUE}.
#' @param max_search_depth Integer. Maximum directory depth when searching for FCS files.
#'   Default 5.
#' @param stop_on_multiple Logical. Stop if multiple FCS files match a sample?
#'   When TRUE, stops on both missing files and multiple matches.
#'   When FALSE, only stops on missing files and continues with multiple matches.
#'   Default TRUE.
#' @param mc.cores Integer. Number of cores for parallel processing. Default 1.
#' @param ... Additional arguments passed to FCS parser.
#'
#' @details
#' The \code{subset} argument can be:
#' \describe{
#'   \item{Numeric vector: }{Sample indices to import}
#'   \item{Character vector: }{FCS filenames to import}
#'   \item{Expression: }{Logical expression to filter pData (requires \code{keywords} argument)}
#' }
#'
#' FlowJo v11 stores data differently than v10:
#' \describe{
#'   \item{populationDefinitions: }{Master gate definitions}
#'   \item{desyncTable: }{Per-sample gate variations (tailored gates)}
#'   \item{groups: }{Sample groupings}
#'   \item{dataSources: }{FCS file references}
#' }
#'
#' @return A named list of \code{GatingSet} objects, one per sample. Each element
#'   contains the full gating hierarchy, gates, compensation, and transformations for
#'   that sample. To obtain a single merged \code{GatingSet} for downstream analysis,
#'   call \code{flowWorkspace::merge_list_to_gs(gsList)} on the returned list.
#'
#' @examples
#' ws_path  <- system.file("extdata", "min_test.flowjo", package = "CyFj11")
#' fcs_path <- system.file("extdata", package = "CyFj11")
#' ws       <- read_flowjo11_workspace(ws_path)
#' gsList   <- fj11_to_gatingset(ws, group_name = 1, path = fcs_path)
#' length(gsList)  # Number of GatingSets created
#'
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter names
#'   when adding populations? Default TRUE. Set to FALSE if compensation has already
#'   been applied and gate names should match the compensated parameter names.
#' @param sanitize_slashes Logical. Replace "/" with "_" in parameter names?
#'   Default TRUE (matches flowCore behavior). Set to FALSE if you want to preserve
#'   "/" in marker names (e.g., "CD3/CD4" stays as-is).
#' @export
#' @importFrom flowWorkspace GatingSet cytoset load_cytoset_from_fcs gs_pop_add gs_get_pop_paths recompute
#' @importFrom dplyr filter enquo
#' @importFrom utils menu
fj11_to_gatingset <- function(fj11_workspace,
                              group_name = NULL,
                              subset = list(),
                              execute = TRUE,
                              path = NULL,
                              cytoset = NULL,
                              backend_dir = tempdir(),
                              backend = c("h5", "tile"),
                              include_gates = TRUE,
                              compensation = NULL,
                              additional.keys = "$TOT",
                              additional.sampleID = FALSE,
                              keywords = character(),
                              keyword.ignore.case = FALSE,
                              channel.ignore.case = FALSE,
                              extend_val = 0,
                              extend_to = -4000,
                              leaf.bool = TRUE,
                              include_empty_tree = FALSE,
                              correct_faulty_gate = 0,
                              transform = TRUE,
                              use_transformed_coords = transform,
                              max_search_depth = 5,
                              stop_on_multiple = TRUE,
                              mc.cores = 1,
                              strip_comp_prefix = TRUE,
                              sanitize_slashes = TRUE,
                              ...) {
  backend <- match.arg(backend)
  # Extract workspace components
  # The workspace already has the components extracted at the top level
  groups <- fj11_workspace$groups
  dataSources <- fj11_workspace$dataSources
  populationDefinitions <- fj11_workspace$populationDefinitions
  populations <- fj11_workspace$populations
  platforms <- fj11_workspace$platforms
  
  # Step 1: Select group ----
  group_info <- get_group_info(groups)
  
  if (is.null(group_name)) {
    group_idx <- menu(group_info$name, graphics = FALSE,
                      title = "Choose which group of samples to import:")
    if (group_idx == 0) stop("No group selected")
  } else if (is.numeric(group_name)) {
    if (group_name > nrow(group_info)) {
      stop("Invalid group index: ", group_name)
    }
    group_idx <- group_name
  } else if (is.character(group_name)) {
    group_idx <- which(group_info$name == group_name)
    if (length(group_idx) == 0) {
      stop("Group not found: ", group_name)
    }
  }
  
  selected_group_uuid <- group_info$uuid[group_idx]
  selected_group <- groups[[selected_group_uuid]]
  
  if (.pkgenv$verbose) {
    message("Selected group:", group_info$name[group_idx], "\n")
    message("Contains", length(selected_group$results$dataSources), "samples\n")
  }
  # Step 2: Get samples in group ----
  sample_uuids <- selected_group$results$dataSources
  
  # Filter samples based on subset argument
  sample_uuids <- filter_samples(sample_uuids, subset, dataSources, keywords)
  
  if (.pkgenv$verbose) message("Processing", length(sample_uuids), "samples\n\n")
  
  # Step 3: Resolve FCS file paths ----
  if (is.null(cytoset)) {
    if (is.null(path)) {
      stop("Either 'path' or 'cytoset' must be provided")
    }
    
    if (.pkgenv$verbose) message("Resolving FCS file paths...\n")
    path_resolution <- resolve_all_fcs_paths(
      dataSources = dataSources,
      root_dir = path,
      stop_on_multiple = stop_on_multiple,
      stop_on_missing = TRUE
    )
    
    # Filter to selected samples
    path_resolution <- path_resolution[path_resolution$sample_id %in% sample_uuids, ]
    
    # Check for resolution failures
    # Only stop on missing files when stop_on_multiple=FALSE
    # When stop_on_multiple=TRUE, stop on both missing and multiple matches
    if (stop_on_multiple) {
      # Strict mode: stop on any failure (missing or multiple)
      failed <- path_resolution[path_resolution$status != "FOUND", ]
      if (nrow(failed) > 0) {
        stop("Could not resolve ", nrow(failed), " FCS file(s). See path_resolution for details.")
      }
    } else {
      # Permissive mode: only stop on truly missing files, allow multiple matches
      failed <- path_resolution[path_resolution$status == "NOT_FOUND" | path_resolution$status == "NO_URI", ]
      if (nrow(failed) > 0) {
        stop("Could not resolve ", nrow(failed), " FCS file(s). Missing files detected. See path_resolution for details.")
      }
    }
    
    # Create sample file map
    if (stop_on_multiple) {
      # Strict mode: only include FOUND samples
      sample_file_map <- get_sample_file_map(path_resolution)
    } else {
      # Permissive mode: include both FOUND and MULTIPLE samples
      # For FOUND samples, use resolved_path directly
      found_rows <- path_resolution$status == "FOUND"
      sample_file_map <- setNames(path_resolution$resolved_path[found_rows], path_resolution$sample_id[found_rows])
      
      # For MULTIPLE samples, extract the first path from the resolved_path field
      multiple_rows <- path_resolution$status == "MULTIPLE"
      if (any(multiple_rows)) {
        # Split the resolved_path by " | " and take the first path
        first_paths <- vapply(
          strsplit(path_resolution$resolved_path[multiple_rows], " \\| "),
          function(p) p[1], character(1)
        )
        names(first_paths) <- path_resolution$sample_id[multiple_rows]
        sample_file_map <- c(sample_file_map, first_paths)
      }
    }
  }
  
  # Step 4: Load data into cytoset ----
  if (is.null(cytoset)) {
    if (.pkgenv$verbose) message("\nLoading FCS files into cytoset...\n")
    fcs_files <- sample_file_map[unlist(sample_uuids)]
    
    # Create cytoset from FCS files
    cytoset <- flowWorkspace::load_cytoset_from_fcs(
      files = fcs_files,
      transformation = FALSE,
      backend_dir = backend_dir,
      backend = backend,
      ...
    )
  }
  
  # Step 5: Build gating hierarchy ----
  if (.pkgenv$verbose) message("\nBuilding gating hierarchy...\n")
  
  # Get root population (usually the ungated data)
  root_pop_uuid <- find_root_population(populations, populationDefinitions, sample_uuids[1])
  # browser()
  # Build hierarchy tree for each sample
  # browser()
  gating_trees <- lapply(sample_uuids, function(sample_uuid) {
    build_gating_tree(
      sample_uuid = sample_uuid,
      populations = populations,
      populationDefinitions = populationDefinitions,
      root_uuid = root_pop_uuid
    )
  })
  # Step 8: Extract transformations ----
  if (.pkgenv$verbose) message("\nExtracting transformations...\n")
  trans_list <- extract_transformations(
    populationDefinitions = populationDefinitions,
    sample_uuids = sample_uuids
  )
  
  # Step 6: Extract gates ----
  if (include_gates) {
    if (.pkgenv$verbose) message("\nExtracting gates...\n")

    # save(file = "extract_all_gates.Rdata", list = ls())
    # load("extract_all_gates.Rdata")
    gates_list <- extract_all_gates(
      populationDefinitions = populationDefinitions,
      sample_uuids = sample_uuids,
      channel.ignore.case = channel.ignore.case,
      extend_val = extend_val,
      extend_to = extend_to,
      correct_faulty_gate = correct_faulty_gate,
      use_transformed_coords = use_transformed_coords
    )
  }
  # browser()
  # Step 7: Extract compensation ----
  if (.pkgenv$verbose) message("\nExtracting compensation matrices...\n")
  comp_list <- extract_compensation(
    dataSources = dataSources,
    sample_uuids = sample_uuids,
    custom_compensation = compensation,
    platforms = platforms
  )
  
  # browser()
  # Step 9: Create per-sample GatingSet list ----
  if (.pkgenv$verbose) message("\nCreating GatingSet list...\n")
  gsList <- create_gatingset_from_cytoset(
    cytoset = cytoset,
    gating_trees = gating_trees,
    gates = if (include_gates) gates_list else NULL,
    compensations = comp_list,
    transformations = if (transform) trans_list else NULL,
    sample_uuids = sample_uuids,
    dataSources = dataSources,
    keywords = keywords,
    additional.keys = additional.keys,
    additional.sampleID = additional.sampleID,
    keyword.ignore.case = keyword.ignore.case,
    strip_comp_prefix = strip_comp_prefix
  )

  # browser()
  # Step 10: Execute gating ----
  if (execute && include_gates) {
    if (.pkgenv$verbose) message("\nExecuting gates...\n")

    lapply(gsList, flowWorkspace::recompute)

    if (.pkgenv$verbose) message("Gating complete\n")
  } else {
    if (.pkgenv$verbose) message("Gating not executed (set execute=TRUE to compute cell counts)\n")
  }

  # # Step 11: Compute boolean gates ----
  if (.pkgenv$verbose) {
    message("\n")
    message("========================================\n")
    message("  GatingSet List Created Successfully\n")
    message("========================================\n")
    message("Samples:     ", length(gsList), "\n")

    # Safely get population count
    if (!is.null(gsList) && length(gsList) > 0 && !is.null(gsList[[1]])) {
      message("Populations: ", length(flowWorkspace::gs_get_pop_paths(gsList[[1]])), "\n")
    } else {
      message("Populations: 0 (GatingSet list is empty)\n")
    }

    message("Execute:     ", execute, "\n")
    message("======================================\n\n")

    # Debug information before returning
    message("DEBUG: About to return GatingSet\n")
    message("DEBUG: gsList is null:", is.null(gsList), "\n")
  }
  return(gsList)
}




