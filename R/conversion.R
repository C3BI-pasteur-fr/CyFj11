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
#' @importFrom flowWorkspace GatingSet cytoset load_cytoset_from_fcs
#'  gs_pop_add gs_get_pop_paths recompute sampleNames gs_pop_get_gate
#' @importFrom dplyr filter enquo
NULL

#' Convert FlowJo v11 Workspace to GatingSet
#'
#' Parses a FlowJo v11 workspace and creates a GatingSet object with gating
#'  hierarchy,
#' compensation, and transformations. Supports biexponential (Biex)
#'  transformations
#' which are mapped to logicle transforms.
#'
#' @param fj11_workspace List returned by \code{read_flowjo11_workspace()}
#'  containing
#'   the complete workspace structure
#' @param group_name Character or numeric. Name or index of the sample group
#'  to import.
#'   If NULL, will list available groups interactively.
#' @param subset Numeric vector, character vector of FCS filenames, or
#'  expression
#'   to filter samples. See Details.
#' @param execute Logical. Should gates be executed immediately? Default TRUE.
#'   If FALSE, creates GatingSet structure without computing cell counts.
#' @param path Character. Root directory to search for FCS files. Required if
#'   FCS files are not in workspace location.
#' @param cytoset A cytoset object providing alternative data source (instead
#'  of FCS files).
#'   Useful for preprocessed data.
#' @param backend_dir Directory for storing h5 backend data. Default is
#'  tempdir().
#' @param backend Character. Storage backend, either "h5" or "tile". Default
#'  "h5".
#' @param include_gates Logical. Should gates be imported? Default TRUE.
#' @param compensation Compensation object, matrix, data.frame, or named list
#'  of these.
#'   Overrides workspace compensation. Names should match sample GUIDs.
#' @param additional.keys Character vector of FCS keywords to include in
#'  sample GUID
#'   (concatenated with filename). Default is "$TOT".
#' @param additional.sampleID Logical. Include FlowJo sample ID in GUID?
#'  Default FALSE.
#' @param keywords Character vector of keywords to extract as pData.
#' @param keyword.ignore.case Logical. Case-insensitive keyword matching?
#'  Default FALSE.
#' @param channel.ignore.case Logical. Case-insensitive channel matching?
#'  Default FALSE.
#' @param extend_val Numeric. Threshold for extending gate coordinates.
#'  Default 0.
#' @param extend_to Numeric. Value to extend gates to. Default -4000.
#' @param leaf.bool Logical. Compute leaf boolean gates? Default TRUE.
#' @param include_empty_tree Logical. Include samples without gates? Default
#'  FALSE.
#' @param correct_faulty_gate Numerical, zero = don't correct, otherwise use
#'   this value as transformation max value and Try to correct faulty gates
#'   instead of erroring? Default TRUE.
#' @param transform Logical. Apply transformations? Default TRUE.
#' @param use_transformed_coords Logical. Keep gate coordinates in FlowJo
#'  display
#' (transformed) space rather than raw space? Default is the value of
#'  \code{transform}.
#'   Only relevant when \code{transform = TRUE}.
#' @param max_search_depth Integer. Maximum directory depth when searching for
#'  FCS files.
#'   Default 5.
#' @param stop_on_multiple Logical. Stop if multiple FCS files match a sample?
#'   When TRUE, stops on both missing files and multiple matches.
#' When FALSE, only stops on missing files and continues with multiple
#'  matches.
#'   Default TRUE.
#' @param mc.cores Integer. Number of cores for parallel processing. Default 1.
#' @param ... Additional arguments passed to FCS parser.
#'
#' @details
#' The \code{subset} argument can be:
#' \describe{
#'   \item{Numeric vector: }{Sample indices to import}
#'   \item{Character vector: }{FCS filenames to import}
#' \item{Expression: }{Logical expression to filter pData (requires
#'  \code{keywords} argument)}
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
#' @return A named list of \code{GatingSet} objects, one per sample. Each
#'  element
#' contains the full gating hierarchy, gates, compensation, and
#'  transformations for
#' that sample. To obtain a single merged \code{GatingSet} for downstream
#'  analysis,
#'   call \code{flowWorkspace::merge_list_to_gs(gsList)} on the returned list.
#'
#' @examples
#' ws_path <- system.file("extdata", "min_test.flowjo", package = "CyFj11")
#' fcs_path <- system.file("extdata", package = "CyFj11")
#' ws <- read_flowjo11_workspace(ws_path)
#' gsList <- fj11_to_gatingset(ws, group_name = 1, path = fcs_path)
#' length(gsList) # Number of GatingSets created
#'
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter
#'  names
#' when adding populations? Default TRUE. Set to FALSE if compensation has
#'  already
#'   been applied and gate names should match the compensated parameter names.
#' @param sanitize_slashes Logical. Replace "/" with "_" in parameter names?
#' Default TRUE (matches flowCore behavior). Set to FALSE if you want to
#'  preserve
#'   "/" in marker names (e.g., "CD3/CD4" stays as-is).
#' @export
#' @importFrom flowWorkspace GatingSet cytoset load_cytoset_from_fcs
#'  gs_pop_add gs_get_pop_paths recompute
#' @importFrom dplyr filter enquo
#' @importFrom utils menu
fj11_to_gatingset <- function(
    fj11_workspace,
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
    ...
) {
    backend <- match.arg(backend)

    prep <- fj11_prepare_data(fj11_workspace, group_name, subset, path,
        cytoset, stop_on_multiple, backend_dir, backend, keywords, ...)
    sample_uuids <- prep$sample_uuids
    # FlowJo >= 11.3 stores the values of algorithm-derived parameters (event
    # quality filter, FlowSOM, t-SNE, UMAP) as CSV files inside the .flowjo
    # archive; gates on them need these columns in the flowFrames.
    cytoset <- append_derived_parameters(
        prep$cytoset, fj11_workspace$path, sample_uuids
    )

    gsList <- fj11_build_gatingset(
        fj11_workspace, sample_uuids, cytoset, include_gates, compensation,
        transform, execute, channel.ignore.case, extend_val, extend_to,
        correct_faulty_gate, use_transformed_coords, keywords,
        additional.keys, additional.sampleID, keyword.ignore.case,
        strip_comp_prefix
    )
    if (.pkgenv$verbose) {
        fj11_report_summary(gsList, execute)
    }
    return(gsList)
}



#' Prepare data for GatingSet conversion
#'
#' Runs Steps 1-4 of the conversion pipeline: selecting the sample
#' group, filtering samples, resolving FCS file paths, and loading
#' the cytoset. Returns NULL paths when a cytoset is supplied.
#'
#' @param fj11_workspace Workspace list from read_flowjo11_workspace
#' @param group_name NULL, numeric index, or character group name
#' @param subset Sample filter passed to filter_samples
#' @param path Root directory for FCS files (NULL when cytoset given)
#' @param cytoset A pre-built cytoset, or NULL
#' @param stop_on_multiple Whether to stop on multiple path matches
#' @param backend_dir Backend directory for the cytoset
#' @param backend Backend type
#' @param ... Extra arguments passed to load_cytoset_from_fcs
#' @return List with \code{sample_uuids} and \code{cytoset}
#' @keywords internal
fj11_prepare_data <- function(fj11_workspace, group_name, subset, path,
                                cytoset, stop_on_multiple, backend_dir,
                                backend, keywords, ...) {
    groups <- fj11_workspace$groups
    dataSources <- fj11_workspace$dataSources

    # Step 1: Select group ----
    group_sel <- fj11_select_group(groups, group_name)
    selected_group <- group_sel$selected_group

    # Step 2: Get samples in group ----
    sample_uuids <- selected_group$results$dataSources

    # Filter samples based on subset argument
    sample_uuids <- filter_samples(
        sample_uuids, subset, dataSources,
        keywords
    )

    if (.pkgenv$verbose) {
        message(
            "Processing", length(sample_uuids),
            "samples\n\n"
        )
    } # nocov

    # Step 3: Resolve FCS file paths ----
    sample_file_map <- fj11_resolve_paths(
        dataSources, sample_uuids, path, cytoset, stop_on_multiple
    )

    # Step 4: Load data into cytoset ----
    if (is.null(cytoset)) {
        cytoset <- fj11_load_cytoset(
            sample_file_map, sample_uuids, backend_dir, backend, ...
        )
    }

    list(sample_uuids = sample_uuids, cytoset = cytoset)
}


#' Prepare data for GatingSet conversion
#'
#' Runs Steps 1-4 of the conversion pipeline: selecting the sample
#' group, filtering samples, resolving FCS file paths, and loading
#' the cytoset. Returns NULL paths when a cytoset is supplied.
#'
#' @param fj11_workspace Workspace list from read_flowjo11_workspace
#' @param group_name NULL, numeric index, or character group name
#' @param subset Sample filter passed to filter_samples
#' @param path Root directory for FCS files (NULL when cytoset given)
#' @param cytoset A pre-built cytoset, or NULL
#' @param stop_on_multiple Whether to stop on multiple path matches
#' @param backend_dir Backend directory for the cytoset
#' @param backend Backend type
#' @param keywords Character vector of keywords used by expression-based
#'   subset filters (passed to filter_samples)
#' @param ... Extra arguments passed to load_cytoset_from_fcs
#' @return List with \code{sample_uuids} and \code{cytoset}

#' Select a sample group from the workspace
#'
#' Resolves the group index from the \code{group_name} argument
#' (NULL prompts interactively, numeric indexes, character matches by
#' name) and returns the group index and the selected group record.
#'
#' @param groups Groups list from the workspace
#' @param group_name NULL, numeric index, or character name
#' @return List with \code{group_idx}, \code{selected_group_uuid}, and
#'   \code{selected_group}
#' @keywords internal
fj11_select_group <- function(groups, group_name) {
    group_info <- get_group_info(groups)

    if (is.null(group_name)) {
        group_idx <- menu(group_info$name,
            graphics = FALSE,
            title = "Choose which group of samples to import:"
        )
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
        message("Contains", length(selected_group$results$dataSources),
            "samples\n")
    }

    list(
        group_idx = group_idx,
        selected_group_uuid = selected_group_uuid,
        selected_group = selected_group
    )
}

#' Stop when FCS path resolution failed
#'
#' Strict mode stops on any non-FOUND status; permissive mode only on
#' truly missing files (NOT_FOUND / NO_URI), allowing MULTIPLE matches.
#'
#' @param path_resolution The path resolution data frame
#' @param stop_on_multiple Whether to stop on multiple matches
#' @return NULL invisibly; stops on failure
#' @keywords internal
fj11_check_path_resolution <- function(path_resolution, stop_on_multiple) {
    # Only stop on missing files when stop_on_multiple=FALSE
    # When stop_on_multiple=TRUE, stop on both missing and multiple matches
    if (stop_on_multiple) {
        # Strict mode: stop on any failure (missing or multiple)
        failed <- path_resolution[path_resolution$status != "FOUND", ]
        if (nrow(failed) > 0) {
            stop(
                "Could not resolve ", nrow(failed),
                " FCS file(s). See path_resolution for details."
            )
        }
    } else {
        # Permissive mode: only stop on truly missing files, allow
        #   multiple matches
        failed <- path_resolution[
            path_resolution$status == "NOT_FOUND" |
                path_resolution$status == "NO_URI",
        ]
        if (nrow(failed) > 0) {
            stop(
                "Could not resolve ", nrow(failed),
                " FCS file(s). Missing files detected. ",
                "See path_resolution for details."
            )
        }
    }
}

#' Build the sample-to-file map from a path resolution
#'
#' Strict mode maps only FOUND samples; permissive mode also includes
#' MULTIPLE samples using their first resolved path.
#'
#' @param path_resolution The path resolution data frame
#' @param stop_on_multiple Whether to stop on multiple matches
#' @return Named character vector mapping sample UUIDs to file paths
#' @keywords internal
fj11_build_file_map <- function(path_resolution, stop_on_multiple) {
    if (stop_on_multiple) {
        # Strict mode: only include FOUND samples
        return(get_sample_file_map(path_resolution))
    }
    # Permissive mode: include both FOUND and MULTIPLE samples
    # For FOUND samples, use resolved_path directly
    found_rows <- path_resolution$status == "FOUND"
    sample_file_map <- setNames(
        path_resolution$resolved_path[found_rows],
        path_resolution$sample_id[found_rows]
    )

    # For MULTIPLE samples, extract the first path from the
    #   resolved_path field
    multiple_rows <- path_resolution$status == "MULTIPLE"
    if (any(multiple_rows)) {
        # Split the resolved_path by " | " and take the first path
        first_paths <- vapply(
            strsplit(
                path_resolution$resolved_path[multiple_rows],
                " \\| "
            ),
            function(p) p[1], character(1)
        )
        names(first_paths) <- path_resolution$sample_id[multiple_rows]
        sample_file_map <- c(sample_file_map, first_paths)
    }
    sample_file_map
}


#' Resolve FCS file paths for the selected samples
#'
#' Runs resolve_all_fcs_paths, filters to the selected samples, checks
#' for resolution failures, and builds the sample file map. Returns
#' NULL when \code{cytoset} is supplied directly.
#'
#' @param dataSources Data sources from the workspace
#' @param sample_uuids Selected sample UUIDs
#' @param path Root directory for FCS files (or NULL when cytoset given)
#' @param cytoset A pre-built cytoset, or NULL
#' @param stop_on_multiple Whether to stop on multiple matches
#' @return Named character vector mapping sample UUIDs to file paths,
#'   or NULL when a cytoset was supplied
#' @keywords internal
fj11_resolve_paths <- function(dataSources, sample_uuids, path, cytoset,
                                stop_on_multiple) {
    if (!is.null(cytoset)) {
        return(NULL)
    }
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
    path_resolution <- path_resolution[
        path_resolution$sample_id %in% sample_uuids,
    ]

    fj11_check_path_resolution(path_resolution, stop_on_multiple)
    fj11_build_file_map(path_resolution, stop_on_multiple)
}


#' Load FCS files into a cytoset
#'
#' @param sample_file_map Named vector mapping sample UUIDs to paths
#' @param sample_uuids Selected sample UUIDs
#' @param backend_dir Backend directory for the cytoset
#' @param backend Backend type ("h5" or "tile")
#' @param ... Extra arguments passed to load_cytoset_from_fcs
#' @return A cytoset object
#' @keywords internal
fj11_load_cytoset <- function(sample_file_map, sample_uuids, backend_dir,
                                backend, ...) {
    if (.pkgenv$verbose) message("\nLoading FCS files into cytoset...\n")
    fcs_files <- sample_file_map[unlist(sample_uuids)]

    # Create cytoset from FCS files
    flowWorkspace::load_cytoset_from_fcs(
        files = fcs_files,
        transformation = FALSE,
        backend_dir = backend_dir,
        backend = backend,
        ...
    )
}

#' Print the GatingSet creation summary banner
#'
#' @param gsList The created list of GatingSet objects
#' @param execute Whether gating was executed
#' @return NULL invisibly
#' @keywords internal
fj11_report_summary <- function(gsList, execute) {
    message("\n")
    message("========================================\n")
    message("  GatingSet List Created Successfully\n")
    message("========================================\n")
    message("Samples:     ", length(gsList), "\n")

    # Safely get population count
    if (!is.null(gsList) && length(gsList) > 0 && !is.null(gsList[[1]])) {
        message(
            "Populations: ",
            length(flowWorkspace::gs_get_pop_paths(gsList[[1]])), "\n"
        )
    } else {
        message("Populations: 0 (GatingSet list is empty)\n")
    }

    message("Execute:     ", execute, "\n")
    message("======================================\n\n")

    # Debug information before returning
    message("DEBUG: About to return GatingSet\n")
    message("DEBUG: gsList is null:", is.null(gsList), "\n")
}


#' Extract workspace content and build the GatingSet
#'
#' Covers the extraction pipeline: building per-sample gating trees
#' (Step 5), extracting transformations (Step 8), gates (Step 6) and
#' compensation (Step 7), creating the GatingSet list (Step 9), and
#' optionally executing gating (Step 10).
#'
#' @param populations Populations list from the workspace
#' @param populationDefinitions Population definitions list
#' @param sample_uuids Selected sample UUIDs
#' @param cytoset The loaded cytoset
#' @param dataSources Data sources from the workspace
#' @param platforms Platforms list from the workspace
#' @param include_gates Whether to extract and apply gates
#' @param compensation Custom compensation matrix or NULL
#' @param transform Whether to apply transformations
#' @param execute Whether to execute gating
#' @param channel.ignore.case Passed to extract_all_gates
#' @param extend_val Passed to extract_all_gates
#' @param extend_to Passed to extract_all_gates
#' @param correct_faulty_gate Passed to extract_all_gates
#' @param use_transformed_coords Passed to extract_all_gates
#' @param keywords Keywords to attach (passed to
#'   create_gatingset_from_cytoset)
#' @param additional.keys Passed to create_gatingset_from_cytoset
#' @param additional.sampleID Passed to
#'   create_gatingset_from_cytoset
#' @param keyword.ignore.case Passed to
#'   create_gatingset_from_cytoset
#' @param strip_comp_prefix Passed to create_gatingset_from_cytoset
#' @return The created list of GatingSet objects
#' @keywords internal
fj11_build_gatingset <- function(fj11_workspace, sample_uuids, cytoset,
                                include_gates, compensation, transform,
                                execute, channel.ignore.case, extend_val,
                                extend_to, correct_faulty_gate,
                                use_transformed_coords, keywords,
                                additional.keys, additional.sampleID,
                                keyword.ignore.case, strip_comp_prefix) {
    # Extract workspace components
    groups <- fj11_workspace$groups
    dataSources <- fj11_workspace$dataSources
    populationDefinitions <- fj11_workspace$populationDefinitions
    populations <- fj11_workspace$populations
    platforms <- fj11_workspace$platforms
    schema_version <- fj11_schema_version(fj11_workspace)

    # Step 5: Build gating hierarchy ----
    if (.pkgenv$verbose) message("\nBuilding gating hierarchy...\n")

    gating_trees <- fj11_build_trees(
        populations, populationDefinitions, sample_uuids, schema_version
    )

    # Steps 6-8: extract transformations, gates, compensation
    comps <- fj11_extract_components(
        fj11_workspace, sample_uuids, include_gates, compensation,
        channel.ignore.case, extend_val, extend_to,
        correct_faulty_gate, use_transformed_coords
    )
    # Step 9: Create per-sample GatingSet list ----
    if (.pkgenv$verbose) message("\nCreating GatingSet list...\n")
    gsList <- fj11_create_gatingsets(
        cytoset, gating_trees, include_gates, comps$gates_list,
        comps$comp_list, transform, comps$trans_list, sample_uuids,
        dataSources, keywords, additional.keys,
        additional.sampleID, keyword.ignore.case,
        strip_comp_prefix, schema_version
    )

    fj11_execute_gating(gsList, execute, include_gates)
    gsList
}


#' Extract transformations, gates, and compensation
#'
#' Runs Steps 6-8 of the conversion pipeline: extracting
#' transformations, gates (when \code{include_gates}), and
#' compensation matrices.
#'
#' @param fj11_workspace Workspace list from read_flowjo11_workspace
#' @param sample_uuids Selected sample UUIDs
#' @param include_gates Whether gates should be extracted
#' @param compensation Custom compensation matrix or NULL
#' @param channel.ignore.case Case-insensitive channel matching
#' @param extend_val Threshold for extending gate coordinates
#' @param extend_to Value to extend gates to
#' @param correct_faulty_gate Correct faulty gates with this max value
#' @param use_transformed_coords Keep gate coordinates in transformed
#'   space
#' @return List with \code{trans_list}, \code{gates_list}, and
#'   \code{comp_list}
#' @keywords internal
fj11_extract_components <- function(fj11_workspace, sample_uuids,
                                    include_gates, compensation,
                                    channel.ignore.case, extend_val,
                                    extend_to, correct_faulty_gate,
                                    use_transformed_coords) {
    populationDefinitions <- fj11_workspace$populationDefinitions
    dataSources <- fj11_workspace$dataSources
    platforms <- fj11_workspace$platforms
    schema_version <- fj11_schema_version(fj11_workspace)

    # Step 8: Extract transformations ----
    if (.pkgenv$verbose) message("\nExtracting transformations...\n")
    trans_list <- extract_transformations(
        populationDefinitions = populationDefinitions,
        sample_uuids = sample_uuids,
        schema_version = schema_version
    )

    # Step 6: Extract gates ----
    gates_list <- if (include_gates) {
        if (.pkgenv$verbose) message("\nExtracting gates...\n")
        fj11_extract_gates(
            populationDefinitions, sample_uuids,
            channel.ignore.case, extend_val, extend_to,
            correct_faulty_gate, use_transformed_coords, schema_version
        )
    }

    # Step 7: Extract compensation ----
    if (.pkgenv$verbose) {
        message("\nExtracting compensation matrices...\n")
    }
    comp_list <- extract_compensation(
        dataSources = dataSources,
        sample_uuids = sample_uuids,
        custom_compensation = compensation,
        platforms = platforms,
        schema_version = schema_version
    )

    list(
        trans_list = trans_list,
        gates_list = gates_list,
        comp_list = comp_list
    )
}

#' Extract gates for the selected samples
#'
#' Thin wrapper around extract_all_gates forwarding the
#' extraction-tuning arguments.
#'
#' @param populationDefinitions Population definitions list
#' @param sample_uuids Selected sample UUIDs
#' @param channel.ignore.case Case-insensitive channel matching
#' @param extend_val Threshold for extending gate coordinates
#' @param extend_to Value to extend gates to
#' @param correct_faulty_gate Correct faulty gates with this max value
#' @param use_transformed_coords Keep gate coordinates in transformed
#'   space
#' @return Gates list from extract_all_gates
#' @keywords internal
fj11_extract_gates <- function(populationDefinitions, sample_uuids,
                                channel.ignore.case, extend_val,
                                extend_to, correct_faulty_gate,
                                use_transformed_coords,
                                schema_version = NULL) {
    extract_all_gates(
        populationDefinitions = populationDefinitions,
        sample_uuids = sample_uuids,
        channel.ignore.case = channel.ignore.case,
        extend_val = extend_val,
        extend_to = extend_to,
        correct_faulty_gate = correct_faulty_gate,
        use_transformed_coords = use_transformed_coords,
        schema_version = schema_version
    )
}

#' Execute gating on the GatingSet list
#'
#' Runs recompute on each GatingSet when executing, or reports the
#' skip when not.
#'
#' @param gsList The created list of GatingSet objects
#' @param execute Whether to execute gating
#' @param include_gates Whether gates were applied
#' @return NULL invisibly
#' @keywords internal
fj11_execute_gating <- function(gsList, execute, include_gates) {
    if (execute && include_gates) {
        if (.pkgenv$verbose) message("\nExecuting gates...\n")
        lapply(gsList, flowWorkspace::recompute)
        if (.pkgenv$verbose) message("Gating complete\n")
    } else {
        if (.pkgenv$verbose) {
            message(
                "Gating not executed ",
                "(set execute=TRUE to compute cell counts)\n"
            )
        } # nocov
    }
    invisible(NULL)
}


#' Build per-sample gating trees
#'
#' Finds the root population and builds a gating tree for each
#' selected sample.
#'
#' @param populations Populations list from the workspace
#' @param populationDefinitions Population definitions list
#' @param sample_uuids Selected sample UUIDs
#' @return Named list of gating trees from build_gating_tree
#' @keywords internal
fj11_build_trees <- function(populations, populationDefinitions,
                            sample_uuids, schema_version = NULL) {
    root_pop_uuid <- find_root_population(
        populations, populationDefinitions, sample_uuids[1]
    )
    lapply(sample_uuids, function(sample_uuid) {
        build_gating_tree(
            sample_uuid = sample_uuid,
            populations = populations,
            populationDefinitions = populationDefinitions,
            root_uuid = root_pop_uuid,
            schema_version = schema_version
        )
    })
}
#' Create the per-sample GatingSet list
#'
#' Thin wrapper around create_gatingset_from_cytoset.
#'
#' @param cytoset The loaded cytoset
#' @param gating_trees Per-sample gating trees
#' @param include_gates Whether gates were extracted
#' @param gates_list Extracted gates (or NULL)
#' @param comp_list Extracted compensations
#' @param transform Whether transformations apply
#' @param trans_list Extracted transformations (or NULL)
#' @param sample_uuids Selected sample UUIDs
#' @param dataSources Data sources from the workspace
#' @param keywords Keywords to attach as pData
#' @param additional.keys Keys for the sample GUID
#' @param additional.sampleID Include FlowJo sample ID in GUID
#' @param keyword.ignore.case Case-insensitive keyword matching
#' @param strip_comp_prefix Strip compensation prefix from parameters
#' @return The created list of GatingSet objects
#' @keywords internal
fj11_create_gatingsets <- function(cytoset, gating_trees,
                                    include_gates, gates_list,
                                    comp_list, transform, trans_list,
                                    sample_uuids, dataSources, keywords,
                                    additional.keys,
                                    additional.sampleID,
                                    keyword.ignore.case,
                                    strip_comp_prefix,
                                    schema_version = NULL) {
    create_gatingset_from_cytoset(
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
        strip_comp_prefix = strip_comp_prefix,
        schema_version = schema_version
    )
}
