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

#' @title Compensation Functions for FlowJo v11
#' @name compensation
#' @keywords internal
#' @importFrom flowCore compensation
NULL

#' Compensation from a custom override, when one exists
#'
#' Handles the named-list override (per-sample) and the single-override
#' (all samples) forms. Returns the validated compensation or NULL when this
#' sample has no override.
#'
#' @param sample_uuid Sample being processed
#' @param custom_compensation Override value passed by the caller
#' @return Validated compensation object or NULL
#' @noRd
comp_custom_for_sample <- function(sample_uuid, custom_compensation) {
    if (is.null(custom_compensation)) {
        return(NULL)
    }

    if (is.list(custom_compensation) && !is.data.frame(custom_compensation)) {
        # Named list - check if this sample has custom comp
        if (sample_uuid %in% names(custom_compensation)) {
            return(validate_compensation(custom_compensation[[sample_uuid]]))
        }
        return(NULL)
    }

    # Single compensation for all samples
    validate_compensation(custom_compensation)
}

#' Compensation for one sample from the configured sources
#'
#' Priority: custom override, then platform matrix, then the workspace
#' data source itself.
#'
#' @param sample_uuid Sample being processed
#' @param custom_comp Validated custom override or NULL
#' @param sample_comp_map Sample-to-platform-compensation map or NULL
#' @param ds Data source entry for the sample
#' @return Compensation object or NULL (with a warning) when none found
#' @noRd
comp_platform_or_legacy <- function(sample_uuid, custom_comp,
    sample_comp_map, ds) {
    if (!is.null(custom_comp)) {
        return(custom_comp)
    }

    # Try compensation from platforms first (sample-specific)
    if (!is.null(sample_comp_map) &&
        !is.null(sample_comp_map[[sample_uuid]])) {
        return(sample_comp_map[[sample_uuid]])
    }

    # Extract from workspace (legacy method)
    comp_matrix <- extract_workspace_compensation(ds)

    if (is.null(comp_matrix)) {
        warning("No compensation found for sample: ", sample_uuid)
    }
    comp_matrix
}

#' Extract Compensation Matrices from FlowJo v11 Workspace
#'
#' Extracts compensation matrices for each sample, with option to override
#' with custom compensation.
#'
#' @param dataSources Data sources from workspace
#' @param sample_uuids Vector of sample UUIDs
#' @param custom_compensation NULL, or compensation object/matrix/data.frame,
#'   or named list of these (names = sample UUIDs)
#' @return Named list of compensation matrices (one per sample)
#' @keywords internal
#' @importFrom flowCore compensation
extract_compensation <- function(dataSources,
    sample_uuids,
    custom_compensation = NULL,
    platforms = NULL,
    schema_version = NULL) {
    # First, build a map of sample UUID to compensation matrix from platforms
    sample_comp_map <- NULL
    if (!is.null(platforms) && !is.null(platforms$spilloverMatrix)) {
        sample_comp_map <- extract_compensation_from_platforms(
            platforms$spilloverMatrix, dataSources, sample_uuids
        )
    }

    comp_list <- list()
    for (sample_uuid in sample_uuids) {
        ds <- dataSources[[sample_uuid]]

        custom_comp <- comp_custom_for_sample(sample_uuid,
            custom_compensation)
        comp_list[[sample_uuid]] <- comp_platform_or_legacy(sample_uuid,
            custom_comp, sample_comp_map, ds)
    }

    comp_list
}


#' Normalize a FlowJo platform matrix to flowCore spillover units
#'
#' FlowJo v11 platform matrices are often scaled so that the diagonal is 100
#' (i.e. they are 100 * the spillover matrix); flowCore::compensation() expects
#' a true spillover matrix with a diagonal of 1, so normalize.
#'
#' @param comp_matrix Numeric matrix with channel dimnames
#' @return flowCore compensation object
#' @noRd
comp_normalize_matrix <- function(comp_matrix) {
    if (all(diag(comp_matrix) > 50)) {
        comp_matrix <- comp_matrix / 100
    }
    flowCore::compensation(comp_matrix)
}

#' Build a compensation object from the compSpec section
#'
#' @param comp_data One entry of the platforms spilloverMatrix list
#' @return flowCore compensation object or NULL when compSpec is absent
#' @noRd
comp_from_comp_spec <- function(comp_data) {
    if (is.null(comp_data$definition$compSpec$CompensationSpec)) {
        return(NULL)
    }
    comp_spec <- comp_data$definition$compSpec$CompensationSpec

    # Get coefficients
    coefficients <- comp_spec$coefficients
    if (is.list(coefficients)) {
        coefficients <- do.call(rbind, coefficients)
    }

    # Get compensated channel names from parameters (e.g.,
    #   "Comp-FITC-A", not "FITC-A")
    # FlowJo 11 stores the compensated names in the parameters list
    comp_names <- vapply(
        comp_spec$parameters, function(p) p$name %||% NA_character_,
        character(1)
    )

    # Create matrix
    comp_matrix <- as.matrix(coefficients)
    rownames(comp_matrix) <- comp_names
    colnames(comp_matrix) <- comp_names
    comp_matrix <- matrix(
        as.numeric(unlist(comp_matrix)),
        nrow      = nrow(comp_matrix),
        ncol      = ncol(comp_matrix),
        dimnames  = dimnames(comp_matrix)
    )
    comp_normalize_matrix(comp_matrix)
}

#' Build a compensation object from the spillover section
#'
#' Fallback used when compSpec did not yield a matrix.
#'
#' @param comp_data One entry of the platforms spilloverMatrix list
#' @return flowCore compensation object or NULL when spillover is absent
#' @noRd
comp_from_spillover <- function(comp_data) {
    if (is.null(comp_data$spillover)) {
        return(NULL)
    }
    spillover <- comp_data$spillover
    values <- spillover$values
    if (is.list(values)) {
        values <- do.call(rbind, values)
    }

    columns <- spillover$columns

    comp_matrix <- as.matrix(values)
    colnames(comp_matrix) <- columns
    comp_matrix <- matrix(
        as.numeric(unlist(comp_matrix)),
        nrow      = nrow(comp_matrix),
        ncol      = ncol(comp_matrix),
        dimnames  = dimnames(comp_matrix)
    )
    comp_normalize_matrix(comp_matrix)
}

#' Map each sample to its compensation via dataSources' parent references
#'
#' When no sample-specific mapping is found but compensations exist, the
#' first compensation is used for all samples.
#'
#' @param comp_map Compensation UUID to compensation map
#' @param dataSources Data sources to map samples to compensation
#' @param sample_uuids Vector of sample UUIDs
#' @return Named list of compensation matrices (one per sample)
#' @noRd
comp_map_samples <- function(comp_map, dataSources, sample_uuids) {
    sample_comp_map <- list()
    for (sample_uuid in sample_uuids) {
        ds <- dataSources[[sample_uuid]]
        if (!is.null(ds) && !is.null(ds$parents) &&
            !is.null(ds$parents$platforms)) {
            # Get the compensation UUID(s) this sample references
            comp_uuids <- ds$parents$platforms
            if (length(comp_uuids) > 0) {
                comp_uuid <- comp_uuids[[1]]
                if (!is.null(comp_map[[comp_uuid]])) {
                    sample_comp_map[[sample_uuid]] <- comp_map[[comp_uuid]]
                }
            }
        }
    }

    # If no sample-specific mappings found but we have compensations, use the
    #   first one for all
    if (length(sample_comp_map) == 0 && length(comp_map) > 0) {
        first_comp <- comp_map[[1]]
        for (sample_uuid in sample_uuids) {
            sample_comp_map[[sample_uuid]] <- first_comp
        }
    }

    sample_comp_map
}

#' Extract Compensation from Platforms Section
#' @param spillover_matrices List of spillover matrices from platforms
#' @param dataSources Data sources to map samples to compensation
#' @param sample_uuids Vector of sample UUIDs
#' @return Named list of compensation matrices (one per sample)
#' @keywords internal
extract_compensation_from_platforms <- function(spillover_matrices,
    dataSources = NULL, sample_uuids = NULL) {
    # Return NULL if no spillover matrices
    if (length(spillover_matrices) == 0) {
        return(NULL)
    }

    # Build a map of compensation UUID to compensation matrix
    comp_map <- list()
    for (comp_uuid in names(spillover_matrices)) {
        comp_data <- spillover_matrices[[comp_uuid]]

        comp_matrix <- comp_from_comp_spec(comp_data)

        # Try spillover section if compSpec didn't work
        if (is.null(comp_matrix)) {
            comp_matrix <- comp_from_spillover(comp_data)
        }

        if (!is.null(comp_matrix)) {
            comp_map[[comp_uuid]] <- comp_matrix
        }
    }

    # If no dataSources provided, return the comp_map
    if (is.null(dataSources) || is.null(sample_uuids)) {
        # If only one compensation, return it for all samples
        if (length(comp_map) == 1) {
            return(comp_map[[1]])
        }
        return(comp_map)
    }

    # Map each sample to its compensation based on dataSources' parent
    #   references
    comp_map_samples(comp_map, dataSources, sample_uuids)
}


#' Extract Compensation from Data Source
#' @keywords internal
extract_workspace_compensation <- function(ds) {
    # FlowJo v11 stores compensation in different places
    comp_data <- NULL

    # Check compensationReference
    if (!is.null(ds$compensationReference)) {
        comp_ref <- ds$compensationReference
        # This would reference a compensation definition elsewhere in workspace
        # For now, we'll look for inline compensation
        }

    # Check definition$compensation
    if (!is.null(ds$definition$compensation)) {
        comp_data <- ds$definition$compensation
        }

    # Check customKeywords for SPILL or SPILLOVER
    if (is.null(comp_data) && !is.null(ds$definition$customKeywords)) {
        kw <- ds$definition$customKeywords
        comp_data <- kw$SPILL %||% kw$SPILLOVER %||% kw$`$SPILLOVER`
        }

    if (is.null(comp_data)) {
        return(NULL)
        }

    # Parse compensation data
    comp_matrix <- parse_compensation_data(comp_data)

    return(comp_matrix)
    }


#' Parse a compensation matrix from FCS-style CSV text
#'
#' Format: "n,channel1,channel2,...,val1,val2,..." (common in FCS files).
#'
#' @param comp_data Character scalar in the CSV compensation format
#' @return flowCore compensation object or NULL when the format is invalid
#' @noRd
parse_compensation_string <- function(comp_data) {
    parts <- strsplit(comp_data, ",")[[1]]
    n <- as.integer(parts[1])

    if (is.na(n) | length(parts) != (n + n * n + 1)) {
        warning("Invalid compensation string format")
        return(NULL)
    }

    channels <- parts[2:(n + 1)]
    values <- as.numeric(parts[(n + 2):length(parts)])

    comp_matrix <- matrix(values, nrow = n, ncol = n, byrow = TRUE)
    rownames(comp_matrix) <- channels
    colnames(comp_matrix) <- channels

    flowCore::compensation(comp_matrix)
}

#' Parse a compensation matrix from the FlowJo v11 JSON format
#'
#' @param comp_data List with matrix and parameters entries
#' @return flowCore compensation object or NULL when the entries are absent
#' @noRd
parse_compensation_list <- function(comp_data) {
    if (is.null(comp_data$matrix) || is.null(comp_data$parameters)) {
        return(NULL)
    }

    matrix_data <- comp_data$matrix
    params <- comp_data$parameters

    if (is.list(matrix_data)) {
        matrix_data <- do.call(rbind, matrix_data)
    }

    comp_matrix <- as.matrix(matrix_data)
    rownames(comp_matrix) <- params
    colnames(comp_matrix) <- params

    flowCore::compensation(comp_matrix)
}

#' Parse Compensation Data
#' @keywords internal
#' @importFrom flowCore compensation
parse_compensation_data <- function(comp_data) {
    if (is.matrix(comp_data)) {
        return(flowCore::compensation(comp_data))
    }

    if (is.data.frame(comp_data)) {
        return(flowCore::compensation(as.matrix(comp_data)))
    }

    if (is.character(comp_data)) {
        # Parse from string format (common in FCS files)
        return(parse_compensation_string(comp_data))
    }

    if (is.list(comp_data)) {
        # FlowJo v11 JSON format
        return(parse_compensation_list(comp_data))
    }

    warning("Unrecognized compensation format")
    return(NULL)
}


#' Validate Compensation Object
#' @keywords internal
#' @importFrom methods is
#' @importFrom flowCore compensation
validate_compensation <- function(comp) {
    if (is(comp, "compensation")) {
        return(comp)
        }

    if (is.matrix(comp)) {
        if (is.null(rownames(comp)) || is.null(colnames(comp))) {
            stop(
                    "Compensation matrix must have row and column names",
                    " (channel names)")
            }
        return(flowCore::compensation(comp))
        }

    if (is.data.frame(comp)) {
        return(flowCore::compensation(as.matrix(comp)))
        }

    stop(
            "Invalid compensation object. Must be compensation, matrix, or",
            " data.frame")
    }


#' Warn about compensation channels that did not map to parameters
#'
#' @param name_mapping Mapping produced by map_param_names()
#' @param comp_names_base Source channel names without the "Comp-" prefix
#' @noRd
comp_warn_unmapped <- function(name_mapping, comp_names_base) {
    for (i in seq_along(comp_names_base)) {
        if (is.null(name_mapping[[comp_names_base[i]]])) {
            warning("Could not map compensation channel '",
                comp_names_base[i], "' to any parameter")
        }
    }
}

#' Map Compensation Channel Names to Cytoframe Parameter Names
#'
#' flowCore's compensation() function sanitizes channel names (e.g., "/" ->
#  "_", spaces -> ".").
#' This function creates a mapping between sanitized compensation names and
#  original cytoframe
#' parameter names to ensure proper compensation application.
#'
#' @param comp_matrix Compensation matrix (may have sanitized names)
#' @param param_names Original parameter names from cytoframe
#' @return Compensation matrix with names matching cytoframe parameters
#' @keywords internal
#' @importFrom flowCore compensation
map_compensation_names <- function(comp_matrix, param_names) {
    if (is.null(comp_matrix)) {
        return(NULL)
    }

    # Convert to matrix if compensation object
    if (methods::is(comp_matrix, "compensation")) {
        comp_matrix <- comp_matrix@spillover
    }

    # Get current compensation channel names
    comp_names <- colnames(comp_matrix)

    if (is.null(comp_names)) {
        warning("Compensation matrix has no column names")
        return(comp_matrix)
    }

    # Strip "Comp-" prefix for matching against cytoframe parameters
    # The compensation matrix now has "Comp-" prefixed names (e.g.,
    #   "Comp-FITC-A")
    # but the cytoframe has original names (e.g., "FITC-A")
    comp_names_base <- sub("^Comp-", "", comp_names)

    # Use unified parameter name mapping
    # Compensation names may have "/" that flowCore converts to "_"
    name_mapping <- map_param_names(
        source_names = comp_names_base,
        target_names = param_names,
        strip_comp_prefix = FALSE,
        case_insensitive = FALSE,
        sanitize_slashes = TRUE
    )

    # Check for unmapped names and warn
    comp_warn_unmapped(name_mapping, comp_names_base)

    # Apply mapping to create new compensation matrix with cytoframe parameter
    #   names
    # (without "Comp-" prefix, so flowCore::compensate can match them)
    mapped_names <- apply_param_mapping(comp_names_base, name_mapping,
        on_no_match = "keep")

    # Create new compensation matrix with mapped names
    mapped_comp <- comp_matrix
    colnames(mapped_comp) <- mapped_names
    rownames(mapped_comp) <- mapped_names

    flowCore::compensation(mapped_comp)
}
