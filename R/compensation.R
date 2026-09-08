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
                                 platforms = NULL) {
  
  comp_list <- list()

  # First, build a map of sample UUID to compensation matrix from platforms
  sample_comp_map <- NULL
  if (!is.null(platforms) && !is.null(platforms$spilloverMatrix)) {
    sample_comp_map <- extract_compensation_from_platforms(platforms$spilloverMatrix, dataSources, sample_uuids)
  }
  
  for (sample_uuid in sample_uuids) {
    ds <- dataSources[[sample_uuid]]
    
    # Check for custom compensation first
    if (!is.null(custom_compensation)) {
      if (is.list(custom_compensation) && !is.data.frame(custom_compensation)) {
        # Named list - check if this sample has custom comp
        if (sample_uuid %in% names(custom_compensation)) {
          comp_list[[sample_uuid]] <- validate_compensation(
            custom_compensation[[sample_uuid]]
          )
          next
        }
      } else {
        # Single compensation for all samples
        comp_list[[sample_uuid]] <- validate_compensation(custom_compensation)
        next
      }
    }
    
    # Try compensation from platforms first (sample-specific)
    if (!is.null(sample_comp_map) && !is.null(sample_comp_map[[sample_uuid]])) {
      comp_list[[sample_uuid]] <- sample_comp_map[[sample_uuid]]
      next
    }
    
    # Extract from workspace (legacy method)
    comp_matrix <- extract_workspace_compensation(ds)
    
    if (!is.null(comp_matrix)) {
      comp_list[[sample_uuid]] <- comp_matrix
    } else {
      warning("No compensation found for sample: ", sample_uuid)
    }
  }
  
  return(comp_list)
}


#' Extract Compensation from Platforms Section
#' @param spillover_matrices List of spillover matrices from platforms
#' @param dataSources Data sources to map samples to compensation
#' @param sample_uuids Vector of sample UUIDs
#' @return Named list of compensation matrices (one per sample)
#' @keywords internal
extract_compensation_from_platforms <- function(spillover_matrices, dataSources = NULL, sample_uuids = NULL) {
  # Return NULL if no spillover matrices
  if (length(spillover_matrices) == 0) {
    return(NULL)
  }
  
  # Build a map of compensation UUID to compensation matrix
  comp_map <- list()
  for (comp_uuid in names(spillover_matrices)) {
    comp_data <- spillover_matrices[[comp_uuid]]
    comp_matrix <- NULL
    
    # Extract from compSpec
    if (!is.null(comp_data$definition$compSpec$CompensationSpec)) {
      comp_spec <- comp_data$definition$compSpec$CompensationSpec

      # Get coefficients
      coefficients <- comp_spec$coefficients
      if (is.list(coefficients)) {
        coefficients <- do.call(rbind, coefficients)
      }

      # Get compensated channel names from parameters (e.g., "Comp-FITC-A", not "FITC-A")
      # FlowJo 11 stores the compensated names in the parameters list
      comp_names <- vapply(comp_spec$parameters, function(p) p$name %||% NA_character_,
                           character(1))

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
      # FlowJo v11 platform matrices are often scaled so that the diagonal is 100
      # (i.e. they are 100 * the spillover matrix). flowCore::compensation()
      # expects a true spillover matrix with a diagonal of 1, so normalize.
      if (all(diag(comp_matrix) > 50)) {
        comp_matrix <- comp_matrix / 100
      }
      comp_matrix <- flowCore::compensation(comp_matrix)
    }

    # Try spillover section if compSpec didn't work
    if (is.null(comp_matrix) && !is.null(comp_data$spillover)) {
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
      # Normalize FlowJo's 100-diagonal platform matrices to true spillover units.
      if (all(diag(comp_matrix) > 50)) {
        comp_matrix <- comp_matrix / 100
      }
      comp_matrix <- flowCore::compensation(comp_matrix)
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
  
  # Map each sample to its compensation based on dataSources' parent references
  sample_comp_map <- list()
  for (sample_uuid in sample_uuids) {
    ds <- dataSources[[sample_uuid]]
    if (!is.null(ds) && !is.null(ds$parents) && !is.null(ds$parents$platforms)) {
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
  
  # If no sample-specific mappings found but we have compensations, use the first one for all
  if (length(sample_comp_map) == 0 && length(comp_map) > 0) {
    first_comp <- comp_map[[1]]
    for (sample_uuid in sample_uuids) {
      sample_comp_map[[sample_uuid]] <- first_comp
    }
  }
  
  sample_comp_map
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
    # Format: "n,channel1,channel2,...,val1,val2,..."
    parts <- strsplit(comp_data, ",")[[1]]
    n <- as.integer(parts[1])
    
    if (is.na(n) | length(parts) != (n + n*n + 1)) {
      warning("Invalid compensation string format")
      return(NULL)
    }
    
    channels <- parts[2:(n+1)]
    values <- as.numeric(parts[(n+2):length(parts)])
    
    comp_matrix <- matrix(values, nrow = n, ncol = n, byrow = TRUE)
    rownames(comp_matrix) <- channels
    colnames(comp_matrix) <- channels
    
    return(flowCore::compensation(comp_matrix))
  }
  
  if (is.list(comp_data)) {
    # FlowJo v11 JSON format
    if (!is.null(comp_data$matrix) && !is.null(comp_data$parameters)) {
      matrix_data <- comp_data$matrix
      params <- comp_data$parameters
      
      if (is.list(matrix_data)) {
        matrix_data <- do.call(rbind, matrix_data)
      }
      
      comp_matrix <- as.matrix(matrix_data)
      rownames(comp_matrix) <- params
      colnames(comp_matrix) <- params
      
      return(flowCore::compensation(comp_matrix))
    }
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
      stop("Compensation matrix must have row and column names (channel names)")
    }
    return(flowCore::compensation(comp))
  }
  
  if (is.data.frame(comp)) {
    return(flowCore::compensation(as.matrix(comp)))
  }
  
  stop("Invalid compensation object. Must be compensation, matrix, or data.frame")
}


#' Map Compensation Channel Names to Cytoframe Parameter Names
#'
#' flowCore's compensation() function sanitizes channel names (e.g., "/" -> "_", spaces -> ".").
#' This function creates a mapping between sanitized compensation names and original cytoframe
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
  # The compensation matrix now has "Comp-" prefixed names (e.g., "Comp-FITC-A")
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
  for (i in seq_along(comp_names_base)) {
    if (is.null(name_mapping[[comp_names_base[i]]])) {
      warning("Could not map compensation channel '", comp_names_base[i], "' to any parameter")
    }
  }

  # Apply mapping to create new compensation matrix with cytoframe parameter names
  # (without "Comp-" prefix, so flowCore::compensate can match them)
  mapped_names <- apply_param_mapping(comp_names_base, name_mapping, on_no_match = "keep")

  # Create new compensation matrix with mapped names
  mapped_comp <- comp_matrix
  colnames(mapped_comp) <- mapped_names
  rownames(mapped_comp) <- mapped_names

  flowCore::compensation(mapped_comp)
}











