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

#' @title Helper Functions for the Package
#' @name helpers
#' @keywords internal
NULL

#' Safe Numeric/Integer Coercion Without Warnings
#'
#' BiocCheck discourages \code{suppressWarnings(as.numeric(...))}; these helpers
#' validate the string form with a regular expression first so that
#' \code{as.numeric}/\code{as.integer} never produce a coercion warning while
#' producing exactly the same result as the suppressed call.
#'
#' @param x Character vector to coerce
#' @return Numeric/integer vector, NA where x does not denote a number
#' @name safe-coercion
#' @keywords internal
NULL

#' @rdname safe-coercion
as_num_quiet <- function(x) {
  # Regex mirrors R's numeric literal grammar (incl. hex, Inf, NaN); the
  # pre-check keeps as.numeric from warning on non-numeric text.
  rx <- paste0("^[[:space:]]*[+-]?(Inf(inity)?|NaN|(0[xX][0-9a-fA-F]+",
               "|[0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?)[[:space:]]*$")
  ok <- grepl(rx, x)
  out <- rep(NA_real_, length(x))
  if (any(ok)) {
    # Values reaching as.numeric() already matched the numeric-literal regex,
    # so the conversion cannot warn.
    out[ok] <- as.numeric(x[ok])
  }
  out
}

#' @rdname safe-coercion
as_int_quiet <- function(x) {
  n <- as_num_quiet(x)
  ok <- !is.na(n)
  # as.integer() only warns when the value is out of integer range
  out <- rep(NA_integer_, length(x))
  in_range <- ok & n >= -.Machine$integer.max & n <= .Machine$integer.max
  out[in_range] <- as.integer(n[in_range])
  out
}

#' @importFrom utils packageVersion
NULL

#' Custom Null Coalescing Operator
#'
#' Provides a null coalescing operator that returns the second argument if the first is NULL.
#'
#' @param x First value to check
#' @param y Second value to return if first is NULL
#' @return Either x if not NULL, or y
#' @name null-coalesce
#' @rdname null-coalesce
#' @keywords internal
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Load Example FlowJo v11 Workspace
#'
#' Loads the example FlowJo v11 workspace from inst/extdata/min_test.flowjo
#' for use in tests and examples. This workspace contains a simple gating
#' hierarchy (Ungated -> Lymphocytes -> CD8, CD3 subset) and references
#' sample04.fcs.
#'
#' @return A flowjo11_workspace object
#' @export
#' @examples
#' ws <- load_example_workspace()
#' length(ws$samples)  # Number of samples
load_example_workspace <- function() {
  test_file <- system.file("extdata", "min_test.flowjo", package = "CyFj11")
  if (!file.exists(test_file)) {
    stop("Example FlowJo v11 file not found. Please ensure the package is installed with inst/extdata/min_test.flowjo")
  }

  return(read_flowjo11_workspace(test_file))
}

#' Unified Parameter Name Mapping Function
#'
#' Creates a mapping between source parameter names (e.g., from compensation, gates, transformations)
#' and target parameter names (e.g., cytoframe flowFrame column names).
#'
#' This function handles common naming discrepancies:
#' - "Comp-" prefix (e.g., "Comp-APC-Ax700-A" vs "APC-Ax700-A")
#' - "/" vs "_" sanitization (flowCore converts "/" to "_")
#' - Case sensitivity (optional)
#' - Detector-to-marker description matching (optional, via \code{target_descriptions})
#'
#' @param source_names Character vector of source parameter names
#' @param target_names Character vector of target parameter names (e.g. flowFrame colnames)
#' @param target_descriptions Optional named character vector of descriptions for
#'   each target name (e.g. the marker names returned by \code{markernames()}).
#'   When supplied, source names are also matched against these descriptions.
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from source names? Default FALSE.
#' @param case_insensitive Logical. Case-insensitive matching? Default FALSE.
#' @param sanitize_slashes Logical. Replace "/" with "_" in parameter names?
#'   Default TRUE (matches flowCore behavior). Set to FALSE if you want to preserve
#'   "/" in marker names (e.g., "CD3/CD4" stays as-is).
#' @return Named list mapping source names to target names. Unmapped names have NULL values.
#' @keywords internal
map_param_names <- function(source_names,
                            target_names,
                            target_descriptions = NULL,
                            strip_comp_prefix = FALSE,
                            case_insensitive = FALSE,
                            sanitize_slashes = TRUE) {

  if (is.null(source_names) || length(source_names) == 0) {
    return(list())
  }

  if (is.null(target_names) || length(target_names) == 0) {
    return(setNames(lapply(source_names, function(x) NULL), source_names))
  }

  # Sanitization function
  sanitize_name <- function(x) {
    # Strip "Comp-" prefix if requested
    if (strip_comp_prefix) {
      x <- sub("^Comp-", "", x)
    }
    # Replace "/" with "_" (flowCore compensation sanitization)
    if (sanitize_slashes) {
      x <- gsub("/", "_", x)
    }
    # Apply case normalization if requested
    if (case_insensitive) {
      x <- tolower(x)
    }
    x
  }

  # Create sanitized versions for matching
  sanitized_source <- vapply(source_names, sanitize_name, character(1))
  sanitized_target <- vapply(target_names, sanitize_name, character(1))

  # Prepare optional description matching.
  # target_descriptions should be a named vector: target_name -> description/marker
  has_descriptions <- !is.null(target_descriptions) && length(target_descriptions) > 0
  if (has_descriptions) {
    # Ensure descriptions are named by target_names; if names are missing,
    # assume target_descriptions is in the same order as target_names.
    desc_names <- names(target_descriptions)
    if (is.null(desc_names)) {
      if (length(target_descriptions) == length(target_names)) {
        names(target_descriptions) <- target_names
      } else {
        has_descriptions <- FALSE
      }
    }
    if (has_descriptions) {
      sanitized_desc <- vapply(target_descriptions, sanitize_name, character(1))
    }
  }

  # Create mapping
  mapping <- list()

  for (i in seq_along(source_names)) {
    source_name <- source_names[i]
    sanitized_s <- sanitized_source[i]

    # Direct match (original names)
    if (source_name %in% target_names) {
      mapping[[source_name]] <- source_name
      next
    }

    # Match via sanitized names
    match_idx <- which(sanitized_target == sanitized_s)

    if (length(match_idx) > 0) {
      # Use the original target name
      mapping[[source_name]] <- target_names[match_idx[1]]
      next
    }

    # Match via descriptions (e.g. detector name against marker name)
    if (has_descriptions) {
      match_idx <- which(sanitized_desc == sanitized_s)
      if (length(match_idx) > 0) {
        mapping[[source_name]] <- names(target_descriptions)[match_idx[1]]
        next
      }
    }

    # No match found
    mapping[[source_name]] <- NULL
  }

  return(mapping)
}

#' Apply Parameter Name Mapping to Names Vector
#'
#' Applies a parameter name mapping to convert source names to target names.
#'
#' @param source_names Character vector of source names
#' @param mapping Named list from \code{map_param_names()}
#' @param on_no_match Character. What to do with unmapped names:
#'   "keep" - keep original name (default)
#'   "drop" - remove from result
#'   "warn" - keep with warning
#' @return Character vector of mapped names
#' @keywords internal
apply_param_mapping <- function(source_names, mapping, on_no_match = "keep") {
  result <- character(length(source_names))
  names(result) <- names(source_names)

  for (i in seq_along(source_names)) {
    source_name <- source_names[i]
    mapped <- mapping[[source_name]]

    if (is.null(mapped)) {
      if (on_no_match == "keep") {
        result[i] <- source_name
      } else if (on_no_match == "drop") {
        result[i] <- NA_character_
      } else if (on_no_match == "warn") {
        warning("Could not map parameter name '", source_name, "' to target names")
        result[i] <- source_name
      }
    } else {
      result[i] <- mapped
    }
  }

  # Remove NA values (dropped names)
  result[!is.na(result)]
}

#' Verify Marker Name Consistency After Gate Creation
#'
#' This function checks that marker names in gates match the parameter names
#' in the GatingSet/flowFrame after compensation and transformations have been applied.
#' It helps detect issues where automated routines may have changed marker names.
#'
#' @param gate_obj A flowCore gate object (rectangleGate, polygonGate, etc.)
#' @param flowframe_params Named character vector of parameter names from the flowFrame.
#'   Typically the output of \code{markernames(gh)}, where names are detector/channel
#'   names and values are descriptive marker names.
#' @param gate_source Character. Description of gate source for error messages
#' @return List with components:
#'   \describe{
#'     \item{valid}{Logical. Are all gate parameters mappable to flowFrame params?}
#'     \item{mapping}{Named list of parameter name mappings}
#'     \item{warnings}{Character vector of warning messages for unmapped params}
#'   }
#' @keywords internal
verify_gate_marker_names <- function(gate_obj, flowframe_params, gate_source = "") {
  # Collect warnings in a list, then combine at the end
  warnings_vec <- character()

  # Get gate parameters
  gate_params <- tryCatch({
    flowCore::parameters(gate_obj)
  }, error = function(e) {
    warnings_vec <- c(warnings_vec, paste("Could not extract parameters from gate:", e$message))
    return(character())
  })

  if (length(gate_params) == 0) {
    return(list(
      valid = FALSE,
      mapping = list(),
      warnings = warnings_vec
    ))
  }

  # Create parameter name mapping.  flowframe_params is a named vector from
  # markernames() with detector/channel names as names and marker names as values.
  # We match gate parameters against both the detector names and the marker names.
  name_mapping <- map_param_names(
    source_names = gate_params,
    target_names = names(flowframe_params),
    target_descriptions = flowframe_params,
    strip_comp_prefix = TRUE,
    case_insensitive = FALSE,
    sanitize_slashes = TRUE
  )

  # Check for unmapped parameters - collect all warnings at once
  unmapped_params <- gate_params[vapply(gate_params, function(param) {
    is.null(name_mapping[[param]])
  }, logical(1))]

  if (length(unmapped_params) > 0) {
    additional_warnings <- vapply(unmapped_params, function(param) {
      msg <- paste0("Gate parameter '", param, "' does not match any flowFrame parameter")
      if (gate_source != "") {
        msg <- paste0(msg, " (gate: ", gate_source, ")")
      }
      msg
    }, character(1))
    warnings_vec <- c(warnings_vec, additional_warnings)
  }

  list(
    valid = length(warnings_vec) == 0,
    mapping = name_mapping,
    warnings = warnings_vec
  )
}

#' Restore Marker Names After Automated Routine Changes
#'
#' When an automated routine (like compensation) changes marker names,
#' this function helps restore the original names for consistency.
#'
#' @param gate_obj A flowCore gate object
#' @param name_mapping Named list mapping current names to original names
#' @return Gate object with restored parameter names
#' @keywords internal
# restore_gate_marker_names <- function(gate_obj, name_mapping) {
#   # Get current gate parameters
#   gate_params <- flowCore::parameters(gate_obj)
# 
#   # Create reverse mapping (target -> source)
#   reverse_mapping <- setNames(names(name_mapping), unlist(name_mapping))
# 
#   # Map gate parameters back to original names
#   original_params <- character(length(gate_params))
#   for (i in seq_along(gate_params)) {
#     original_params[i] <- reverse_mapping[gate_params[i]] %||% gate_params[i]
#   }
# 
#   # Set the original parameter names on the gate
#   flowCore::parameters(gate_obj) <- original_params
# 
#   gate_obj
# }
