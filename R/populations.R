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

#' @title Population Functions for FlowJo v11
#' @name populations
#' @keywords internal
#' @importFrom flowWorkspace gs_pop_add gs_get_pop_paths gs_pop_get_gate recompute markernames
NULL

#' Sanitize Population Name for flowWorkspace
#'
#' flowWorkspace uses '/' as the path separator in population paths.  Population
#' names imported from FlowJo may contain '/' (e.g. "CD45, l/d subset"), which
#' causes flowWorkspace to rewrite them inconsistently.  This helper normalizes
#' such names consistently before nodes are added to a GatingSet.
#'
#' @param name Character vector or list of population names.
#' @return Sanitized name of the same structure as the input.
#' @keywords internal
sanitize_population_name <- function(name) {
  sanitize_one <- function(x) {
    if (is.character(x)) {
      gsub("/", ":", x, fixed = TRUE)
    } else {
      x
    }
  }

  if (is.list(name)) {
    lapply(name, sanitize_one)
  } else {
    sanitize_one(name)
  }
}

#' Sanitize Path While Preserving flowWorkspace Separator
#'
#' Population paths use '/' as the separator, but individual population names may
#' contain '/' which must be replaced by ':'.  This helper sanitizes only the
#' name components, leaving the path separator '/' intact.
#'
#' @param path Character vector of path strings (e.g. "Ungated/Lymphocytes/CD45, l/d subset").
#' @return Path strings with name-level '/' replaced by ':'.
#' @keywords internal
sanitize_path_with_separator <- function(path) {
  if (is.null(path)) return(NULL)
  # Split on the path separator, sanitize each component, then rejoin.
  vapply(path, function(p) {
    parts <- strsplit(p, "/", fixed = TRUE)[[1]]
    parts <- sanitize_population_name(parts)
    paste(parts, collapse = "/")
  }, character(1))
}

#' Add Populations to GatingSet
#'
#' Adds population nodes and gates to GatingSet based on gating trees
#'
#' @param gs GatingSet object
#' @param gating_trees List of gating trees (one per sample)
#' @param gates List of gate objects
#' @param sample_uuids Vector of sample UUIDs
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter names
#'   when matching to GatingSet parameters? Default TRUE.
#' @param verbose Logical. Print verbose messages? Default FALSE.
#' @keywords internal
#' @importFrom flowWorkspace gs_pop_add
add_populations_to_gatingset <- function(gs, gating_trees, gates, sample_uuids,
                                         strip_comp_prefix = TRUE, verbose = FALSE) {

  for (i in seq_along(sample_uuids)) {
    sample_uuid <- sample_uuids[i]
    tree        <- gating_trees[[i]]
    gh          <- gs[[i]]

    deferred <- new.env(parent = emptyenv())
    deferred$gates <- list()

    add_population_node(gh, tree, gates, sample_uuid,
                        strip_comp_prefix = strip_comp_prefix,
                        verbose = verbose, deferred = deferred)

    # Retry loop: keep trying until no more progress (handles chained logic gates)
    repeat {
      if (length(deferred$gates) == 0) break
      n_before     <- length(deferred$gates)
      new_deferred <- new.env(parent = emptyenv())
      new_deferred$gates <- list()

      for (lg in deferred$gates) {
        add_population_node(gh, lg$node, gates, sample_uuid,
                            parent = lg$parent,
                            strip_comp_prefix = strip_comp_prefix,
                            verbose = verbose, deferred = new_deferred)
      }

      deferred <- new_deferred
      if (length(deferred$gates) == n_before) break   # no progress — give up
    }

    # Warn only about permanently unresolvable logical gates
    for (lg in deferred$gates) {
      node_name <- if (is.list(lg$node$name)) unlist(lg$node$name) else lg$node$name
      component_names <- lg$node$logical_gate_info$combined_populations
      all_paths <- tryCatch(flowWorkspace::gs_get_pop_paths(gh), error = function(e) character(0))
      missing <- setdiff(component_names, sub(".*/", "", all_paths))
      warning("Could not resolve all component paths for logical gate: ", node_name, "\n",
              "  Missing components: ", paste(missing, collapse = ", "))
    }
  }
}

#' Compare and Adjust Gate Parameter Names to Match a GatingHierarchy
#'
#' Ensures that the parameter names stored inside a \code{flowCore} gate object
#' are consistent with the parameter names present in a \code{GatingHierarchy}
#' after compensation and transformations have been applied.
#'
#' The function performs two distinct tasks:
#' \enumerate{
#'   \item \strong{Parameter-name remapping.}  Gate parameter names extracted
#'     from the FlowJo workspace may differ from the names actually present in
#'     the \code{GatingHierarchy} (e.g., a gate may reference
#'     \code{"Comp-APC-Ax700-A"} while the compensated \code{cytoframe} column
#'     is named \code{"APC-Ax700-A"}, or \code{"/"} characters may have been
#'     sanitized to \code{"_"} by \code{flowCore}).  The function resolves these
#'     discrepancies via \code{\link{map_gate_params_to_gh}} and then rewrites
#'     the gate's internal parameter slots with
#'     \code{\link{update_gate_param_names}}.
#'   \item \strong{Transformation-space adjustment (currently disabled).}  A
#'     code path exists—guarded by \code{if (FALSE)}—that would convert gate
#'     coordinates from raw data space to transformed (display) space using the
#'     transformers registered in the \code{GatingHierarchy}.  This path is
#'     disabled because the current pipeline calls
#'     \code{\link[=extract_all_gates]{extract_all_gates}} with
#'     \code{use_transformed_coords = FALSE}, so gate coordinates are already
#'     expressed in the same space as the (untransformed) \code{cytoframe} data.
#' }
#'
#' The set of parameter names considered for matching is the union of:
#' \itemize{
#'   \item names returned by \code{\link[flowWorkspace]{gh_get_transformations}}
#'     (non-linear channels only), and
#'   \item all column names of the root \code{cytoframe} (so that linear /
#'     scatter parameters such as \code{FSC-A} and \code{SSC-A} are also
#'     available for matching).
#' }
#'
#' If no gate parameter can be mapped to any \code{GatingHierarchy} parameter
#' the original gate object is returned unchanged.
#'
#' @param gh A \code{GatingHierarchy} object (single sample).  Used to retrieve
#'   both the registered transformers (\code{gh_get_transformations}) and the
#'   full list of column names from the root population data frame.
#' @param gate_obj A \code{flowCore} gate object (e.g., \code{rectangleGate},
#'   \code{polygonGate}, \code{ellipsoidGate}, \code{quadGate}).  Its internal
#'   parameter names will be remapped to match \code{gh}.
#' @param strip_comp_prefix Logical. When \code{TRUE} (default), the
#'   \code{"Comp-"} prefix is stripped from gate parameter names before
#'   attempting to match them against the \code{GatingHierarchy} parameter
#'   names.  Set to \code{FALSE} when gate names are expected to already include
#'   the \code{"Comp-"} prefix (i.e., compensation was applied externally with
#'   matching names).
#'
#' @return The input \code{gate_obj} with its internal parameter names updated
#'   to match those in \code{gh}.  All gate coordinate slots (\code{@min},
#'   \code{@max}, \code{@boundaries}, \code{@boundary}, \code{@mean},
#'   \code{@cov}) are updated consistently via
#'   \code{\link{update_gate_param_names}}.  If no mapping is found, the
#'   original object is returned unmodified.
#'
#' @seealso
#'   \code{\link{map_gate_params_to_gh}} for the underlying name-mapping logic,
#'   \code{\link{update_gate_param_names}} for the slot-level renaming,
#'   \code{\link{apply_transforms_to_gate}} for the (currently disabled)
#'   coordinate-transformation step,
#'   \code{\link{add_population_node}} which calls this function before
#'   passing a gate to \code{\link[flowWorkspace]{gs_pop_add}}.
#'
#' @keywords internal
adjust_gate_transformations <- function(gh, gate_obj, strip_comp_prefix = TRUE) {

  # Get transformations from gating hierarchy
  gh_trans <- flowWorkspace::gh_get_transformations(gh)
  gh_param_names <- names(gh_trans)

  # Get parameters from the gate
  gate_params <- flowCore::parameters(gate_obj)

  # Linear transforms are filtered out before applying permanent transforms,
  # so gh_get_transformations() only returns channels with non-linear transforms.
  # Always include all flowFrame column names so that linear/scatter parameters
  # (e.g. FSC-A) are available for matching gate parameters.
  gh_param_names <- unique(c(
    gh_param_names,
    colnames(flowWorkspace::gh_pop_get_data(gh, "root"))
  ))

  # Map gate parameter names to GatingSet parameter names
  # Gate params may have "Comp-" prefix or different sanitization
  mapped_params <- map_gate_params_to_gh(gate_params, gh_param_names, strip_comp_prefix = strip_comp_prefix)

  # Check if we need to adjust transformations
  # With the current pipeline the cytoframe is kept in raw space and gate
  # coordinates are already converted to raw space by display_to_raw().  There
  # is therefore no need to forward-transform gate coordinates before gating.
  needs_adjustment <- FALSE
  trans_to_apply <- list()

  if (FALSE) { # nocov start
    # Disabled: gate coordinates are converted to raw data space by display_to_raw().
    for (i in seq_along(gate_params)) {
      gate_param <- gate_params[i]
      gh_param <- mapped_params[[gate_param]]

      # Skip if no mapping found
      if (is.null(gh_param)) {
        warning("Could not map gate parameter '", gate_param, "' to GatingSet parameters")
        next
      }

      # Get transformation from hierarchy using mapped name
      gh_trans_func <- gh_trans[[gh_param]]

      # Skip if no transformation in hierarchy
      if (is.null(gh_trans_func)) {
        next
      }

      # Get transformation type from hierarchy
      gh_trans_type <- attr(gh_trans_func, "type")

      # If hierarchy has a transformation (not "none"), we need to apply it to gate coords
      if (!is.null(gh_trans_type)) {
        # Gate coordinates are in raw space, need to convert to transformed space
        needs_adjustment <- TRUE
        trans_to_apply[[gh_param]] <- gh_trans_func

        if (.pkgenv$verbose) {
          message("Parameter '", gate_param, "' -> '", gh_param, "' needs adjustment:")
          message("  Gate coords in: raw data space")
          message("  Hierarchy expects: ", gh_trans_type, " transformed space")
        }
      }
    }
  } # nocov end

  # If no mapping found at all, return original gate
  if (length(mapped_params) == 0 || all(vapply(mapped_params, is.null, logical(1)))) {
    return(gate_obj)
  }

  # Always update parameter names to match flowFrame (critical for compensated data)
  gate_obj <- update_gate_param_names(gate_obj, mapped_params)

  # If no transformation adjustment needed, return gate with updated names
  if (!needs_adjustment) {
    return(gate_obj)
  }

  # Apply transformations to convert from raw space to transformed space
  gate_obj_adjusted <- apply_transforms_to_gate(gate_obj, trans_to_apply)

  return(gate_obj_adjusted)
}


#' Map gate parameter names to GatingSet parameter names
#'
#' Gate parameter names may have "Comp-" prefix or different sanitization.
#' This function maps them to the actual GatingSet parameter names.
#'
#' @param gate_params Character vector of gate parameter names
#' @param gh_param_names Character vector of GatingSet parameter names
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter names?
#'   Default TRUE. Set to FALSE if compensation has already been applied with
#'   matching parameter names.
#' @param sanitize_slashes Logical. Replace "/" with "_" in parameter names?
#'   Default TRUE (matches flowCore behavior). Set to FALSE if you want to preserve
#'   "/" in marker names (e.g., "CD3/CD4" stays as-is).
#' @return Named list mapping gate params to GatingSet params
#' @keywords internal
map_gate_params_to_gh <- function(gate_params, gh_param_names, strip_comp_prefix = TRUE, sanitize_slashes = TRUE) {
  # Use unified parameter name mapping
  # Gate params may have "Comp-" prefix from flowCore compensation
  # and "/" may be sanitized to "_"
  map_param_names(
    source_names = gate_params,
    target_names = gh_param_names,
    strip_comp_prefix = strip_comp_prefix,
    case_insensitive = FALSE,
    sanitize_slashes = sanitize_slashes
  )
}

#' Update gate parameter names to match flowFrame
#'
#' After compensation, flowFrame parameter names may have "Comp-" prefix.
#' This function updates gate parameter names to match.
#'
#' @param gate_obj flowCore gate object
#' @param mapped_params Named list mapping old param names to new param names
#' @return Gate object with updated parameter names
#' @keywords internal
update_gate_param_names <- function(gate_obj, mapped_params) {
  
  # Filter out NULL mappings
  valid_mapping <- mapped_params[!vapply(mapped_params, is.null, logical(1))]
  
  if (length(valid_mapping) == 0) {
    return(gate_obj)
  }
  
  # Handle different gate types
  if (inherits(gate_obj, "rectangleGate")) {
    # Update min/max names
    old_names <- names(gate_obj@min)
    new_names <- vapply(old_names, function(n) {
      if (n %in% names(valid_mapping)) valid_mapping[[n]] else n
    }, character(1))
    names(gate_obj@min) <- new_names
    names(gate_obj@max) <- new_names
    
  } else if (inherits(gate_obj, "quadGate")) {
    # Update boundary names - critical for quadGate!
    old_names <- names(gate_obj@boundary)
    new_names <- vapply(old_names, function(n) {
      if (n %in% names(valid_mapping)) valid_mapping[[n]] else n
    }, character(1))
    names(gate_obj@boundary) <- new_names
    
    if (.pkgenv$verbose) {
      message("  Updated quadGate boundary names: ", paste(old_names, collapse = ", "),
              " -> ", paste(new_names, collapse = ", "))
    }
    
  } else if (inherits(gate_obj, "polygonGate")) {
    # Update boundary column names
    old_names <- colnames(gate_obj@boundaries)
    new_names <- vapply(old_names, function(n) {
      if (n %in% names(valid_mapping)) valid_mapping[[n]] else n
    }, character(1))
    colnames(gate_obj@boundaries) <- new_names
    
  } else if (inherits(gate_obj, "ellipsoidGate")) {
    # Update mean and covariance names
    old_names <- names(gate_obj@mean)
    new_names <- vapply(old_names, function(n) {
      if (n %in% names(valid_mapping)) valid_mapping[[n]] else n
    }, character(1))
    names(gate_obj@mean) <- new_names
    colnames(gate_obj@cov) <- new_names
    rownames(gate_obj@cov) <- new_names
  }
  
  # Update the parameters slot
  flowCore::parameters(gate_obj) <- unlist(valid_mapping)
  
  return(gate_obj)
}

#' Apply Transformation Functions to Gate Coordinates
#'
#' Converts the coordinate values stored inside a \code{flowCore} gate object
#' from one data space to another by evaluating a set of named transformation
#' functions.  This is typically used to move gate coordinates from raw
#' (instrument) space into a transformed (display) space so that they align with
#' the transformed \code{cytoframe} data held in a \code{GatingHierarchy}.
#'
#' Each gate type stores its defining coordinates in different slots; the
#' function handles all four common types:
#' \describe{
#'   \item{\code{rectangleGate}}{The \code{@min} and \code{@max} named numeric
#'     vectors are transformed channel-by-channel for every parameter that
#'     appears in \code{trans_list}.}
#'   \item{\code{quadGate}}{The \code{@boundary} named numeric vector (the two
#'     divider positions, one per axis) is transformed for every matching
#'     parameter.}
#'   \item{\code{polygonGate}}{Each column of the \code{@boundaries} matrix that
#'     matches a parameter in \code{trans_list} is transformed element-wise,
#'     preserving the polygon vertex order.}
#'   \item{\code{ellipsoidGate}}{The \code{@mean} vector is transformed
#'     parameter-by-parameter.  A \code{warning} is emitted because applying a
#'     non-linear transform to the mean without also adjusting the covariance
#'     matrix \code{@cov} does not correctly preserve the ellipsoid shape; the
#'     covariance matrix is left unchanged.}
#' }
#'
#' Parameters present in the gate but absent from \code{trans_list} (e.g.,
#' scatter channels with a linear or identity transform) are left unchanged.
#' If \code{trans_list} is empty the original gate object is returned
#' immediately without modification.
#'
#' @note
#' This function is currently \strong{not called} by the active pipeline.  Gate
#' coordinates are extracted in transformed (display) space by
#' \code{\link{extract_all_gates}} when \code{use_transformed_coords = TRUE}
#' (the default), making an additional in-place transformation of gate
#' coordinates unnecessary.  The function is retained for use cases where gate
#' coordinates must be converted programmatically after extraction.
#'
#' @param gate_obj A \code{flowCore} gate object whose coordinate slots will be
#'   transformed.  Supported classes: \code{rectangleGate},
#'   \code{quadGate}, \code{polygonGate}, \code{ellipsoidGate}.  Objects of
#'   other classes are returned unmodified.
#' @param trans_list A named list of transformation functions, where each name
#'   is a \code{GatingHierarchy} parameter name (e.g., \code{"APC-Ax700-A"})
#'   and each value is a single-argument function that maps numeric raw values
#'   to transformed values.  Typically a subset of the list returned by
#'   \code{\link[flowWorkspace]{gh_get_transformations}}.  An empty list causes
#'   the function to return \code{gate_obj} unchanged.
#'
#' @return The input \code{gate_obj} with coordinate slots updated in place for
#'   all parameters found in \code{trans_list}.  The object class and all other
#'   slots are preserved.
#'
#' @seealso
#'   \code{\link{adjust_gate_transformations}} which determines which
#'   transformations need to be applied and calls this function,
#'   \code{\link{extract_all_gates}} for gate extraction with optional
#'   coordinate-space control via \code{use_transformed_coords},
#'   \code{\link[flowWorkspace]{gh_get_transformations}} for retrieving
#'   the transformation functions registered in a \code{GatingHierarchy}.
#'
#' @keywords internal
apply_transforms_to_gate <- function(gate_obj, trans_list) {
  
  if (length(trans_list) == 0) {
    return(gate_obj)
  }
  
  # Handle rectangleGate
  if (inherits(gate_obj, "rectangleGate")) {
    
    for (param in names(trans_list)) {
      trans_func <- trans_list[[param]]
      
      # Apply transformation to min and max values
      if (param %in% names(gate_obj@min)) {
        old_val <- gate_obj@min[param]
        gate_obj@min[param] <- trans_func(old_val)
        if (.pkgenv$verbose) {
          message("  Transformed ", param, " min: ", old_val, " -> ", gate_obj@min[param])
        }
      }
      
      if (param %in% names(gate_obj@max)) {
        old_val <- gate_obj@max[param]
        gate_obj@max[param] <- trans_func(old_val)
        if (.pkgenv$verbose) {
          message("  Transformed ", param, " max: ", old_val, " -> ", gate_obj@max[param])
        }
      }
    }
    
  } else if (inherits(gate_obj, "quadGate")) {
    # browser() # nocov
    # Quadrant gates have boundary (divider positions) that need transformation
    # Note: slot is "boundary" (singular), not "boundaries"
    boundary <- gate_obj@boundary
    
    for (param in names(trans_list)) {
      if (param %in% names(boundary)) {
        trans_func <- trans_list[[param]]
        old_val <- boundary[param]
        boundary[param] <- trans_func(old_val)
        
        if (.pkgenv$verbose) {
          message("  Transformed ", param, " quad divider: ", old_val, " -> ", boundary[param])
        }
      }
    }
    
    gate_obj@boundary <- boundary
    
  } else if (inherits(gate_obj, "polygonGate")) {
    
    boundaries <- gate_obj@boundaries
    
    for (param in names(trans_list)) {
      if (param %in% colnames(boundaries)) {
        trans_func <- trans_list[[param]]
        old_vals <- boundaries[, param]
        boundaries[, param] <- trans_func(old_vals)
        
        if (.pkgenv$verbose) {
          message("  Transformed ", param, " polygon boundaries")
          message("    Range: ", min(old_vals), "-", max(old_vals), 
                  " -> ", min(boundaries[, param]), "-", max(boundaries[, param]))
        }
      }
    }
    
    gate_obj@boundaries <- boundaries
    
  } else if (inherits(gate_obj, "ellipsoidGate")) {
    
    # For ellipsoid gates, need to transform mean and possibly adjust covariance
    for (param in names(trans_list)) {
      param_idx <- which(names(gate_obj@mean) == param)
      
      if (length(param_idx) > 0) {
        trans_func <- trans_list[[param]]
        old_val <- gate_obj@mean[param_idx]
        gate_obj@mean[param_idx] <- trans_func(old_val)
        
        if (.pkgenv$verbose) {
          message("  Transformed ", param, " ellipse mean: ", old_val, " -> ", gate_obj@mean[param_idx])
        }
        
        # Note: transforming an ellipsoid properly requires transforming the covariance matrix
        # This is complex and may not preserve the ellipsoid shape
        warning("Ellipsoid gate transformation may not preserve exact shape")
      }
    }
    
  }
  
  return(gate_obj)
}

#' Add Population Node Recursively to a GatingHierarchy
#'
#' Recursively walks a gating tree node and adds each population—together with
#' its associated gate—to a single-sample \code{GatingHierarchy}.  Both regular
#' flow-cytometry gates (rectangle, polygon, ellipsoid, quadrant) and logical
#' (boolean) gates (\code{and}/\code{or}/\code{not}) are supported.
#'
#' For logical gates, all component populations must already exist in the
#' \code{GatingHierarchy} before the boolean expression can be evaluated.  When
#' a component is missing the gate is pushed onto \code{deferred} rather than
#' silently dropped, allowing \code{\link{add_populations_to_gatingset}} to retry
#' in a subsequent pass once earlier populations have been resolved.
#'
#' The function performs the following steps for each node:
#' \enumerate{
#'   \item Sanitizes the population name (replaces \code{/} with \code{:}) to
#'     avoid conflicts with the flowWorkspace path separator.
#'   \item Skips the virtual root/ungated node (already present in every new
#'     \code{GatingSet}).
#'   \item Looks up the gate object in \code{gates} using the key
#'     \code{paste0(definition_uuid, "_", sample_uuid)}.
#'   \item For \strong{regular gates}: verifies marker-name consistency, adjusts
#'     parameter names via \code{\link{adjust_gate_transformations}}, and adds
#'     the gate with \code{\link[flowWorkspace]{gs_pop_add}}.
#'   \item For \strong{logical gates}: resolves each component population to its
#'     path in the hierarchy, builds a \code{booleanFilter} expression, and adds
#'     it with \code{\link[flowWorkspace]{gs_pop_add}}.  If any component path
#'     cannot be resolved the node is pushed to \code{deferred}.
#'   \item Recurses into \code{node$children}.
#' }
#'
#' @param gh A \code{GatingHierarchy} object (single sample) to which populations
#'   are added.
#' @param node A named list representing one node of the gating tree produced by
#'   \code{\link{build_gating_tree}}.  Expected fields:
#'   \describe{
#'     \item{\code{name}}{Character or list. Population name.}
#'     \item{\code{type}}{Character. Gate type, e.g. \code{"rectangle"},
#'       \code{"polygon"}, \code{"and"}, \code{"or"}, \code{"not"},
#'       \code{"root"}.}
#'     \item{\code{definition_uuid}}{Character or list. UUID of the corresponding
#'       population definition—used to look up the gate object.}
#'     \item{\code{logical_gate_info}}{Optional list describing a boolean gate;
#'       present only when \code{type} is \code{"and"}, \code{"or"}, or
#'       \code{"not"}.  Contains \code{operator} and
#'       \code{combined_populations}.}
#'     \item{\code{children}}{Optional list of child nodes.}
#'     \item{\code{parent}}{Character. Full path of the parent population, used
#'       when recursing into children.}
#'   }
#' @param gates Named list of \code{flowCore} gate objects keyed by
#'   \code{paste0(definition_uuid, "_", sample_uuid)}.  Typically the output of
#'   \code{\link{extract_all_gates}}.
#' @param sample_uuid Character scalar. UUID of the sample being processed.
#'   Used to look up the correct gate from \code{gates} and to provide
#'   informative warning messages.
#' @param parent Character scalar. Path of the parent population in \code{gh}
#'   under which the current node will be added.  Use \code{"root"} (the
#'   default) for top-level populations.
#' @param strip_comp_prefix Logical. When \code{TRUE} (default), the
#'   \code{"Comp-"} prefix is stripped from gate parameter names before
#'   matching against the \code{GatingHierarchy} parameters.  Set to
#'   \code{FALSE} when compensation has already been applied with names that
#'   include the prefix.
#' @param verbose Logical. When \code{TRUE}, diagnostic messages are emitted
#'   via \code{\link{message}} at each step (gate lookup, boolean expression
#'   construction, name mapping, etc.).  Default \code{FALSE}.
#' @param deferred An \code{environment} (or \code{NULL}) used to collect
#'   logical gates that could not be added because one or more component
#'   populations are not yet present in \code{gh}.  Callers should pass the
#'   same environment object on all recursive calls so that
#'   \code{\link{add_populations_to_gatingset}} can retry the deferred gates in
#'   a subsequent pass.  When \code{NULL}, unresolvable logical gates raise a
#'   \code{warning} immediately instead of being deferred.
#'
#' @return \code{NULL} invisibly.  The function modifies \code{gh} and
#'   \code{deferred} in place as side effects.
#'
#' @seealso
#'   \code{\link{add_populations_to_gatingset}} for the top-level caller that
#'   manages the deferred-retry loop,
#'   \code{\link{build_gating_tree}} for tree construction,
#'   \code{\link{adjust_gate_transformations}} for parameter-name remapping,
#'   \code{\link{extract_all_gates}} for gate extraction.
#'
#' @keywords internal
#' @importFrom flowWorkspace gs_pop_add gs_get_pop_paths
#' @importFrom flowCore parameters
#' @importFrom magrittr %>%
add_population_node <- function(gh, node, gates, sample_uuid, parent = "root",
                                strip_comp_prefix = TRUE, verbose = FALSE,
                                deferred = NULL) {
  if (.pkgenv$verbose) message(node$name, "\n")
  # if(stringr::str_starts(node$name, "TNF")) {
  # browser() # nocov
  # }
  # Extract node name correctly (it's a list in FlowJo v11)
  node_name <- if (is.list(node$name) && length(node$name) > 0) {
    node$name %>% unlist()
  } else if (is.character(node$name)) {
    node$name
  } else {
    "Unnamed"
  }
  # flowWorkspace uses '/' as the path separator, so population names containing
  # '/' get rewritten inconsistently.  Sanitize once here so parent/child paths
  # built later match the node names actually created in the GatingSet.
  node_name <- sanitize_population_name(node_name)
  if (.pkgenv$verbose) message(node$type)
  # browser() # nocov
  # Skip root node (already exists).  The incoming parent may be a sanitized
  # path (with '/' as path separator) or the literal string 'root'.
  parent_sanitized <- if (identical(parent, "root")) parent else sanitize_path_with_separator(parent)
  if (.pkgenv$verbose) message("parent: ", paste(parent_sanitized, collapse = " : "),"\n")
  if (parent_sanitized[1] == "root" && (node_name[1] == "root" || node_name[1] == "Ungated")) {
    if (!is.null(node$children)) {
      for (child in node$children) {
        add_population_node(
          gh = gh,
          node = child,
          gates = gates,
          sample_uuid = sample_uuid,
          parent = "root",
          strip_comp_prefix = strip_comp_prefix,
          verbose = verbose,
          deferred = deferred
        )
      }
    }
    return()
  }
  
  # Get gate for this population-sample combination
  # Extract definition_uuid correctly (it's also a list)
  definition_uuid <- if (is.list(node$definition_uuid) && length(node$definition_uuid) > 0) {
    node$definition_uuid[[1]]
  } else if (is.character(node$definition_uuid)) {
    node$definition_uuid
  } else {
    NULL
  }
  
  if (is.null(definition_uuid)) {
    warning("No definition UUID found for population: ", node_name)
    return()
  }
  
  gate_key <- paste0(definition_uuid, "_", sample_uuid)
  gate_obj <- gates[[gate_key]]
  is_logica_gate <- !is.null(node$logical_gate_info)
  if (.pkgenv$verbose) message(node_name, " (sample: ", sample_uuid, ")\n")
  # browser() # nocov
  if (is.null(gate_obj) && ! is_logica_gate) {
    # Get available population paths for debugging
    all_pop_paths <- tryCatch(flowWorkspace::gs_get_pop_paths(gh), error = function(e) NULL)
    pop_paths_str <- if (!is.null(all_pop_paths)) {
      paste("\n  Available populations:", paste(all_pop_paths, collapse = "\n  "))
    } else {
      ""
    }
    
    warning("No gate found for population: ", node_name, " (sample: ", sample_uuid, ")\n",
            "  Parent: ", paste(parent, collapse = " : "),
            "  Definition UUID: ", definition_uuid,
            pop_paths_str)
    return()
  }
  if(is_logica_gate){
    # Create boolean filter for logical gate
    # browser() # nocov
    # Create boolean filter for logical gate
    operator <- node$logical_gate_info$operator
    component_names <- node$logical_gate_info$combined_populations
    
    if (.pkgenv$verbose) message("Building logical gate: ", node_name, " (", operator, ")")
    if (.pkgenv$verbose) message("  Components: ", paste(component_names, collapse = ", "))
    
    # Get all existing populations in the gating hierarchy
    all_paths <- flowWorkspace::gs_get_pop_paths(gh)
    if (.pkgenv$verbose) message("  Available paths: ", paste(all_paths, collapse = ", "))
    
    # For each component, find its path in the hierarchy
    component_refs <- vapply(component_names, function(comp_name) {
      # Search for exact match in population names
      matching_paths <- grep(paste0("/", comp_name, "$"), all_paths, value = TRUE)

      if (length(matching_paths) == 0) {
        if (is.null(deferred)) {
          warning("Could not find population '", comp_name, "' for logical gate '", node_name, "'\n",
                  "  Available populations: ", paste(all_paths, collapse = ", "))
        }
        return(NA_character_)
      }

      # Use the first matching path
      path <- matching_paths[1]
      if (.pkgenv$verbose) message("    Found: ", comp_name, " -> ", path)
      return(path)
    }, character(1))

    # Remove unresolved components
    component_refs <- component_refs[!is.na(component_refs)]
    
    if (length(component_refs) < length(component_names)) {
      if (!is.null(deferred)) {
        deferred$gates <- c(deferred$gates,
                            list(list(node = node, parent = parent)))
      } else {
        missing_components <- setdiff(component_names, names(component_refs))
        warning("Could not resolve all component paths for logical gate: ", node_name, "\n",
                "  Missing components: ", paste(missing_components, collapse = ", "), "\n",
                "  Resolved: ", paste(names(component_refs), collapse = ", "))
      }
      return()
    }
    
    # Clean paths - remove leading slash as per booleanFilter examples
    clean_refs <- gsub("^/", "", component_refs)
    
    # Build boolean expression - NO SPACES!
    if (operator == "and") {
      bool_expr <- paste(clean_refs, collapse = "&")  # No spaces
    } else if (operator == "or") {
      bool_expr <- paste(clean_refs, collapse = "|")  # No spaces
    } else if (operator == "not") {
      bool_expr <- paste0("!", clean_refs[1])  # No spaces
    } else {
      warning("Unknown logical operator '", operator, "' for gate: ", node_name,
              ". Expected 'and', 'or', or 'not'")
      return()
    }
    
    if (.pkgenv$verbose) message("  Boolean expression: ", bool_expr)
    tryCatch({
      # Use the programmatic approach from the documentation
      # Create as symbol and substitute into booleanFilter call
      call_expr <- substitute(booleanFilter(v), list(v = as.symbol(bool_expr)))
      bool_filter <- eval(call_expr)
      
      if (.pkgenv$verbose) message("  Created filter: ", class(bool_filter))
      
      flowWorkspace::gs_pop_add(
        gh,
        bool_filter,
        parent = parent,
        name = node_name
      )
      
      if (.pkgenv$verbose) message("  Successfully added logical gate: ", node_name)
      
      # Recompute immediately to verify it works
      # flowWorkspace::recompute(gh)
      # message("  Recomputed successfully")
      
    }, error = function(e) {
      warning("Failed to add logical gate '", node_name, "': ", e$message, "\n",
              "  Parent: ", paste(parent, collapse = " : "), "\n",
              "  Components: ", paste(component_refs, collapse = ", "), "\n",
              "  Expression: ", bool_expr)
    })
    
    # Recursively add children
    if (!is.null(node$children)) {
      for (child in node$children) {
        # Use the path stored in the tree when available; otherwise build it
        # consistently from sanitized names.
        child_parent <- if (!is.null(child$parent) && is.character(child$parent)) {
          sanitize_population_name(child$parent)
        } else {
          paste0(parent, "/", node_name)
        }
        add_population_node(
          gh = gh,
          node = child,
          gates = gates,
          sample_uuid = sample_uuid,
          parent = child_parent,
          strip_comp_prefix = strip_comp_prefix,
          verbose = verbose,
          deferred = deferred
        )
      }
    }
  }  else{
    # And gate:
    # those are the populations where this definition should be applied to 
    node$pop_def$children$populations
    # Add population node
    #
    tryCatch({
      # for now we ignore tsne gates
      if(!startsWith(node_name[1], "tsne")){
        # name for quadrant should be of length 4
        if (.pkgenv$verbose) message("parent: ", parent, " ", node_name[1], "\n")
        # browser() # nocov
        # this seems to be working for the current case but should
        if(inherits(gate_obj, "quadGate")){
          node_name <- node_name[c(3,4,2,1)]
        }
        # quad gate is tried to be added multiple times.
        # Get flowFrame parameter names from the GatingHierarchy
        flowframe_params <- markernames(gh)
        gate_params <- flowCore::parameters(gate_obj)
        
        if (.pkgenv$verbose) {
          message("Verifying marker name consistency for gate: ", node_name[1])
          message("  Gate params: ", paste(gate_params, collapse = ", "))
          message("  FlowFrame params: ", paste(flowframe_params, collapse = ", "))
        }
        
        # Verify marker name consistency
        verification <- verify_gate_marker_names(
          gate_obj = gate_obj,
          flowframe_params = flowframe_params,
          gate_source = paste(node_name, collapse = "/")
        )
        
        if (!verification$valid && .pkgenv$verbose) {
          for (warn in verification$warnings) {
            message("  Note: ", warn)
          }
        }
        
        # Adjust gate transformations if needed
        gate_obj_adjusted <- adjust_gate_transformations(gh, gate_obj, strip_comp_prefix = strip_comp_prefix)
        
        # Verify marker names again after adjustment
        if (.pkgenv$verbose) {
          gate_params_adjusted <- flowCore::parameters(gate_obj_adjusted)
          message("  Gate params after adjustment: ", paste(gate_params_adjusted, collapse = ", "))
        }
        
        tryCatch({
          flowWorkspace::gs_pop_add(
            gh,
            gate = gate_obj_adjusted,
            parent = parent,
            name = node_name
          )
          if (.pkgenv$verbose) {
            message("  Successfully added gate: ", paste(node_name, collapse = "/"))
          }
        }, error = function(e) {
          # Ignore "already exists" errors - this can happen with quadGates
          # when the population was already added in a previous step
          if (grepl("already exists", e$message, ignore.case = TRUE)) {
            if (.pkgenv$verbose) {
              message("  Population already exists, skipping: ", paste(node_name, collapse = "/"))
            }
          } else {
            warning("Failed to add population '", paste(node_name, collapse = "/"), "': ", e$message, "\n",
                    "  Parent: ", paste(parent, collapse = " : "), "\n",
                    "  Gate type: ", class(gate_obj_adjusted)[1])
          }
        })
      }
      
      # Recursively add children
      if (!is.null(node$children)) {
        for (child in node$children) {
          # child$parent is a full path built from sanitized names with '/' as the
          # flowWorkspace separator.  Strip the leading "Ungated/" or "root" prefix
          # so we get the parent path relative to the root of the GatingSet.
          child_parent <- sub("^(root|Ungated)/", "", child$parent)
          # The result is already sanitized; call sanitize_population_name only to
          # handle any stray name-level '/' that might remain.
          child_parent <- sanitize_path_with_separator(child_parent)
          add_population_node(
            gh = gh,
            node = child,
            gates = gates,
            sample_uuid = sample_uuid,
            parent = child_parent,
            strip_comp_prefix = strip_comp_prefix,
            verbose = verbose,
            deferred = deferred
          )
        }
      }
    }, error = function(e) {
      warning("Failed to process gate for population '", paste(node_name, collapse = "/"), "': ", e$message, "\n",
              "  Parent: ", paste(parent, collapse = " : "), "\n",
              "  Sample UUID: ", sample_uuid)
    })
  }
}

