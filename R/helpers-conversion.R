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

#' @title Helper Functions for FlowJo v11 Conversion
#' @name helpers-conversion
#' @keywords internal
NULL

#' Get Group Information
#' @param groups Groups list from FlowJo v11 workspace
#' @return Data frame with group info
#' @keywords internal
get_group_info <- function(groups) {
  group_data <- lapply(names(groups), function(uuid) {
    g <- groups[[uuid]]
    data.frame(
      uuid = uuid,
      name = g$definition$name %||% "Unnamed Group",
      n_samples = length(g$results$dataSources %||% list()),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, group_data)
}


#' Filter Samples Based on Subset Argument
#' @keywords internal
filter_samples <- function(sample_uuids, subset, dataSources, keywords) {
  
  # Parse subset argument
  subset_parsed <- try(eval(substitute(subset)), silent = TRUE)
  
  if (inherits(subset_parsed, "try-error")) {
    # Try as expression filter
    keys_df <- extract_keywords_for_samples(sample_uuids, dataSources, keywords)
    subset_parsed <- try({
      filtered <- filter(keys_df, !!enquo(subset))
      filtered$sample_uuid
    }, silent = TRUE)
    
    if (inherits(subset_parsed, "try-error")) {
      stop("Invalid 'subset' argument: ", attr(subset_parsed, "condition")$message)
    }
  }
  
  # Handle different subset types
  if (is.numeric(subset_parsed)) {
    # Numeric indices
    return(sample_uuids[subset_parsed])
  } else if (is.character(subset_parsed)) {
    # Filenames - match to sample UUIDs
    matched <- sapply(subset_parsed, function(fname) {
      idx <- which(sapply(dataSources, function(ds) {
        basename(ds$definition$uri) == fname ||
          ds$definition$customKeywords$`File Name` == fname
      }))
      if (length(idx) > 0) names(dataSources)[idx[1]] else NA
    })
    matched <- matched[!is.na(matched)]
    return(intersect(sample_uuids, matched))
  } else if (is.list(subset_parsed) && "name" %in% names(subset_parsed)) {
    # List with name element
    return(filter_samples(sample_uuids, subset_parsed$name, dataSources, keywords))
  }
  
  # Default: return all
  return(sample_uuids)
}


#' Extract Keywords for Samples
#' @keywords internal
extract_keywords_for_samples <- function(sample_uuids, dataSources, keywords) {
  if (length(keywords) == 0) {
    stop("'keywords' must be specified when using expression-based subset")
  }
  
  keys_list <- lapply(sample_uuids, function(uuid) {
    ds <- dataSources[[uuid]]
    kw <- ds$definition$customKeywords[keywords]
    kw$sample_uuid <- uuid
    kw$filename <- basename(ds$definition$uri)
    as.data.frame(kw, stringsAsFactors = FALSE)
  })
  
  do.call(rbind, keys_list)
}


#' Find Root Population
#' @keywords internal
find_root_population <- function(populations, populationDefinitions, sample_uuid) {
  # Find root population definition (type = "root")
  root_popdefs <- Filter(function(pd) {
    !is.null(pd$definition$type) && pd$definition$type == "root"
  }, populationDefinitions)
  
  if (length(root_popdefs) == 0) {
    stop("Could not find root population definition")
  }
  
  root_popdef_uuid <- names(root_popdefs)[1]
  root_popdef <- root_popdefs[[root_popdef_uuid]]
  
  # In FlowJo v11, we need to find populations that reference this population definition
  # and belong to our sample
  sample_pop_uuid <- NULL
  
  # Iterate through all populations to find one that:
  # 1. References our root population definition
  # 2. Belongs to our sample (through parent relationships)
  # TODO: there is only compoundPopulations but no populationReference, this needs to be verified
  for (pop_uuid in names(populations)) {
    pop <- populations[[pop_uuid]]
    if (!is.null(pop)) {
      # Check if this population references our root population definition
      if (!is.null(pop$populationReference) && pop$populationReference == root_popdef_uuid) {
        # Check if this population belongs to our sample
        parents <- pop$parents
        if (!is.null(parents)) {
          ds_uuid <- parents[['_dataSource']]
          if (is.null(ds_uuid)) {
            ds_uuid <- parents[['dataSource']]
          }
          if (!is.null(ds_uuid) && ds_uuid[[1]] == sample_uuid) {
            sample_pop_uuid <- pop_uuid
            break
          }
        }
      }
    }
  }
  
  if (is.null(sample_pop_uuid)) {
    # If we couldn't find by direct reference, try to find by parentPopulation = NULL
    # which indicates a root population
    for (pop_uuid in names(populations)) {
      pop <- populations[[pop_uuid]]
      if (!is.null(pop) && is.null(pop$parentPopulation)) {
        # Check if this population belongs to our sample
        parents <- pop$parents
        if (!is.null(parents)) {
          ds_uuid <- parents[['_dataSource']]
          if (is.null(ds_uuid)) {
            ds_uuid <- parents[['dataSource']]
          }
          if (!is.null(ds_uuid) && ds_uuid[[1]] == sample_uuid) {
            sample_pop_uuid <- pop_uuid
            break
          }
        }
      }
    }
  }
  
  if (is.null(sample_pop_uuid)) {
    stop("Could not find root population for sample: ", sample_uuid)
  }
  
  sample_pop_uuid
}


#' Build Gating Tree for given Sample
#' @keywords internal
build_gating_tree <- function(sample_uuid, populations, populationDefinitions, root_uuid) {
  
  # First, identify all logical gates before building the tree
  if (.pkgenv$verbose) message("Identifying logical gates...") # nocov
  gate_result <- identify_logical_gates(populations, populationDefinitions)
  logical_gates_info <- gate_result$gates
  populationDefinitions <- gate_result$populationDefinitions  # Use updated version!
  
  if (length(logical_gates_info) > 0) {
    if (.pkgenv$verbose) message(sprintf("Found %d logical gates", length(logical_gates_info))) # nocov
    if (.pkgenv$verbose) { # nocov
      summary_df <- create_logical_gate_summary(logical_gates_info)
      print(summary_df)
    }
  } else {
    if (.pkgenv$verbose) message("No logical gates found") # nocov
  }
  
  # Track visited population UUIDs to prevent infinite loops
  visited <- new.env(parent = emptyenv())
  
  # Recursive function to build tree
  build_node <- function(pop_uuid, parent_path = NULL, parent_pop_uuid = NULL) {
    
    # Prevent infinite recursion
    if (!is.null(visited[[pop_uuid]])) {
      if (.pkgenv$verbose) warning("Circular reference detected for population: ", pop_uuid) # nocov
      return(NULL)
    }
    visited[[pop_uuid]] <- TRUE
    
    pop <- populations[[pop_uuid]]
    if (is.null(pop)) {
      warning("Population not found: ", pop_uuid)
      return(NULL)
    }
    
    # Get population definition
    pop_parents <- pop$parents
    if (is.null(pop_parents)) {
      warning("Population has no parents: ", pop_uuid)
      return(NULL)
    }
    pop_def_uuid <- pop_parents[['populationDefinitions']]
    if (is.null(pop_def_uuid)) {
      warning("Population has no populationDefinitions parent: ", pop_uuid)
      return(NULL)
    }
    
    pop_def <- populationDefinitions[[pop_def_uuid[[1]]]]
    if (is.null(pop_def)) {
      warning("Population definition not found: ", pop_def_uuid)
      return(NULL)
    }
    
    # Extract definition details safely
    node_name <- if (!is.null(pop_def$definition$name)) pop_def$definition$name else "Unnamed"
    # flowWorkspace uses '/' as the path separator, so population names containing
    # '/' get rewritten (e.g. "CD45, l/d subset" -> "CD45, l:d subset").  Sanitize
    # names here so parent/child paths stay consistent.
    node_name <- sanitize_population_name(node_name)
    
    node <- list(
      uuid = pop_uuid,
      name = node_name,
      parent = parent_path,
      data_uuid = pop_parents$dataSources,
      parentDef = pop_parents$populationDefinitions,
      pop_def = pop_def,
      definition_uuid = pop_def_uuid,
      type = pop_def$definition$type,
      kind = pop_def$definition$kind
    )
    
    # Add logical gate information if this is a logical gate
    gate_info <- find_gate_info(pop_uuid, logical_gates_info)
    if (!is.null(gate_info)) {
      if (.pkgenv$verbose) message(sprintf("Adding gate info for: %s", unlist(node_name))) # nocov
      
      node$logical_gate_info <- list(
        operator = gate_info$gate_type,
        combined_populations = gate_info$combined_populations,
        combined_population_uuids = gate_info$combined_population_uuids,
        combined_definition_uuids = gate_info$combined_definition_uuids
      )
      
      # The gateDefinition should already be in pop_def from identify_logical_gates
      # But verify and add if missing
      if (is.null(node$pop_def$definition$gateDefinition)) {
        if (.pkgenv$verbose) message(sprintf("  -> Adding missing gateDefinition for %s", unlist(node_name))) # nocov
        node$pop_def$definition$gateDefinition <- list(
          type = "logical",
          operator = gate_info$gate_type,
          components = gate_info$combined_populations,
          component_uuids = gate_info$combined_population_uuids
        )
      } else {
        if (.pkgenv$verbose) message(sprintf("  -> gateDefinition already exists for %s", unlist(node_name))) # nocov
      }
    }
    
    # Build current path.  Both parts are already sanitized, but keep the path
    # separator as '/' (flowWorkspace convention) while ensuring no stray '/' from
    # names leaks in.
    current_path <- if (is.null(parent_path)) {
      node_name
    } else {
      paste(parent_path, node_name, sep = "/")
    }
    
    # Find children populations
    child_populations <- list()
    
    if (!is.null(pop$children) && !is.null(pop$children$populations)) {
      child_uuids <- pop$children$populations
      
      if (length(child_uuids) > 0) {
        # get order
        # here we need to use populationNumber if available.
        pop_order = c()
        idx=1
        for (child_pop_uuid in child_uuids) {
          child_pop <- populations[[child_pop_uuid]]
          if(is.null(child_pop$definition$populationNumber)){
            pop_order <- c(pop_order, idx)
          } else {
            pop_order <- c(pop_order,child_pop$definition$populationNumber)
          }
          idx = idx + 1
        }
        # browser() # nocov
        for (child_pop_uuid in child_uuids[order(pop_order)]) {
          child_pop <- populations[[child_pop_uuid]]
          pop_num <- child_pop$definition$populationNumber
          if (!is.null(child_pop)) {
            child_sample_uuid <- child_pop$parents[['_dataSource']] %||% child_pop$parents[['dataSource']]
            if (!is.null(child_sample_uuid)) {
              child_populations[[child_pop_uuid]] <- child_pop
            }
          }
        }
      }
    }
    
    # Build children nodes
    if (length(child_populations) > 0) {
      children_list <- list()
      for (child_pop_uuid in names(child_populations)) {
        child_node <- build_node(child_pop_uuid, current_path, pop_uuid)
        if (!is.null(child_node)) {
          children_list[[child_pop_uuid]] <- child_node
        }
      }
      
      children_list <- Filter(Negate(is.null), children_list)
      
      if (length(children_list) > 0) {
        node$children <- unname(children_list)
      }
    }
    
    # Don't remove from visited - this prevents circular references
    # rm(list = pop_uuid, envir = visited)
    
    return(node)
  }
  
  # Build the tree
  if (.pkgenv$verbose) message("Building tree...") # nocov
  tree <- build_node(root_uuid)
  
  # Move logical gates up one level (if any exist)
  if (length(logical_gates_info) > 0) {
    if (.pkgenv$verbose) message("\nMoving logical gates up in hierarchy...") # nocov
    tree <- move_logical_gates_up(tree)
    
    # Remove duplicates at each level (after moving)
    if (.pkgenv$verbose) message("\nRemoving duplicates...") # nocov
    tree <- deduplicate_tree(tree)
    
    # Get updated summary
    summary_after <- summarize_logical_gates(tree)
    if (!is.null(summary_after)) {
      if (.pkgenv$verbose) message("\nLogical gates in final tree:") # nocov
      if (.pkgenv$verbose) print(summary_after) # nocov
    }
  }
  
  # Clean up visited environment
  rm(list = ls(envir = visited), envir = visited)
  
  return(tree)
}
# Main function to identify logical gates from populations and definitions
identify_logical_gates <- function(populations, populationDefinitions) {
  
  # Helper to find population definition by UUID
  find_pop_def_by_uuid <- function(uuid) {
    for (i in seq_along(populationDefinitions)) {
      if (!is.null(populationDefinitions[[i]]$uuid) && 
          populationDefinitions[[i]]$uuid == uuid) {
        return(populationDefinitions[[i]])
      }
    }
    return(NULL)
  }
  
  # Helper to find population by UUID
  find_pop_by_uuid <- function(uuid) {
    for (i in seq_along(populations)) {
      if (!is.null(populations[[i]]$uuid) && 
          populations[[i]]$uuid == uuid) {
        return(populations[[i]])
      }
    }
    return(NULL)
  }
  
  # Find what populations are combined by a logical gate
  find_combined_populations <- function(pop) {
    # The populations combined by the logical gate are in parents$populations
    parent_pop_uuids <- unlist(pop$parents$populations)
    
    if (is.null(parent_pop_uuids) || length(parent_pop_uuids) == 0) {
      return(NULL)
    }
    
    # For each parent population UUID, get its name
    combined_pops <- list()
    
    for (parent_uuid in parent_pop_uuids) {
      parent_pop <- find_pop_by_uuid(parent_uuid)
      
      if (!is.null(parent_pop)) {
        # Get the definition for this population to get the name
        pop_def_uuid <- unlist(parent_pop$parents$populationDefinitions)
        
        if (!is.null(pop_def_uuid) && length(pop_def_uuid) > 0) {
          parent_def <- find_pop_def_by_uuid(pop_def_uuid[1])
          
          if (!is.null(parent_def) && !is.null(parent_def$definition$name)) {
            combined_pops <- c(combined_pops, list(list(
              population_uuid = parent_pop$uuid,
              definition_uuid = parent_def$uuid,
              name = sanitize_population_name(unlist(parent_def$definition$name))
            )))
          }
        }
      }
    }
    
    return(combined_pops)
  }
  
  # Process each population to find logical gates
  # First, identify which populations are logical gates
  logical_gate_indices <- which(vapply(populations, function(pop) {
    pop_def_uuid <- unlist(pop$parents$populationDefinitions)
    if (is.null(pop_def_uuid) || length(pop_def_uuid) == 0) return(FALSE)
    pop_def <- find_pop_def_by_uuid(pop_def_uuid[1])
    if (is.null(pop_def)) return(FALSE)
    !is.null(pop_def$definition$type) && pop_def$definition$type %in% c("and", "or", "not")
  }, logical(1)))

  # Process logical gates and collect results
  results_list <- lapply(logical_gate_indices, function(i) {
    pop <- populations[[i]]
    pop_def_uuid <- unlist(pop$parents$populationDefinitions)[1]
    pop_def <- find_pop_def_by_uuid(pop_def_uuid)
    gate_type <- pop_def$definition$type
    gate_name <- unlist(pop_def$definition$name)

    if (.pkgenv$verbose) message(sprintf("Found logical gate: %s (type: %s, uuid: %s)", # nocov
                                         gate_name, gate_type, pop$uuid))

    # Find combined populations
    combined <- find_combined_populations(pop)

    if (!is.null(combined) && length(combined) > 0) {
      combined_names <- sapply(combined, function(x) x$name)
      if (.pkgenv$verbose) message(sprintf("  - Combines: %s", paste(combined_names, collapse = ", "))) # nocov

      # Add gateDefinition if missing
      if (is.null(pop_def$definition$gateDefinition)) {
        if (.pkgenv$verbose) message(sprintf("  - Adding gateDefinition to populationDefinitions")) # nocov
        populationDefinitions[[pop_def$uuid]]$definition$gateDefinition <<- list(
          type = "logical",
          operator = gate_type,
          components = combined_names,
          component_uuids = sapply(combined, function(x) x$population_uuid)
        )
      }

      list(
        population_uuid = pop$uuid,
        definition_uuid = pop_def$uuid,
        gate_name = gate_name,
        gate_type = gate_type,
        combined_populations = combined_names,
        combined_population_uuids = sapply(combined, function(x) x$population_uuid),
        combined_definition_uuids = sapply(combined, function(x) x$definition_uuid),
        num_components = length(combined)
      )
    } else {
      message("  - Warning: No combined populations found!")
      message(sprintf("  - parents$populations: %s",
                      paste(unlist(pop$parents$populations), collapse = ", ")))
      NULL
    }
  })

  # Remove NULL entries
  results <- results_list[!vapply(results_list, is.null, logical(1))]
  
  # Return both the gate info AND the updated populationDefinitions
  return(list(
    gates = results,
    populationDefinitions = populationDefinitions
  ))
}

# Remove duplicate children at each level based on UUID
deduplicate_tree <- function(tree) {

  deduplicate_node <- function(node) {
    if (!is.list(node)) {
      return(node)
    }

    # Process children if they exist
    if (!is.null(node$children) && length(node$children) > 0) {

      # Use Filter to keep only first occurrence of each UUID
      # For non-list children or children with unique UUIDs, keep them
      seen_uuids <- character()
      keep_indices <- logical(length(node$children))

      for (i in seq_along(node$children)) {
        child <- node$children[[i]]
        if (!is.list(child)) {
          # Keep non-list children
          keep_indices[i] <- TRUE
        } else {
          child_uuid <- child$uuid
          if (is.null(child_uuid) || !(child_uuid %in% seen_uuids)) {
            # First time seeing this UUID - keep it
            if (!is.null(child_uuid)) {
              seen_uuids <- c(seen_uuids, child_uuid)
            }
            keep_indices[i] <- TRUE
          }
        }
      }

      # Filter to unique children, then recursively deduplicate
      unique_children <- node$children[keep_indices]
      unique_children <- lapply(unique_children, function(child) {
        if (is.list(child)) {
          deduplicate_node(child)
        } else {
          child
        }
      })

      node$children <- unique_children
      if (length(node$children) == 0) {
        node$children <- NULL
      }
    }

    return(node)
  }

  deduplicate_node(tree)
}

# Helper function to find gate info for a specific population UUID
find_gate_info <- function(pop_uuid, logical_gates_info) {
  for (gate in logical_gates_info) {
    if (gate$population_uuid == pop_uuid) {
      return(gate)
    }
  }
  return(NULL)
}



# Create summary dataframe
create_logical_gate_summary <- function(logical_gates) {
  if (length(logical_gates) == 0) {
    return(NULL)
  }
  
  df <- data.frame(
    gate_name = sapply(logical_gates, function(x) x$gate_name),
    gate_type = sapply(logical_gates, function(x) x$gate_type),
    population_uuid = sapply(logical_gates, function(x) x$population_uuid),
    num_components = sapply(logical_gates, function(x) x$num_components),
    stringsAsFactors = FALSE
  )
  
  df$combined_populations <- sapply(logical_gates, function(x) {
    paste(x$combined_populations, collapse = " | ")
  })
  
  return(df)
}


# Summarize logical gates in tree
summarize_logical_gates <- function(tree) {

  collect_logical_gates <- function(node, path = "") {
    gates_list <- list()

    if (!is.list(node)) {
      return(gates_list)
    }

    # Build current_path handling vector cases
    current_path <- if (!is.null(node$name)) {
      node_names <- unlist(node$name)
      # Handle when path is empty string(s) or has values
      if (all(path == "")) {
        node_names
      } else {
        # Create all combinations of paths and names
        as.vector(outer(path, node_names, paste, sep = "/"))
      }
    } else {
      path
    }


    if (!is.null(node$type) && node$type %in% c("and", "or", "not")) {
      gate_info <- list(
        path = current_path,
        name = unlist(node$name),
        type = node$type,
        uuid = node$uuid
      )

      if (!is.null(node$logical_gate_info)) {
        gate_info$combined_populations <- node$logical_gate_info$combined_populations
        gate_info$num_components <- length(node$logical_gate_info$combined_populations)
      } else {
        gate_info$num_components <- 0
      }

      gates_list <- list(gate_info)
    }

    if (!is.null(node$children) && is.list(node$children)) {
      child_gates <- lapply(node$children, function(child) {
        collect_logical_gates(child, current_path)
      })
      gates_list <- c(gates_list, unlist(child_gates, recursive = FALSE))
    }

    return(gates_list)
  }
  
  gates <- collect_logical_gates(tree)
  
  if (length(gates) > 0) {
    df <- data.frame(
      path = sapply(gates, function(g) g$path),
      name = sapply(gates, function(g) g$name),
      type = sapply(gates, function(g) g$type),
      num_components = sapply(gates, function(g) g$num_components),
      stringsAsFactors = FALSE
    )
    
    df$combined_populations <- sapply(gates, function(g) {
      if (!is.null(g$combined_populations) && length(g$combined_populations) > 0) {
        paste(g$combined_populations, collapse = " | ")
      } else {
        NA
      }
    })
    
    return(df)
  } else {
    return(NULL)
  }
}

# Move logical gates up to the nearest non-logical ancestor
#
# FlowJo v11 stores boolean (logical) populations as children of the populations
# that appear in their definition.  In flowWorkspace, however, a boolean gate is
# evaluated within its parent population.  This means a gate such as
# "bothNOT = NOT both" must be attached at the same level as "both" (i.e. under
# root) in order to count the complement of "both".  Keeping it as a child of
# "both" would always yield zero events.
#
# This helper walks the raw tree and pulls every logical gate up until its
# parent is a non-logical gate (or root).  Non-logical children are processed
# recursively so that logical gates nested under logical gates are also moved
# up to the correct level.
move_logical_gates_up <- function(tree) {
  
  # Recursive function to process each node.
  # target_ancestor_name is the name of the nearest non-logical ancestor that
  # logical descendants should be attached to.
  process_node <- function(node, target_ancestor_name = NULL) {
    if (!is.list(node) || is.null(node$children) || length(node$children) == 0) {
      return(node)
    }
    
    node_is_logical <- !is.null(node$type) && node$type %in% c("and", "or", "not")
    current_target <- if (!node_is_logical) unlist(node$name) else target_ancestor_name
    
    children_to_keep <- list()
    gates_to_move_up <- list()
    
    for (child in node$children) {
      if (!is.list(child)) {
        children_to_keep <- c(children_to_keep, list(child))
        next
      }
      
      is_logical_gate <- !is.null(child$type) && child$type %in% c("and", "or", "not")
      
      if (is_logical_gate) {
        # Process children first so nested logical gates can bubble up.
        processed_child <- process_node(child, current_target)
        
        # Any logical grandchildren of a logical gate belong at the target
        # ancestor level, not under this gate.
        if (!is.null(processed_child$children)) {
          for (grandchild in processed_child$children) {
            if (is.list(grandchild) &&
                !is.null(grandchild$type) &&
                grandchild$type %in% c("and", "or", "not")) {
              grandchild$moved_from <- unlist(processed_child$name)
              grandchild$parent <- current_target
              gates_to_move_up <- c(gates_to_move_up, list(grandchild))
            }
          }
          
          processed_child$children <- Filter(
            function(g) {
              !(is.list(g) && !is.null(g$type) && g$type %in% c("and", "or", "not"))
            },
            processed_child$children
          )
          if (length(processed_child$children) == 0) {
            processed_child$children <- NULL
          }
        }
        
        processed_child$moved_from <- unlist(processed_child$parent)
        processed_child$parent <- current_target
        
        if (.pkgenv$verbose) { # nocov
          message(sprintf("Moving logical gate '%s' from '%s' to '%s'",
                          unlist(processed_child$name),
                          processed_child$moved_from,
                          processed_child$parent))
        }
        
        if (!node_is_logical) {
          # Keep the logical gate at this non-logical level.
          children_to_keep <- c(children_to_keep, list(processed_child))
        } else {
          # Bubble it up further.
          gates_to_move_up <- c(gates_to_move_up, list(processed_child))
        }
      } else {
        # Non-logical child: recurse, then pull its logical grandchildren up.
        processed_child <- process_node(child, current_target)
        
        child_logical_gates <- list()
        child_regular_children <- list()
        
        if (!is.null(processed_child$children)) {
          for (grandchild in processed_child$children) {
            if (is.list(grandchild) &&
                !is.null(grandchild$type) &&
                grandchild$type %in% c("and", "or", "not")) {
              grandchild$moved_from <- unlist(processed_child$name)
              grandchild$parent <- current_target
              child_logical_gates <- c(child_logical_gates, list(grandchild))
            } else {
              child_regular_children <- c(child_regular_children, list(grandchild))
            }
          }
        }
        
        processed_child$children <- child_regular_children
        if (length(processed_child$children) == 0) {
          processed_child$children <- NULL
        }
        
        children_to_keep <- c(children_to_keep, list(processed_child))
        gates_to_move_up <- c(gates_to_move_up, child_logical_gates)
      }
    }
    
    node$children <- c(children_to_keep, gates_to_move_up)
    if (length(node$children) == 0) {
      node$children <- NULL
    }
    
    return(node)
  }
  
  process_node(tree)
}








get_uuids <- function(tree, uuids=c()){
  uuids = c(uuids, tree$uuid)
  for (idx in seq_along(tree$children)){
    uuids = get_uuids(tree$children[[idx]], uuids)
  }
  return(uuids)
}


#' Create GatingSet from Components - Transform First Approach
#' @keywords internal
#' @importFrom flowWorkspace GatingSet
#' @importFrom magrittr %>%
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter names
#'   when adding populations? Default TRUE. Set to FALSE if compensation has already
#'   been applied and gate names should match the compensated parameter names.
create_gatingset_from_cytoset <- function(cytoset,
                                          gating_trees,
                                          gates,
                                          compensations,
                                          transformations,
                                          sample_uuids,
                                          dataSources,
                                          keywords,
                                          additional.keys,
                                          additional.sampleID,
                                          keyword.ignore.case,
                                          strip_comp_prefix = TRUE) {
  
  actual_sample_count <- length(cytoset)
  # browser() # nocov
  # Create individual GatingHierarchy objects with transformations
  gsList = list()
  
  for (i in seq_len(min(length(sample_uuids), actual_sample_count))) {
    sample_uuid <- sample_uuids[i][[1]]

    if (.pkgenv$verbose) cat("Processing sample", i, "of", actual_sample_count, "\n") # nocov

    # Extract single cytoframe
    cf <- cytoset[[i]]

    # Apply compensation (if available)
    if (!is.null(compensations[[sample_uuid]])) {
      # Store original compensation names (with "Comp-" prefix) before mapping
      comp_orig <- compensations[[sample_uuid]]
      if (methods::is(comp_orig, "compensation")) {
        comp_prefix_names <- colnames(comp_orig@spillover)
      } else {
        comp_prefix_names <- colnames(comp_orig)
      }

      # Map compensation channel names to cytoframe parameter names
      # This handles cases where flowCore sanitizes names (e.g., "/" -> "_")
      comp_mapped <- map_compensation_names(
        compensations[[sample_uuid]],
        colnames(cf)
      )
      cf_comp <- compensate(cf, comp_mapped)

      # After compensation, rename the compensated channels to have "Comp-" prefix.
      # FlowJo 10/11 expects compensated channels to be named "Comp-<channel>"
      # (e.g., "Comp-FITC-A" instead of "FITC-A").
      if (!is.null(comp_prefix_names) && length(comp_prefix_names) > 0) {
        # Get current column names and replace compensated ones with Comp- prefix
        new_names <- as.character(colnames(cf_comp))
        for (comp_name in comp_prefix_names) {
          # Strip "Comp-" prefix to find the base name in the cytoframe
          base_name <- sub("^Comp-", "", comp_name)
          if (base_name %in% new_names) {
            new_names[new_names == base_name] <- comp_name
          }
        }
        # Use flowCore::colnames<- to avoid any S4 method issues
        flowCore::colnames(cf_comp) <- new_names
        cf <- cf_comp
      } else {
        cf <- cf_comp
      }
    }
    
    cs = cytoset()
    cs_add_cytoframe(cs, identifier(cf), cf)
    # Create GatingSet from single sample
    gs_single <- GatingSet(cs)
    
    # Apply transformations (sample-specific)
    # sample_transformations <- transformations[[sample_uuid]] %||% transformations[[1]]
    sample_transformations <- transformations[[sample_uuid]]
    if (is.null(sample_transformations) && length(transformations) > 0) {
      sample_transformations <- transformations[[1]]
    }
    
    # When transform=TRUE, data is transformed to display space and gates are
    # extracted in the same transformed space (via use_transformed_coords).
    # This keeps gating evaluation, visualization, and gate coordinates
    # consistent.
    if (!is.null(sample_transformations) && length(sample_transformations) > 0) {
      # Map transformation channel names to flowFrame parameter names
      # This handles cases where transformation names don't match (e.g., "Comp-APC-Ax700-A" vs "APC-Ax700-A")
      trans_mapped <- map_transformation_names(
        sample_transformations,
        colnames(cf)
      )

      # Keep only non-NULL, non-linear transforms (flowWorkspace-compatible objects)
      trans_apply <- Filter(function(t) {
        !is.null(t) &&
          (is.null(attr(t, "type")) || attr(t, "type") != "linear")
      }, trans_mapped)

      if (length(trans_apply) > 0) {
        tryCatch({
          transList <- flowWorkspace::transformerList(
            from  = names(trans_apply),
            trans = trans_apply
          )
          gs_single <- flowWorkspace::transform(gs_single, transList)
        }, error = function(e) {
          warning(sprintf("Failed to apply transformations for sample %d: %s",
                          i, e$message))
        })
      }
    }
    
    # Extract the GatingHierarchy
    gsList[[i]] <- gs_single
  }
  
  # Combine GatingHierarchy objects into a GatingSet
  if (.pkgenv$verbose) cat("Combining", length(gsList), "GatingHierarchy objects into GatingSet...\n") # nocov
  # this will permanately transformt the data and loose transformation information
  gs = merge_list_to_gs(gsList)
  
  
  
  # Set sample names
  sample_names <- create_sample_names(
    sample_uuids,
    dataSources,
    additional.keys,
    additional.sampleID
  )
  names(gsList) = sample_names
  # tryCatch({
  #   flowWorkspace::sampleNames(gs) <- sample_names
  # }, error = function(e) {
  #   warning("Could not set sample names: ", e$message)
  # })
  
  # Add pData (keywords)
  if (length(keywords) > 0) {
    actual_sample_uuids <- sample_uuids[seq_len(actual_sample_count)]
    pdata <- extract_pdata(actual_sample_uuids, dataSources, keywords, keyword.ignore.case)
    
    if (!is.data.frame(pdata)) {
      pdata <- as.data.frame(pdata, stringsAsFactors = FALSE)
    }
    
    if (nrow(pdata) == actual_sample_count) {
      tryCatch({
        current_rownames <- rownames(flowWorkspace::pData(gs))
        if (length(current_rownames) == nrow(pdata)) {
          rownames(pdata) <- current_rownames
          flowWorkspace::pData(gs) <- pdata
        }
      }, error = function(e) {
        warning("Failed to set pData: ", e$message)
      })
    }
  }
  
  # Add gates and populations
  if (!is.null(gates)) {
    for(idx in seq_len(actual_sample_count)){
      actual_sample_uuids <- sample_uuids[idx] %>% unlist()
      actual_gating_trees <- gating_trees[idx]
      
      add_populations_to_gatingset(gs = gsList[[idx]],
                                   gating_trees = actual_gating_trees,
                                   gates = gates,
                                   sample_uuids = actual_sample_uuids,
                                   strip_comp_prefix = strip_comp_prefix,
                                   verbose = .pkgenv$verbose) # nocov
      
    }
  }
  return(gsList)
}


#' Create Sample Names
#' @keywords internal
create_sample_names <- function(sample_uuids, dataSources, additional.keys, additional.sampleID) {
  sapply(sample_uuids, function(uuid) {
    ds <- dataSources[[uuid]]
    
    # Start with filename
    name_parts <- basename(ds$definition$uri %||% uuid)
    
    # Add additional keywords
    if (!is.null(additional.keys) && length(additional.keys) > 0) {
      for (key in additional.keys) {
        val <- ds$definition$customKeywords[[key]]
        if (!is.null(val)) {
          name_parts <- c(name_parts, as.character(val))
        }
      }
    }
    
    # Add sample ID if requested
    if (additional.sampleID) {
      name_parts <- c(name_parts, uuid)
    }
    
    paste(name_parts, collapse = "_")
  })
}


#' Extract pData from Keywords
#' @keywords internal
extract_pdata <- function(sample_uuids, dataSources, keywords, ignore.case) {
  pdata_list <- lapply(sample_uuids, function(uuid) {
    ds <- dataSources[[uuid]]
    kw <- ds$definition$customKeywords
    
    # Extract requested keywords
    if (ignore.case) {
      kw_names_lower <- tolower(names(kw))
      keywords_lower <- tolower(keywords)
      vals <- lapply(keywords, function(k) {
        idx <- which(kw_names_lower == tolower(k))
        if (length(idx) > 0) kw[[idx[1]]] else NA
      })
      names(vals) <- keywords
    } else {
      vals <- kw[keywords]
    }
    
    as.data.frame(vals, stringsAsFactors = FALSE)
  })
  
  do.call(rbind, pdata_list)
}