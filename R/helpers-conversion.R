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
        keys_df <- extract_keywords_for_samples(
            sample_uuids, dataSources,
            keywords
        )
        subset_parsed <- try(
            {
                filtered <- filter(keys_df, !!enquo(subset))
                filtered$sample_uuid
            },
            silent = TRUE
        )

        if (inherits(subset_parsed, "try-error")) {
            stop("Invalid 'subset' argument: ", attr(
                subset_parsed,
                "condition"
            )$message)
        }
    }

    # Handle different subset types
    if (is.numeric(subset_parsed)) {
        # Numeric indices
        return(sample_uuids[subset_parsed])
    } else if (is.character(subset_parsed)) {
        # Filenames - match to sample UUIDs
        return(filter_sample_filenames(
            subset_parsed, sample_uuids, dataSources
        ))
    } else if (is.list(subset_parsed) && "name" %in% names(subset_parsed)) {
        # List with name element
        return(filter_samples(
            sample_uuids, subset_parsed$name, dataSources,
            keywords
        ))
    }

    # Default: return all
    return(sample_uuids)
}

#' Match filename subsets to sample UUIDs
#'
#' @param subset_parsed Character vector of filenames
#' @param sample_uuids Sample UUIDs to intersect with
#' @param dataSources Data sources searched by basename and File Name keyword
#' @return Sample UUIDs matching any of the filenames
#' @noRd
filter_sample_filenames <- function(subset_parsed, sample_uuids, dataSources) {
    matched <- vapply(subset_parsed, function(fname) {
        idx <- which(vapply(dataSources, function(ds) {
            isTRUE(basename(ds$definition$uri) == fname) ||
                isTRUE(ds$definition$customKeywords$`File Name` == fname)
        }, logical(1)))
        if (length(idx) > 0) names(dataSources)[idx[1]] else NA_character_
    }, character(1))
    matched <- matched[!is.na(matched)]
    intersect(sample_uuids, matched)
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
find_root_population <- function(
    populations, populationDefinitions,
    sample_uuid
) {
    # Find root population definition (type = "root")
    root_popdefs <- Filter(function(pd) {
        !is.null(pd$definition$type) && pd$definition$type == "root"
    }, populationDefinitions)

    if (length(root_popdefs) == 0) {
        stop("Could not find root population definition")
    }

    root_popdef_uuid <- names(root_popdefs)[1]

    # In FlowJo v11, we need to find populations that reference this
    #   population definition
    # and belong to our sample.  First by direct reference, then any
    #   population with
    # parentPopulation = NULL (which indicates a root population).
    # TODO: there is only compoundPopulations but no populationReference, this
    #   needs to be verified
    sample_pop_uuid <- find_root_by_reference(
        populations, root_popdef_uuid, sample_uuid
    )

    if (is.null(sample_pop_uuid)) {
        sample_pop_uuid <- find_root_by_parentless(
            populations, sample_uuid
        )
    }

    if (is.null(sample_pop_uuid)) {
        stop("Could not find root population for sample: ", sample_uuid)
    }

    sample_pop_uuid
}

#' Data source UUID attached to a population, or NULL
#'
#' @param pop Population entry from the workspace
#' @return Data source UUID string or NULL
#' @noRd
root_pop_ds_uuid <- function(pop) {
    parents <- pop$parents
    if (is.null(parents)) {
        return(NULL)
    }
    ds_uuid <- parents[["_dataSource"]]
    if (is.null(ds_uuid)) {
        ds_uuid <- parents[["dataSource"]]
    }
    ds_uuid
}

#' Find the root population referencing the root definition
#'
#' @param populations Populations list from the workspace
#' @param root_popdef_uuid UUID of the root population definition
#' @param sample_uuid Sample the population must belong to
#' @return Population UUID or NULL
#' @noRd
find_root_by_reference <- function(populations, root_popdef_uuid,
    sample_uuid) {
    for (pop_uuid in names(populations)) {
        pop <- populations[[pop_uuid]]
        if (is.null(pop)) {
            next
        }
        # Check if this population references our root population definition
        # and belongs to our sample (through parent relationships)
        if (!is.null(pop$populationReference) &&
            pop$populationReference == root_popdef_uuid) {
            ds_uuid <- root_pop_ds_uuid(pop)
            if (!is.null(ds_uuid) && ds_uuid[[1]] == sample_uuid) {
                return(pop_uuid)
            }
        }
    }
    NULL
}

#' Find a root population by absent parentPopulation
#'
#' @param populations Populations list from the workspace
#' @param sample_uuid Sample the population must belong to
#' @return Population UUID or NULL
#' @noRd
find_root_by_parentless <- function(populations, sample_uuid) {
    for (pop_uuid in names(populations)) {
        pop <- populations[[pop_uuid]]
        if (is.null(pop) || !is.null(pop$parentPopulation)) {
            next
        }
        ds_uuid <- root_pop_ds_uuid(pop)
        if (!is.null(ds_uuid) && ds_uuid[[1]] == sample_uuid) {
            return(pop_uuid)
        }
    }
    NULL
}


#' Identify logical gates and report the summary
#'
#' Runs identify_logical_gates and, when verbose, prints the count and
#' a formatted summary table of the logical gates found.
#'
#' @param populations Populations list from the workspace
#' @param populationDefinitions Population definitions list
#' @return List with \code{gates} (logical gate info) and the updated
#'   \code{populationDefinitions}
#' @keywords internal
tree_identify_logical_gates <- function(populations, populationDefinitions) {
    # First, identify all logical gates before building the tree
    if (.pkgenv$verbose) message("Identifying logical gates...") # nocov
    gate_result <- identify_logical_gates(populations, populationDefinitions)
    logical_gates_info <- gate_result$gates
    populationDefinitions <- gate_result$populationDefinitions
    # Use updated version!

    if (length(logical_gates_info) > 0) {
        if (.pkgenv$verbose) {
            message(sprintf(
                "Found %d logical gates",
                length(logical_gates_info)
            ))
        } # nocov
        if (.pkgenv$verbose) { # nocov
            summary_df <- create_logical_gate_summary(logical_gates_info)
            message(paste0(capture.output(format(summary_df)),
                collapse = "\n"
            ))
        }
    } else {
        if (.pkgenv$verbose) message("No logical gates found") # nocov
    }

    list(
        gates = logical_gates_info,
        populationDefinitions = populationDefinitions
    )
}

#' Post-process the gating tree after construction
#'
#' Moves logical gates up one level, deduplicates the tree, and prints
#' the updated summary when verbose. Only applied when logical gates
#' exist.
#'
#' @param tree The built gating tree
#' @param logical_gates_info Logical gate info list
#' @return The post-processed tree
#' @keywords internal
tree_postprocess <- function(tree, logical_gates_info) {
    if (length(logical_gates_info) == 0) {
        return(tree)
    }

    if (.pkgenv$verbose) {
        message(
            "\nMoving logical gates up in",
            " hierarchy..."
        )
    } # nocov
    tree <- move_logical_gates_up(tree)

    # Remove duplicates at each level (after moving)
    if (.pkgenv$verbose) message("\nRemoving duplicates...") # nocov
    tree <- deduplicate_tree(tree)

    # Get updated summary
    summary_after <- summarize_logical_gates(tree)
    if (!is.null(summary_after)) {
        if (.pkgenv$verbose) {
            message("\nLogical gates in final tree:")
        } # nocov
        if (.pkgenv$verbose) {
            message(paste0(capture.output(format(summary_after)),
                collapse = "\n"
            ))
        } # nocov
    }

    tree
}

#' Check recursion state and resolve the population
#'
#' Guards against circular references and missing populations,
#' then returns the population, its parents, and its definition.
#' Returns NULL when the node must be skipped.
#'
#' @param pop_uuid Population identifier
#' @param visited Environment of visited population UUIDs
#' @param populations Populations list from the workspace
#' @param populationDefinitions Population definitions list
#' @return List with \code{pop}, \code{pop_parents},
#'   \code{pop_def}, and \code{pop_def_uuid}, or NULL when the node
#'   must be skipped
#' @keywords internal
tree_resolve_population <- function(pop_uuid, visited, populations,
                                    populationDefinitions) {
    # Prevent infinite recursion
    if (!is.null(visited[[pop_uuid]])) {
        if (.pkgenv$verbose) {
            warning(
                "Circular reference detected for population: ",
                pop_uuid
            )
        } # nocov
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
    pop_def_uuid <- pop_parents[["populationDefinitions"]]
    if (is.null(pop_def_uuid)) {
        warning(
            "Population has no populationDefinitions parent: ",
            pop_uuid
        )
        return(NULL)
    }

    pop_def <- populationDefinitions[[pop_def_uuid[[1]]]]
    if (is.null(pop_def)) {
        warning("Population definition not found: ", pop_def_uuid)
        return(NULL)
    }

    list(
        pop = pop,
        pop_parents = pop_parents,
        pop_def = pop_def,
        pop_def_uuid = pop_def_uuid
    )
}

#' Construct a tree node from a population definition
#'
#' Sanitizes the node name and builds the base node list.
#'
#' @param pop_uuid Population identifier
#' @param pop_parents Parents list of the population
#' @param pop_def Population definition
#' @param pop_def_uuid Definition UUID
#' @param parent_path Path of the parent node (or NULL)
#' @return The constructed node list
#' @keywords internal
tree_make_node <- function(pop_uuid, pop_parents, pop_def,
                            pop_def_uuid, parent_path) {
    # Extract definition details safely
    node_name <- if (!is.null(pop_def$definition$name)) {
        pop_def$definition$name
    } else "Unnamed"
    # flowWorkspace uses '/' as the path separator, so population names
    #   containing
    # '/' get rewritten (e.g. "CD45, l/d subset" -> "CD45, l:d subset").
    #   Sanitize
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
}

#' Attach logical gate information to a tree node
#'
#' Looks up logical gate info for the population and fills the
#' node's \code{logical_gate_info} and gateDefinition fields.
#'
#' @param node The tree node to modify
#' @param pop_uuid Population identifier
#' @param logical_gates_info Logical gate info list
#' @return The updated node
#' @keywords internal
tree_attach_logical_gate <- function(node, pop_uuid,
                                    logical_gates_info) {
    # Add logical gate information if this is a logical gate
    gate_info <- find_gate_info(pop_uuid, logical_gates_info)
    if (!is.null(gate_info)) {
        if (.pkgenv$verbose) {
            message(sprintf(
                "Adding gate info for: %s",
                unlist(node$name)
            ))
        } # nocov

        node$logical_gate_info <- list(
            operator = gate_info$gate_type,
            combined_populations = gate_info$combined_populations,
            combined_population_uuids = gate_info$combined_population_uuids,
            combined_definition_uuids = gate_info$combined_definition_uuids
        )

        # The gateDefinition should already be in pop_def from
        #   identify_logical_gates
        # But verify and add if missing
        if (is.null(node$pop_def$definition$gateDefinition)) {
            if (.pkgenv$verbose) {
                message(
                    sprintf(
                        "  -> Adding missing gateDefinition for %s",
                        unlist(node$name)
                    )
                )
            } # nocov
            node$pop_def$definition$gateDefinition <- list(
                type = "logical",
                operator = gate_info$gate_type,
                components = gate_info$combined_populations,
                component_uuids = gate_info$combined_population_uuids
            )
        } else {
            if (.pkgenv$verbose) {
                message(
                    sprintf(
                        "  -> gateDefinition already exists for %s",
                        unlist(node$name)
                    )
                )
            } # nocov
        }
    }
    node
}

#' Collect valid child populations of a node
#'
#' Orders children by populationNumber and keeps only those with
#' a data source parent.
#'
#' @param pop The parent population
#' @param populations Populations list from the workspace
#' @return Named list of child populations
#' @keywords internal
tree_collect_children <- function(pop, populations) {
    # Find children populations
    child_populations <- list()

    if (!is.null(pop$children) && !is.null(pop$children$populations)) {
        child_uuids <- pop$children$populations

        if (length(child_uuids) > 0) {
            # get order
            # here we need to use populationNumber if available.
            pop_order <- c()
            idx <- 1
            for (child_pop_uuid in child_uuids) {
                child_pop <- populations[[child_pop_uuid]]
                if (is.null(child_pop$definition$populationNumber)) {
                    pop_order <- c(pop_order, idx)
                } else {
                    pop_order <- c(
                        pop_order,
                        child_pop$definition$populationNumber
                    )
                }
                idx <- idx + 1
            }
            # browser() # nocov
            for (child_pop_uuid in child_uuids[order(pop_order)]) {
                child_pop <- populations[[child_pop_uuid]]
                pop_num <- child_pop$definition$populationNumber
                if (!is.null(child_pop)) {
                    child_sample_uuid <-
                        child_pop$parents[["_dataSource"]] %||%
                        child_pop$parents[["dataSource"]]
                    if (!is.null(child_sample_uuid)) {
                        child_populations[[child_pop_uuid]] <- child_pop
                    }
                }
            }
        }
    }
    child_populations
}

#' Build child nodes for a populated node
#'
#' Recursively builds child nodes via \code{build_node} and
#' attaches them to the parent node when non-empty.
#'
#' @param node The tree node to modify
#' @param child_populations Named list of child populations
#' @param current_path Path of the current node
#' @param pop_uuid Population identifier of the current node
#' @param build_node Recursive node builder function
#' @return The updated node
#' @keywords internal
tree_build_children <- function(node, child_populations,
                                current_path, pop_uuid, build_node) {
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
    node
}

#' Create a recursive tree node builder
#'
#' Returns the \code{build_node} function used by
#' \code{build_gating_tree}, closed over the workspace state.
#'
#' @param populations Populations list from the workspace
#' @param populationDefinitions Population definitions list
#' @param visited Environment of visited population UUIDs
#' @param logical_gates_info Logical gate info list
#' @return The \code{build_node} function
#' @keywords internal
tree_make_builder <- function(populations, populationDefinitions,
                                visited, logical_gates_info) {
    build_node <- function(pop_uuid, parent_path = NULL,
                            parent_pop_uuid = NULL) {
        resolved <- tree_resolve_population(
            pop_uuid, visited, populations, populationDefinitions
        )
        if (is.null(resolved)) {
            return(NULL)
        }
        pop <- resolved$pop
        pop_parents <- resolved$pop_parents
        pop_def <- resolved$pop_def
        pop_def_uuid <- resolved$pop_def_uuid

        node <- tree_make_node(
            pop_uuid, pop_parents, pop_def, pop_def_uuid,
            parent_path
        )
        node <- tree_attach_logical_gate(
            node, pop_uuid, logical_gates_info
        )

        # Build current path. Both parts are already sanitized, but
        #   keep the path separator as '/' (flowWorkspace convention)
        #   while ensuring no stray '/' from names leaks in.
        current_path <- if (is.null(parent_path)) {
            node$name
        } else {
            paste(parent_path, node$name, sep = "/")
        }

        child_populations <- tree_collect_children(pop, populations)
        node <- tree_build_children(
            node, child_populations, current_path, pop_uuid,
            build_node
        )

        return(node)
    }
    build_node
}

#' Build Gating Tree for given Sample
#' @keywords internal
build_gating_tree <- function(
    sample_uuid, populations,
    populationDefinitions, root_uuid,
    schema_version = NULL
) {
    # First, identify all logical gates before building the tree
    if (.pkgenv$verbose) message("Identifying logical gates...") # nocov
    gate_result <- tree_identify_logical_gates(
        populations, populationDefinitions
    )
    logical_gates_info <- gate_result$gates
    populationDefinitions <- gate_result$populationDefinitions
    # Use updated version!


    # Track visited population UUIDs to prevent infinite loops
    visited <- new.env(parent = emptyenv())

    # Recursive function to build tree
    build_node <- tree_make_builder(
        populations, populationDefinitions, visited,
        logical_gates_info
    )

    # Build the tree
    if (.pkgenv$verbose) message("Building tree...") # nocov
    tree <- build_node(root_uuid)

    tree <- tree_postprocess(tree, logical_gates_info)


    # Clean up visited environment
    rm(list = ls(envir = visited), envir = visited)

    return(tree)
}
# Main function to identify logical gates from populations and definitions
#' Find a population definition by UUID
#'
#' @param populationDefinitions Population definitions list
#' @param uuid Definition UUID to find
#' @return Definition entry or NULL
#' @noRd
logical_find_pop_def <- function(populationDefinitions, uuid) {
    for (i in seq_along(populationDefinitions)) {
        if (!is.null(populationDefinitions[[i]]$uuid) &&
            populationDefinitions[[i]]$uuid == uuid) {
            return(populationDefinitions[[i]])
        }
    }
    NULL
}

#' Find a population by UUID
#'
#' @param populations Populations list from the workspace
#' @param uuid Population UUID to find
#' @return Population entry or NULL
#' @noRd
logical_find_pop <- function(populations, uuid) {
    for (i in seq_along(populations)) {
        if (!is.null(populations[[i]]$uuid) &&
            populations[[i]]$uuid == uuid) {
            return(populations[[i]])
        }
    }
    NULL
}

#' Find what populations are combined by a logical gate
#'
#' The populations combined by the logical gate are in parents$populations.
#'
#' @param pop Population entry of the logical gate
#' @param populations Populations list from the workspace
#' @param populationDefinitions Population definitions list
#' @return Named list entries with population_uuid, definition_uuid, name
#' @noRd
logical_find_combined <- function(pop, populations, populationDefinitions) {
    parent_pop_uuids <- unlist(pop$parents$populations)

    if (is.null(parent_pop_uuids) || length(parent_pop_uuids) == 0) {
        return(NULL)
    }

    combined_pops <- list()

    for (parent_uuid in parent_pop_uuids) {
        parent_pop <- logical_find_pop(populations, parent_uuid)

        if (is.null(parent_pop)) {
            next
        }

        # Get the definition for this population to get the name
        pop_def_uuid <- unlist(parent_pop$parents$populationDefinitions)

        if (is.null(pop_def_uuid) || length(pop_def_uuid) == 0) {
            next
        }

        parent_def <- logical_find_pop_def(populationDefinitions,
            pop_def_uuid[1])

        if (!is.null(parent_def) && !is.null(parent_def$definition$name)) {
            combined_pops <- c(combined_pops, list(list(
                population_uuid = parent_pop$uuid,
                definition_uuid = parent_def$uuid,
                name = sanitize_population_name(
                    unlist(parent_def$definition$name)
                )
            )))
        }
    }

    combined_pops
}

#' One logical-gate result entry, adding gateDefinition when missing
#'
#' @param pop Population entry of the logical gate
#' @param pop_def Population definition of the logical gate
#' @param populationDefinitions Population definitions list (mutated when the
#'   definition lacks a gateDefinition)
#' @param combined Result of logical_find_combined()
#' @return list(gate = entry, populationDefinitions = definitions)
#' @noRd
logical_gate_entry <- function(pop, pop_def, populationDefinitions, combined) {
    combined_names <- vapply(
        combined, function(x) unlist(x$name),
        character(1)
    )

    # Add gateDefinition if missing
    if (is.null(pop_def$definition$gateDefinition)) {
        populationDefinitions[[pop_def$uuid]]$definition$
            gateDefinition <- list(
            type = "logical",
            operator = pop_def$definition$type,
            components = combined_names,
            component_uuids = vapply(
                combined,
                function(x) unlist(x$population_uuid), character(1)
            )
        )
    }

    entry <- list(
        population_uuid = pop$uuid,
        definition_uuid = pop_def$uuid,
        gate_name = unlist(pop_def$definition$name),
        gate_type = pop_def$definition$type,
        combined_populations = combined_names,
        combined_population_uuids = vapply(
            combined,
            function(x) unlist(x$population_uuid), character(1)
        ),
        combined_definition_uuids = vapply(
            combined,
            function(x) unlist(x$definition_uuid), character(1)
        ),
        num_components = length(combined)
    )

    list(gate = entry, populationDefinitions = populationDefinitions)
}

#' Is a population a logical gate?
#'
#' @param pop Population entry
#' @param populationDefinitions Population definitions from the workspace
#' @return TRUE when the population resolves to an and/or/not definition
#' @noRd
logical_is_gate_pop <- function(pop, populationDefinitions) {
    pop_def_uuid <- unlist(pop$parents$populationDefinitions)
    if (is.null(pop_def_uuid) || length(pop_def_uuid) == 0) {
        return(FALSE)
    }
    pop_def <- logical_find_pop_def(populationDefinitions, pop_def_uuid[1])
    if (is.null(pop_def)) {
        return(FALSE)
    }
    !is.null(pop_def$definition$type) &&
        pop_def$definition$type %in% c("and", "or", "not")
}

#' Process one logical-gate population
#'
#' Builds the gate entry, or emits diagnostics when no combined populations
#' were found, and returns the possibly-updated populationDefinitions.
#'
#' @param pop Population entry
#' @param populations All populations
#' @param populationDefinitions Population definitions from the workspace
#' @return list(gate = entry or NULL, populationDefinitions = updated list)
#' @noRd
logical_process_pop <- function(pop, populations, populationDefinitions) {
    pop_def_uuid <- unlist(pop$parents$populationDefinitions)[1]
    pop_def <- logical_find_pop_def(populationDefinitions, pop_def_uuid)
    gate_type <- pop_def$definition$type
    gate_name <- unlist(pop_def$definition$name)

    if (.pkgenv$verbose) {
        message(sprintf(
            "Found logical gate: %s (type: %s, uuid: %s)", # nocov
            gate_name, gate_type, pop$uuid
        ))
    }

    # Find combined populations
    combined <- logical_find_combined(
        pop, populations, populationDefinitions
    )

    if (!is.null(combined) && length(combined) > 0) {
        combined_names <- vapply(
            combined, function(x) unlist(x$name),
            character(1)
        )
        if (.pkgenv$verbose) {
            message(sprintf(
                "  - Combines: %s",
                paste(combined_names, collapse = ", ")
            ))
        } # nocov

        entry <- logical_gate_entry(
            pop, pop_def, populationDefinitions, combined
        )
        populationDefinitions <- entry$populationDefinitions

        list(gate = entry$gate, populationDefinitions = populationDefinitions)
    } else {
        message("  - No combined populations found!")
        message(sprintf(
            "  - parents$populations: %s",
            paste(unlist(pop$parents$populations), collapse = ", ")
        ))
        list(gate = NULL, populationDefinitions = populationDefinitions)
    }
}

identify_logical_gates <- function(populations, populationDefinitions) {
    # Process each population to find logical gates
    # First, identify which populations are logical gates
    logical_gate_indices <- which(vapply(populations, function(pop) {
        logical_is_gate_pop(pop, populationDefinitions)
    }, logical(1)))

    # Process logical gates and collect results.  Use an explicit loop so the
    # gateDefinition mutation of populationDefinitions is visible here (a plain
    # <- inside lapply would only bind in the closure, and BiocCheck discourages
    # the <<- needed otherwise).
    results_list <- list()
    for (i in logical_gate_indices) {
        out <- logical_process_pop(
            populations[[i]], populations, populationDefinitions
        )
        populationDefinitions <- out$populationDefinitions

        if (!is.null(out$gate)) {
            results_list[[length(results_list) + 1]] <- out$gate
        }
    }

    # Remove NULL entries
    results <- results_list[!vapply(results_list, is.null, logical(1))]

    # Return both the gate info AND the updated populationDefinitions
    return(list(
        gates = results,
        populationDefinitions = populationDefinitions
    ))
}

#' Indices of first-occurrence children by UUID
#'
#' Non-list children and children without a UUID are always kept.
#'
#' @param children Children list of a tree node
#' @return Logical vector marking children to keep
#' @noRd
tree_dedupe_keep_indices <- function(children) {
    # Use Filter to keep only first occurrence of each UUID
    # For non-list children or children with unique UUIDs, keep them
    seen_uuids <- character()
    keep_indices <- logical(length(children))

    for (i in seq_along(children)) {
        child <- children[[i]]
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

    keep_indices
}

# Remove duplicate children at each level based on UUID
deduplicate_tree <- function(tree) {
    deduplicate_node <- function(node) {
        if (!is.list(node)) {
            return(node)
        }

        # Process children if they exist
        if (!is.null(node$children) && length(node$children) > 0) {
            # Filter to unique children, then recursively deduplicate
            unique_children <- node$children[tree_dedupe_keep_indices(
                node$children
            )]
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
        gate_name = vapply(
            logical_gates, function(x) unlist(x$gate_name),
            character(1)
        ),
        gate_type = vapply(
            logical_gates, function(x) unlist(x$gate_type),
            character(1)
        ),
        population_uuid = vapply(
            logical_gates,
            function(x) unlist(x$population_uuid), character(1)
        ),
        num_components = vapply(
            logical_gates,
            function(x) length(x$combined_populations), integer(1)
        ),
        stringsAsFactors = FALSE
    )

    df$combined_populations <- vapply(logical_gates, function(x) {
        paste(unlist(x$combined_populations), collapse = " | ")
    }, character(1))

    return(df)
}


#' Build the logical-gate summary dataframe
#'
#' @param gates Collected logical gate info entries
#' @return Dataframe or NULL when gates is empty
#' @noRd
logical_summary_df <- function(gates) {
    if (length(gates) == 0) {
        return(NULL)
    }

    df <- data.frame(
        path = vapply(gates, function(g) {
            paste(unlist(g$path),
                collapse = "/"
            )
        }, character(1)),
        name = vapply(gates, function(g) {
            paste(unlist(g$name),
                collapse = "/"
            )
        }, character(1)),
        type = vapply(gates, function(g) unlist(g$type), character(1)),
        num_components = vapply(
            gates,
            function(g) length(g$combined_populations), integer(1)
        ),
        stringsAsFactors = FALSE
    )

    df$combined_populations <- vapply(gates, function(g) {
        if (!is.null(g$combined_populations) &&
            length(g$combined_populations) > 0) {
            paste(unlist(g$combined_populations), collapse = " | ")
        } else {
            NA_character_
        }
    }, character(1))

    df
}

#' Build the traversal path for a tree node
#'
#' Handles vector names by producing all combinations of paths and names.
#'
#' @param node Tree node
#' @param path Parent path vector
#' @return Path vector for this node
#' @noRd
logical_node_path <- function(node, path) {
    if (is.null(node$name)) {
        return(path)
    }

    node_names <- unlist(node$name)
    # Handle when path is empty string(s) or has values
    if (all(path == "")) {
        node_names
    } else {
        # Create all combinations of paths and names
        as.vector(outer(path, node_names, paste, sep = "/"))
    }
}

#' Walk a tree node collecting logical gate info
#'
#' @param node Tree node
#' @param path Current path vector
#' @return List of gate info entries
#' @noRd
logical_collect_gates <- function(node, path = "") {
    gates_list <- list()

    if (!is.list(node)) {
        return(gates_list)
    }

    current_path <- logical_node_path(node, path)

    if (!is.null(node$type) && node$type %in% c("and", "or", "not")) {
        gate_info <- list(
            path = current_path,
            name = unlist(node$name),
            type = node$type,
            uuid = node$uuid
        )

        if (!is.null(node$logical_gate_info)) {
            gate_info$combined_populations <-
                node$logical_gate_info$combined_populations
            gate_info$num_components <-
                length(node$logical_gate_info$combined_populations)
        } else {
            gate_info$num_components <- 0
        }

        gates_list <- list(gate_info)
    }

    if (!is.null(node$children) && is.list(node$children)) {
        child_gates <- lapply(node$children, function(child) {
            logical_collect_gates(child, current_path)
        })
        gates_list <- c(gates_list, unlist(child_gates, recursive = FALSE))
    }

    gates_list
}

# Summarize logical gates in tree
summarize_logical_gates <- function(tree) {
    gates <- logical_collect_gates(tree)
    logical_summary_df(gates)
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
#' Split a processed child's grandchildren into logical gates and the rest
#'
#' Marks the logical grandchildren with moved_from/parent and returns both
#' groups.
#'
#' @param processed_child Child node after process_node()
#' @param current_target Target ancestor name for moved gates
#' @return list(logical, regular) grandchild lists
#' @noRd
logical_collect_grandchildren <- function(processed_child, current_target) {
    child_logical_gates <- list()
    child_regular_children <- list()

    if (!is.null(processed_child$children)) {
        for (grandchild in processed_child$children) {
            if (is.list(grandchild) &&
                !is.null(grandchild$type) &&
                grandchild$type %in% c("and", "or", "not")) {
                grandchild$moved_from <-
                    unlist(processed_child$name)
                grandchild$parent <- current_target
                child_logical_gates <- c(
                    child_logical_gates,
                    list(grandchild)
                )
            } else {
                child_regular_children <- c(
                    child_regular_children,
                    list(grandchild)
                )
            }
        }
    }

    list(logical = child_logical_gates, regular = child_regular_children)
}

#' Prepare a logical child for moving to the target ancestor level
#'
#' Processes the child (so nested logical gates can bubble up), removes any
#' remaining logical children, and tags moved_from/parent.
#'
#' @param child Logical child node
#' @param current_target Target ancestor name
#' @param process_node Recursive node processor from move_logical_gates_up()
#' @return The processed child node
#' @noRd
logical_prep_logical_child <- function(child, current_target, process_node) {
    processed_child <- process_node(child, current_target)

    # Any logical grandchildren of a logical gate belong at the target
    # ancestor level, not under this gate.
    if (!is.null(processed_child$children)) {
        processed_child$children <- Filter(
            function(g) {
                !(is.list(g) && !is.null(g$type) &&
                    g$type %in% c("and", "or", "not"))
            },
            processed_child$children
        )
        if (length(processed_child$children) == 0) {
            processed_child$children <- NULL
        }
    }

    processed_child$moved_from <- unlist(processed_child$parent)
    processed_child$parent <- current_target

    processed_child
}

#' Prepare a non-logical child, pulling its logical grandchildren up
#'
#' @param child Non-logical child node
#' @param current_target Target ancestor name
#' @param process_node Recursive node processor from move_logical_gates_up()
#' @return list(child = processed child without logical grandchildren,
#'   moved = logical grandchildren)
#' @noRd
logical_prep_regular_child <- function(child, current_target, process_node) {
    # Non-logical child: recurse, then pull its logical grandchildren up.
    processed_child <- process_node(child, current_target)

    grandchildren <- logical_collect_grandchildren(processed_child,
        current_target)

    processed_child$children <- grandchildren$regular
    if (length(processed_child$children) == 0) {
        processed_child$children <- NULL
    }

    list(child = processed_child, moved = grandchildren$logical)
}

#' Handle one list child while moving logical gates up
#'
#' @param child Child node (a list)
#' @param current_target Target ancestor name
#' @param node_is_logical Is the parent node itself a logical gate?
#' @param process_node Recursive node processor from move_logical_gates_up()
#' @return list(child = child to keep at this level or NULL,
#'   moved = gates bubbling up)
#' @noRd
logical_handle_child <- function(child, current_target, node_is_logical,
    process_node) {
    is_logical_gate <- !is.null(child$type) &&
        child$type %in% c("and", "or", "not")

    if (!is_logical_gate) {
        prep <- logical_prep_regular_child(child, current_target,
            process_node)
        return(list(child = prep$child, moved = prep$moved))
    }

    processed_child <- logical_prep_logical_child(child, current_target,
        process_node)

    if (.pkgenv$verbose) { # nocov
        message(sprintf(
            "Moving logical gate '%s' from '%s' to '%s'",
            unlist(processed_child$name),
            processed_child$moved_from,
            processed_child$parent
        ))
    }

    if (!node_is_logical) {
        # Keep the logical gate at this non-logical level.
        list(child = processed_child, moved = list())
    } else {
        # Bubble it up further.
        list(child = NULL, moved = list(processed_child))
    }
}

move_logical_gates_up <- function(tree) {
    # Recursive function to process each node.
    # target_ancestor_name is the name of the nearest non-logical ancestor that
    # logical descendants should be attached to.
    process_node <- function(node, target_ancestor_name = NULL) {
        if (!is.list(node) || is.null(node$children) ||
            length(node$children) == 0) {
            return(node)
        }

        node_is_logical <- !is.null(node$type) && node$type %in% c(
            "and",
            "or", "not"
        )
        current_target <- if (!node_is_logical) {
            unlist(node$name)
        } else target_ancestor_name

        children_to_keep <- list()
        gates_to_move_up <- list()

        for (child in node$children) {
            if (!is.list(child)) {
                children_to_keep <- c(children_to_keep, list(child))
                next
            }

            prep <- logical_handle_child(child, current_target,
                node_is_logical, process_node)
            if (!is.null(prep$child)) {
                children_to_keep <- c(children_to_keep, list(prep$child))
            }
            gates_to_move_up <- c(gates_to_move_up, prep$moved)
        }

        node$children <- c(children_to_keep, gates_to_move_up)
        if (length(node$children) == 0) {
            node$children <- NULL
        }

        return(node)
    }

    process_node(tree)
}


get_uuids <- function(tree, uuids = c()) {
    uuids <- c(uuids, tree$uuid)
    for (idx in seq_along(tree$children)) {
        uuids <- get_uuids(tree$children[[idx]], uuids)
    }
    return(uuids)
}


#' Create GatingSet from Components - Transform First Approach
#' @keywords internal
#' @importFrom flowWorkspace GatingSet
#' @importFrom magrittr %>%
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter
#  names
#' when adding populations? Default TRUE. Set to FALSE if compensation has
#  already
#'   been applied and gate names should match the compensated parameter names.
#' Apply compensation to one cytoframe
#'
#' Maps compensation channel names to cytoframe parameters, compensates, and
#' renames the compensated channels to their original "Comp-" prefixed names
#' (FlowJo 10/11 expects compensated channels as "Comp-<channel>").
#'
#' @param cf Cytoframe to compensate
#' @param comp Compensation object for the sample
#' @return Compensated cytoframe
#' @noRd
gs_compensate_frame <- function(cf, comp) {
    # Store original compensation names (with "Comp-" prefix) before mapping
    if (methods::is(comp, "compensation")) {
        comp_prefix_names <- colnames(comp@spillover)
    } else {
        comp_prefix_names <- colnames(comp)
    }

    # Map compensation channel names to cytoframe parameter names
    # This handles cases where flowCore sanitizes names (e.g., "/" -> "_")
    comp_mapped <- map_compensation_names(comp, colnames(cf))
    cf_comp <- compensate(cf, comp_mapped)

    if (!is.null(comp_prefix_names) && length(comp_prefix_names) > 0) {
        # Get current column names and replace compensated ones with
        #   Comp- prefix
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
        cf_comp
    } else {
        cf_comp
    }
}

#' Apply sample transformations to a single-sample GatingSet
#'
#' Maps transformation channel names to flowFrame parameter names, keeps only
#' non-NULL, non-linear transforms (flowWorkspace-compatible objects), and
#' transforms the GatingSet.
#'
#' @param gs_single GatingSet holding the single sample
#' @param sample_transformations Transformations for the sample
#' @param cf Cytoframe whose colnames parameter names map against
#' @param i Sample index used in warning messages
#' @return The transformed GatingSet
#' @noRd
gs_transform_single <- function(gs_single, sample_transformations, cf, i) {
    trans_mapped <- map_transformation_names(
        sample_transformations,
        colnames(cf)
    )

    # Keep only non-NULL, non-linear transforms
    #   (flowWorkspace-compatible objects)
    trans_apply <- Filter(function(t) {
        !is.null(t) &&
            (is.null(attr(t, "type")) || attr(t, "type") != "linear")
    }, trans_mapped)

    if (length(trans_apply) > 0) {
        tryCatch(
            {
                transList <- flowWorkspace::transformerList(
                    from  = names(trans_apply),
                    trans = trans_apply
                )
                gs_single <- flowWorkspace::transform(
                    gs_single,
                    transList
                )
            },
            error = function(e) {
                warning(sprintf(
                    "Failed to apply transformations for sample %d: %s",
                    i, e$message
                ))
            }
        )
    }

    gs_single
}

#' Build one single-sample GatingSet
#'
#' Applies compensation and sample transformations to the sample's cytoframe,
#' then wraps it in a one-sample GatingSet.
#'
#' @param i Sample index (for messages)
#' @param actual_sample_count Total sample count (for messages)
#' @param cf Cytoframe for the sample
#' @param sample_uuid UUID of the sample
#' @param compensations Named compensation list (by sample UUID)
#' @param transformations Named transformation list (by sample UUID)
#' @return One single-sample GatingSet
#' @noRd
gs_build_single_hierarchy <- function(i, actual_sample_count, cf,
    sample_uuid, compensations, transformations) {
    # Apply compensation (if available)
    if (!is.null(compensations[[sample_uuid]])) {
        cf <- gs_compensate_frame(cf, compensations[[sample_uuid]])
    }

    cs <- cytoset()
    cs_add_cytoframe(cs, identifier(cf), cf)
    # Create GatingSet from single sample
    gs_single <- GatingSet(cs)

    # Apply transformations (sample-specific).  Falls back to the first
    #   transformation entry when the sample has none of its own.
    sample_transformations <- transformations[[sample_uuid]]
    if (is.null(sample_transformations) && length(transformations) > 0) {
        sample_transformations <- transformations[[1]]
    }

    # When transform=TRUE, data is transformed to display space and gates
    #   are
    # extracted in the same transformed space (via use_transformed_coords).
    # This keeps gating evaluation, visualization, and gate coordinates
    # consistent.
    if (!is.null(sample_transformations) &&
        length(sample_transformations) > 0) {
        gs_single <- gs_transform_single(
            gs_single, sample_transformations, cf, i
        )
    }

    gs_single
}

#' Build the per-sample GatingHierarchy list
#'
#' @param cytoset Cyotoset with all samples
#' @param sample_uuids Sample UUIDs aligned with cytoset
#' @param compensations Named compensation list (by sample UUID)
#' @param transformations Named transformation list (by sample UUID)
#' @return List of single-sample GatingSet objects
#' @noRd
gs_build_hierarchy_list <- function(cytoset, sample_uuids, compensations,
    transformations) {
    actual_sample_count <- length(cytoset)
    # Create individual GatingHierarchy objects with transformations
    gsList <- list()

    for (i in seq_len(min(length(sample_uuids), actual_sample_count))) {
        sample_uuid <- sample_uuids[i][[1]]

        if (.pkgenv$verbose) {
            message(
                "Processing sample ", i, " of ",
                actual_sample_count
            )
        } # nocov

        # Extract single cytoframe
        cf <- cytoset[[i]]

        gsList[[i]] <- gs_build_single_hierarchy(
            i, actual_sample_count, cf, sample_uuid,
            compensations, transformations
        )
    }

    gsList
}

#' Attach pData (keywords) to the merged GatingSet
#'
#' @param gs Merged GatingSet
#' @param sample_uuids Sample UUIDs actually present in the GatingSet
#' @param dataSources Data sources from the workspace
#' @param keywords Keyword names to expose in pData
#' @param keyword.ignore.case Case-insensitive keyword matching?
#' @param actual_sample_count Number of samples in the GatingSet
#' @return The GatingSet (pData set when dimensions match)
#' @noRd
gs_set_pdata <- function(gs, sample_uuids, dataSources, keywords,
    keyword.ignore.case, actual_sample_count) {
    actual_sample_uuids <- sample_uuids[seq_len(actual_sample_count)]
    pdata <- extract_pdata(
        actual_sample_uuids, dataSources, keywords,
        keyword.ignore.case
    )

    if (!is.data.frame(pdata)) {
        pdata <- as.data.frame(pdata, stringsAsFactors = FALSE)
    }

    if (nrow(pdata) == actual_sample_count) {
        tryCatch(
            {
                current_rownames <- rownames(flowWorkspace::pData(gs))
                if (length(current_rownames) == nrow(pdata)) {
                    rownames(pdata) <- current_rownames
                    flowWorkspace::pData(gs) <- pdata
                }
            },
            error = function(e) {
                warning("Failed to set pData: ", e$message)
            }
        )
    }

    gs
}

#' Add gates and populations to each GatingHierarchy
#'
#' @param gsList List of GatingHierarchy objects
#' @param sample_uuids Sample UUIDs aligned with gsList
#' @param gating_trees Gating trees for the samples
#' @param gates Gates list from the workspace
#' @param strip_comp_prefix Strip "Comp-" prefix from gate parameter names?
#' @param actual_sample_count Number of samples (as in the caller)
#' @noRd
gs_add_gates <- function(gsList, sample_uuids, gating_trees, gates,
    strip_comp_prefix, actual_sample_count, schema_version = NULL) {
    for (idx in seq_len(actual_sample_count)) {
        actual_sample_uuids <- sample_uuids[idx] %>% unlist()
        actual_gating_trees <- gating_trees[idx]

        add_populations_to_gatingset(
            gs = gsList[[idx]],
            gating_trees = actual_gating_trees,
            gates = gates,
            sample_uuids = actual_sample_uuids,
            strip_comp_prefix = strip_comp_prefix,
            verbose = .pkgenv$verbose,
            schema_version = schema_version
        ) # nocov
    }

    gsList
}

#' Create GatingSet from Components - Transform First Approach
#' @keywords internal
#' @importFrom flowWorkspace GatingSet
#' @importFrom magrittr %>%
#' @param strip_comp_prefix Logical. Strip "Comp-" prefix from gate parameter
#  names
#' when adding populations? Default TRUE. Set to FALSE if compensation has
#  already
#'   been applied and gate names should match the compensated parameter names.
#' Merge the hierarchy list, name samples, and attach pData
#'
#' @param gsList List of single-sample GatingSets
#' @param gs Merged GatingSet
#' @param sample_uuids Sample UUIDs aligned with gsList
#' @param dataSources Data sources from the workspace
#' @param keywords Keyword names to expose in pData
#' @param additional.keys Additional keywords to build sample names from
#' @param additional.sampleID Include the UUID in sample names?
#' @param keyword.ignore.case Case-insensitive keyword matching?
#' @return list(gsList = named list, gs = merged GatingSet with pData)
#' @noRd
gs_merge_and_name <- function(gsList, gs, sample_uuids, dataSources,
    keywords, additional.keys, additional.sampleID,
    keyword.ignore.case) {
    # Set sample names
    sample_names <- create_sample_names(
        sample_uuids,
        dataSources,
        additional.keys,
        additional.sampleID
    )
    names(gsList) <- sample_names
    # tryCatch({
    #   flowWorkspace::sampleNames(gs) <- sample_names
    # }, error = function(e) {
    #   warning("Could not set sample names: ", e$message)
    # })

    # Add pData (keywords)
    if (length(keywords) > 0) {
        gs <- gs_set_pdata(
            gs, sample_uuids, dataSources, keywords,
            keyword.ignore.case, length(gsList)
        )
    }

    list(gsList = gsList, gs = gs)
}

create_gatingset_from_cytoset <- function(
    cytoset,
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
    strip_comp_prefix = TRUE,
    schema_version = NULL
) {
    actual_sample_count <- length(cytoset)

    gsList <- gs_build_hierarchy_list(
        cytoset, sample_uuids, compensations, transformations
    )

    # Combine GatingHierarchy objects into a GatingSet
    if (.pkgenv$verbose) {
        message(
            "Combining ", length(gsList),
            " GatingHierarchy objects into GatingSet..."
        ) # nocov
    }
    # this will permanately transformt the data and loose transformation
    #   information
    gs <- merge_list_to_gs(gsList)

    merged <- gs_merge_and_name(
        gsList, gs, sample_uuids, dataSources, keywords,
        additional.keys, additional.sampleID, keyword.ignore.case
    )

    # Add gates and populations
    if (!is.null(gates)) {
        gs_add_gates(
            merged$gsList, sample_uuids, gating_trees, gates,
            strip_comp_prefix, actual_sample_count, schema_version
        )
    }
    merged$gsList
}


#' Create Sample Names
#' @keywords internal
create_sample_names <- function(
    sample_uuids, dataSources, additional.keys,
    additional.sampleID
) {
    vapply(sample_uuids, function(uuid) {
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
    }, character(1))
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
