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

#' @title Gate Conversion Functions for FlowJo v11
#' @name gates
#' @keywords internal
#' @importFrom flowCore rectangleGate polygonGate ellipsoidGate compensation logicleTransform arcsinhTransform logTransform linearTransform
#' @importFrom flowWorkspace booleanFilter
NULL

#' Extract All Gates from FlowJo v11 Workspace
#'
#' Converts FlowJo v11 gate definitions to flowCore gate objects
#'
#' @param populationDefinitions Population definitions from workspace
#' @param sample_uuids Vector of sample UUIDs to extract gates for
#' @param channel.ignore.case Logical. Case-insensitive channel matching?
#' @param extend_val Numeric. Threshold for extending gate coordinates
#' @param extend_to Numeric. Value to extend gates to
#' @return Named list of gates for each population-sample combination
#' @keywords internal
extract_all_gates <- function(populationDefinitions,
                              sample_uuids,
                              channel.ignore.case = FALSE,
                              extend_val = 0,
                              extend_to = -4000,
                              correct_faulty_gate = 0,
                              use_transformed_coords = FALSE) {
  
  gates_list <- list()
  
  for (pop_uuid in names(populationDefinitions)) {
    pop_def <- populationDefinitions[[pop_uuid]]
    
    # Skip if no gate definition
    if (is.null(pop_def$definition)) next
    if(length(pop_def$definition$name) == 1 && pop_def$definition$name == "Ungated") next
    
    gate_def <- pop_def$definition$gateDefinition
    desync_table <- pop_def$definition$desyncTable
    
    # Determine if gate has per-sample variations
    has_desync <- !is.null(desync_table) && length(desync_table) > 0
    
    for (sample_uuid in sample_uuids) {
      # Use desync gate if available, otherwise master gate
      if (has_desync && sample_uuid %in% names(desync_table)) {
        gate_to_use <- desync_table[[sample_uuid]]
      } else {
        gate_to_use <- gate_def
      }
      if (.pkgenv$verbose) message("pop_def: ", pop_def, " ", sample_uuid) # nocov
      # browser() # nocov
      if (is.null(gate_to_use)) next
      
      # Convert to flowCore gate
      gate_obj <- convert_flowjo_gate(
        gate = gate_to_use,
        pop_name = pop_def$definition$name,
        pop_type = pop_def$definition$type,
        channel.ignore.case = channel.ignore.case,
        extend_val = extend_val,
        extend_to = extend_to,
        correct_faulty_gate = correct_faulty_gate,
        use_transformed_coords = use_transformed_coords
      )
      
      if (!is.null(gate_obj)) {
        key <- paste0(pop_uuid, "_", sample_uuid)
        gates_list[[key]] <- gate_obj
      } else {
        warning("Failed to convert gate for population: ", pop_def$definition$name, 
                " (", pop_uuid, "), sample: ", sample_uuid)
      }
    }
  }
  
  return(gates_list)
}

#' Convert FlowJo Gate to flowCore Gate Object
#'
#' @keywords internal
convert_flowjo_gate <- function(gate,
                                pop_name,
                                pop_type,
                                channel.ignore.case = FALSE,
                                extend_val = 0,
                                extend_to = -4000,
                                correct_faulty_gate = 0,
                                use_transformed_coords = FALSE) {
  
  gate_type <- gate$type %||% pop_type
  
  # Infer gate type from structure if needed
  if (is.null(gate_type) || gate_type == "gate") {
    if (!is.null(gate$xVertices) && !is.null(gate$yVertices)) {
      gate_type <- "PolygonGate"
    } else if (!is.null(gate$xMin) || !is.null(gate$x$max) || !is.null(gate$yMin) || !is.null(gate$y$max)) {
      gate_type <- "RectangleGate"
    } else if (!is.null(gate$centerX) || !is.null(gate$centerY)) {
      gate_type <- "EllipsoidGate"
    }
  }
  
  if (.pkgenv$verbose) message("Converting gate: ", gate_type, " - ", pop_name) # nocov
  
  tryCatch({
    switch(gate_type,
           "RectangleGate" = ,
           "rectangle" = convert_rectangle_gate(gate, pop_name, extend_val, extend_to, correct_faulty_gate, use_transformed_coords),

           "PolygonGate" = ,
           "polygon" = convert_polygon_gate(gate, pop_name, extend_val, extend_to, correct_faulty_gate, use_transformed_coords),

           "EllipsoidGate" = ,
           "ellipse" = convert_ellipse_gate(gate, pop_name, extend_val, extend_to, correct_faulty_gate, use_transformed_coords),

           "RangeGate" = ,
           "range" = convert_range_gate(gate, pop_name, extend_val, extend_to, correct_faulty_gate, use_transformed_coords),

           "QuadrantGate" = ,
           "quad" = convert_quadrant_gate(gate, pop_name, extend_val, extend_to, correct_faulty_gate, use_transformed_coords),

           "BooleanGate" = convert_boolean_gate(gate, pop_name),
           
           {
             warning("Unsupported gate type: ", gate_type, " for population: ", pop_name)
             NULL
           })
  }, error = function(e) {
    warning("Failed to convert gate ", gate_type, " for ", pop_name, ": ", e$message)
    NULL
  })
}


#' Transform coordinates from FlowJo display space to raw data space
#'
#' @description
#' FlowJo stores gate coordinates in display space \(0 to gateResolution/vectorLength\).
#' This function converts them to raw data space for use with flowCore gates.
#'
#' Transformation logic:
#' - Linear: scale from \[0, vectorLength\] to \[minRange, maxRange\]
#' - Biex/Logicle: apply inverse transform (display is in transformed space)
#'
#' @param display_coords Numeric vector of coordinates in display space
#' @param transform_spec Transform specification from gate axis
#' @param gate_resolution Gate resolution (overrides vectorLength if provided)
#' @param correct_faulty_gate Fallback maxRange value if maxRange=0
#' @return Numeric vector of coordinates in raw data space
#' @keywords internal
display_to_raw <- function(display_coords, transform_spec, gate_resolution = NULL, correct_faulty_gate = 0, use_transformed_coords = FALSE) {
  
  # Handle NULL or empty input
  if (is.null(display_coords) || length(display_coords) == 0) {
    return(numeric(0))
  }
  
  # Unlist if needed
  if (is.list(display_coords)) {
    display_coords <- unlist(display_coords)
  }
  display_coords <- as.numeric(display_coords)
  
  # If no transform spec, return as-is
  if (is.null(transform_spec)) {
    return(display_coords)
  }
  
  # Get transform type
  trans_type <- transform_spec$transformType %||% "Linear"
  
  # Get vector length (display space range)
  vector_length <- gate_resolution %||% transform_spec$vectorLength %||% 256
  
  # Transform based on type
  if (trans_type == "Linear") {
    # Linear: scale from [0, vectorLength] to [minRange, maxRange]
    min_range <- transform_spec$minRange %||% 0
    max_range <- transform_spec$maxRange %||% 262144
    
    # Handle faulty gates with maxRange=0
    if (max_range == 0 && correct_faulty_gate != 0) {
      max_range <- correct_faulty_gate
    }
    
    if (max_range == 0) {
      stop("Linear transform has maxRange=0. Set correct_faulty_gate parameter or fix workspace.")
    }
    
    # Linear scaling
    raw_coords <- (display_coords / vector_length) * (max_range - min_range) + min_range
    return(raw_coords)
    
  } else if (trans_type == "Biex") {

    if (isTRUE(use_transformed_coords)) {
      # Data will be in transformed space (0-channelRange).
      # FlowJo display coords ARE the transformed coords. Return as-is.
      return(display_coords)
    }

    # Default: convert to raw space (no transform on data)
    trans_spec <- parse_transformation_info(transform_spec)
    trans_obj <- create_biexponential_transform(trans_spec)
    raw_coords <- trans_obj$inverse(display_coords)
    return(raw_coords)

  } else if (trans_type == "Log") {

    decades_offset <- transform_spec$decadesOffset %||% 1
    number_decades <- transform_spec$numberDecades %||% 4
    shift          <- transform_spec$shift         %||% 0
    vector_length  <- gate_resolution %||% transform_spec$vectorLength %||% 256

    if (is.null(vector_length) || length(vector_length) == 0 || vector_length == 0) {
      warning("Log transform has invalid vectorLength (", vector_length, "). Using 256.")
      vector_length <- 256
    }

    if (isTRUE(use_transformed_coords)) {
      # Data will be in transformed (display) space, but the scale applied to the
      # data may differ from the gate's native gateResolution. Rescale the gate
      # display coordinates to the transformation's vectorLength so they line up.
      target_scale <- transform_spec$vectorLength %||% 256
      source_scale <- gate_resolution %||% target_scale
      scale_factor <- target_scale / source_scale
      return(display_coords * scale_factor)
    }

    # FlowJo Log transform: display coords are log-scaled.
    # Inverse: raw = 10^(display * numberDecades / vectorLength + decadesOffset - 1) - shift
    raw_coords <- 10^(display_coords * number_decades / vector_length + decades_offset - 1) - shift
    return(raw_coords)

  } else {
    warning("Unsupported transform type: ", trans_type, ". Returning coordinates as-is.")
    return(display_coords)
  }
}

#' Apply extension to coordinates
#'
#' @param coords Numeric vector of coordinates
#' @param extend_val Threshold value
#' @param extend_to Replacement value
#' @return Extended coordinates
#' @keywords internal
apply_extension <- function(coords, extend_val = 0, extend_to = -4000) {
  if (extend_val == 0 && extend_to == -4000) {
    return(coords)  # No extension requested
  }
  
  coords[!is.infinite(coords) & coords < extend_val] <- extend_to
  return(coords)
}


#' Convert Rectangle Gate
#' @keywords internal
#' @importFrom flowCore rectangleGate
convert_rectangle_gate <- function(gate, pop_name, extend_val, extend_to, correct_faulty_gate = 0, use_transformed_coords = FALSE) {
  # browser() # nocov
  # Extract parameters
  x_param <- gate$xAxis$parameterSpec$name %||% gate$xParameter
  y_param <- gate$yAxis$parameterSpec$name %||% gate$yParameter
  
  if (is.null(x_param)) {
    stop("Rectangle gate missing x parameter for: ", pop_name)
  }
  
  # Get gate resolution
  gate_resolution <- gate$gateResolution %||% gate$resolution
  
  # Transform X coordinates
  x_display <- unlist(gate$xVertices)
  x_raw <- display_to_raw(x_display, gate$xAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
  x_raw <- apply_extension(x_raw, extend_val, extend_to)

  x_min <- min(x_raw)
  x_max <- max(x_raw)

  if (is.null(y_param)) {
    # 1D gate
    gate_obj <- flowCore::rectangleGate(
      filterId = pop_name,
      .gate = matrix(c(x_min, x_max), nrow = 2, ncol = 1,
                     dimnames = list(c("min", "max"), x_param))
    )
  } else {
    # 2D gate
    y_display <- unlist(gate$yVertices)
    y_raw <- display_to_raw(y_display, gate$yAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
    y_raw <- apply_extension(y_raw, extend_val, extend_to)
    
    y_min <- min(y_raw)
    y_max <- max(y_raw)
    
    gate_obj <- flowCore::rectangleGate(
      filterId = pop_name[[1]],
      .gate = matrix(c(x_min, y_min, x_max, y_max), nrow = 2, ncol = 2, byrow = TRUE,
                     dimnames = list(c("min", "max"), c(x_param, y_param)))
    )
  }
  
  return(gate_obj)
}

#' Convert Polygon Gate
#' @keywords internal
#' @importFrom flowCore polygonGate
convert_polygon_gate <- function(gate, pop_name, extend_val, extend_to, correct_faulty_gate = 0, use_transformed_coords = FALSE) {
  
  # Extract parameters
  x_param <- gate$xParameter %||%
    gate$xAxis$parameterSpec$name %||%
    gate$xAxis
  
  y_param <- gate$yParameter %||%
    gate$yAxis$parameterSpec$name %||%
    gate$yAxis
  
  # Validate parameters
  if (is.null(x_param) || is.null(y_param) ||
      (is.character(x_param) && nchar(x_param) == 0) ||
      (is.character(y_param) && nchar(y_param) == 0)) {
    stop("Polygon gate missing valid parameters for: ", pop_name)
  }
  
  # Get vertices
  x_display <- unlist(gate$xVertices)
  y_display <- unlist(gate$yVertices)
  
  if (length(x_display) < 3 || length(y_display) < 3) {
    stop("Polygon must have at least 3 vertices for: ", pop_name)
  }
  
  if (length(x_display) != length(y_display)) {
    stop("X and Y coordinates must have same length for: ", pop_name)
  }
  
  # Get gate resolution
  gate_resolution <- gate$gateResolution %||% gate$resolution
  
  # Transform coordinates to raw data space
  x_raw <- display_to_raw(x_display, gate$xAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
  y_raw <- display_to_raw(y_display, gate$yAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
  
  # Apply extension
  x_raw <- apply_extension(x_raw, extend_val, extend_to)
  y_raw <- apply_extension(y_raw, extend_val, extend_to)
  
  # Ensure polygon is closed
  if (x_raw[1] != x_raw[length(x_raw)] || y_raw[1] != y_raw[length(y_raw)]) {
    x_raw <- c(x_raw, x_raw[1])
    y_raw <- c(y_raw, y_raw[1])
  }
  
  # Create boundary matrix
  boundaries <- matrix(
    c(x_raw, y_raw),
    ncol = 2,
    dimnames = list(NULL, c(x_param, y_param))
  )
  
  # Create polygon gate
  gate_obj <- flowCore::polygonGate(
    filterId = unlist(pop_name)[1],
    .gate = boundaries
  )
  
  return(gate_obj)
}

#' Convert Ellipse Gate
#' @keywords internal
#' @importFrom flowCore ellipsoidGate
convert_ellipse_gate <- function(gate, pop_name, extend_val, extend_to, correct_faulty_gate = 0, use_transformed_coords = FALSE) {
  
  # Extract parameters
  x_param <- gate$xAxis$parameterSpec$name %||% gate$xParameter
  y_param <- gate$yAxis$parameterSpec$name %||% gate$yParameter
  
  # Get vertices (should be 2 for ellipse: major and minor axis endpoints)
  x_display <- unlist(gate$xVertices)
  y_display <- unlist(gate$yVertices)
  
  gate_resolution <- gate$gateResolution %||% gate$resolution
  
  if (length(x_display) == 2 && length(y_display) == 2) {
    # Calculate ellipse parameters in DISPLAY space first
    center_x_display <- mean(x_display)
    center_y_display <- mean(y_display)
    
    a_display <- abs(diff(x_display)) / 2  # semi-major axis
    b_display <- abs(diff(y_display)) / 2  # semi-minor axis
    
    # Get rotation angle. FlowJo stores rotationAngle in degrees.
    angle <- gate$rotationAngle %||% 0
    angle_rad <- angle * pi / 180
    
    # Build covariance matrix in display space
    cos_a <- cos(angle_rad)
    sin_a <- sin(angle_rad)
    
    cov_display <- matrix(c(
      a_display^2 * sin_a^2 + b_display^2 * cos_a^2,
      (a_display^2 - b_display^2) * sin_a * cos_a,
      (a_display^2 - b_display^2) * sin_a * cos_a,
      a_display^2 * cos_a^2 + b_display^2 * sin_a^2
    ), nrow = 2, ncol = 2)
    
    # Transform center to raw data space
    center_x_raw <- display_to_raw(center_x_display, gate$xAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
    center_y_raw <- display_to_raw(center_y_display, gate$yAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
    
    # Calculate scale factors for covariance transformation
    # This depends on the transform type
    vector_length <- gate_resolution %||% gate$xAxis$transform$vectorLength %||% 256
    
    # For Linear transform: scale from display to data
    # For Biex: we need the derivative of the inverse transform at the center point
    
    x_trans_type <- gate$xAxis$transform$transformType %||% "Linear"
    y_trans_type <- gate$yAxis$transform$transformType %||% "Linear"
    
    if (x_trans_type == "Linear") {
      max_range_x <- gate$xAxis$transform$maxRange %||% 262144
      if (max_range_x == 0 && correct_faulty_gate != 0) max_range_x <- correct_faulty_gate
      scale_x <- max_range_x / vector_length
    } else {
      # For Biex, scaling is approximately 1 near the center (simplified)
      scale_x <- 1
    }
    
    if (y_trans_type == "Linear") {
      max_range_y <- gate$yAxis$transform$maxRange %||% 262144
      if (max_range_y == 0 && correct_faulty_gate != 0) max_range_y <- correct_faulty_gate
      scale_y <- max_range_y / vector_length
    } else {
      scale_y <- 1
    }
    
    # Scale covariance matrix: Cov_raw = S * Cov_display * S^T
    scale_matrix <- diag(c(scale_x, scale_y))
    cov_raw <- scale_matrix %*% cov_display %*% t(scale_matrix)
    
  } else {
    stop("Ellipse gate has unexpected number of vertices for: ", pop_name)
  }
  
  # Set names for covariance matrix
  colnames(cov_raw) <- c(x_param, y_param)
  rownames(cov_raw) <- c(x_param, y_param)
  
  # Get distance parameter
  distance <- gate$distance %||% gate$radius %||% 1
  
  # Create ellipsoid gate
  gate_obj <- flowCore::ellipsoidGate(
    filterId = pop_name[[1]],
    .gate = cov_raw,
    mean = c(center_x_raw, center_y_raw),
    distance = distance
  )
  
  flowCore::parameters(gate_obj) <- c(x_param, y_param)
  
  return(gate_obj)
}

#' Convert Range Gate (1D)
#' @keywords internal
#' @importFrom flowCore rectangleGate
convert_range_gate <- function(gate, pop_name, extend_val, extend_to, correct_faulty_gate = 0, use_transformed_coords = FALSE) {
  
  # Extract parameter
  param <- gate$parameter %||% 
    gate$xParameter %||%
    gate$xAxis$parameterSpec$name %||%
    gate$xAxis %||% 
    NULL
  
  if (is.null(param) || (is.character(param) && nchar(param) == 0)) {
    stop("Range gate missing valid parameter for: ", pop_name)
  }
  
  # Get vertices
  x_display <- gate$xVertices %||% 
    c(gate$min, gate$max) %||%
    c(gate$xMin, gate$xMax) %||%
    gate$vertices
  
  if (is.null(x_display) || length(x_display) == 0) {
    stop("Range gate missing vertices for: ", pop_name)
  }
  
  # Ensure numeric and unlisted - IMPORTANT ORDER
  if (is.list(x_display)) {
    x_display <- unlist(x_display)
  }
  x_display <- as.numeric(x_display)
  
  # Get gate resolution
  gate_resolution <- gate$gateResolution %||% gate$resolution
  
  # Get transform
  transform_spec <- gate$xAxis$transform %||% gate$transform %||% gate$axis$transform
  
  # Transform to raw space
  x_raw <- display_to_raw(x_display, transform_spec, gate_resolution, correct_faulty_gate, use_transformed_coords)
  
  # Get min/max BEFORE extension
  min_val <- min(x_raw, na.rm = TRUE)
  max_val <- max(x_raw, na.rm = TRUE)
  
  # Apply extension to the min/max values
  if (!is.infinite(min_val) && min_val < extend_val) {
    min_val <- extend_to
  }
  if (!is.infinite(max_val) && max_val < extend_val) {
    max_val <- extend_to
  }
  
  # Create rectangle gate for 1D range
  gate_obj <- flowCore::rectangleGate(
    filterId = pop_name[[1]],
    .gate = matrix(c(min_val, max_val), nrow = 2, ncol = 1,
                   dimnames = list(c("min", "max"), param))
  )
  # message("convert_range_gate: ", pop_name, "  ", min_val, "  ", max_val)
  
  return(gate_obj)
}
#' Convert Quadrant Gate
#' @keywords internal
#' @importFrom flowCore quadGate
convert_quadrant_gate <- function(gate, pop_name, extend_val, extend_to, correct_faulty_gate = 0, use_transformed_coords = FALSE) {
  
  # Extract parameters
  x_param <- gate$xAxis$parameterSpec$name %||% gate$xParameter %||% gate$xAxis
  y_param <- gate$yAxis$parameterSpec$name %||% gate$yParameter %||% gate$yAxis
  
  if (is.null(x_param) || is.null(y_param) ||
      (is.character(x_param) && nchar(x_param) == 0) ||
      (is.character(y_param) && nchar(y_param) == 0)) {
    stop("Quadrant gate missing valid parameters for: ", paste(pop_name, collapse = ", "))
  }
  
  # Extract divider position
  x_div_display <- gate$xDivider %||% gate$divider$x %||%
    (if (!is.null(gate$xVertices)) unlist(gate$xVertices)[[1]] else NULL) %||%
    gate$x %||% 0
  
  y_div_display <- gate$yDivider %||% gate$divider$y %||%
    (if (!is.null(gate$yVertices)) unlist(gate$yVertices)[[1]] else NULL) %||%
    gate$y %||% 0
  
  x_div_display <- as.numeric(x_div_display)
  y_div_display <- as.numeric(y_div_display)
  
  # Get gate resolution
  gate_resolution <- gate$gateResolution %||% gate$resolution
  
  # Transform dividers to raw space
  x_div_raw <- display_to_raw(x_div_display, gate$xAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
  y_div_raw <- display_to_raw(y_div_display, gate$yAxis$transform, gate_resolution, correct_faulty_gate, use_transformed_coords)
  
  # Apply extension
  x_div_raw <- apply_extension(x_div_raw, extend_val, extend_to)
  y_div_raw <- apply_extension(y_div_raw, extend_val, extend_to)
  
  # Validate population names
  if (length(pop_name) != 4) {
    stop("Quadrant gate must define exactly 4 populations, got: ", length(pop_name))
  }
  
  # Create boundary
  boundary <- c(x_div_raw, y_div_raw)
  names(boundary) <- c(x_param, y_param)
  
  # Create quad gate
  base_name <- paste0("quad_", x_param, "_", y_param)
  
  gate_obj <- flowCore::quadGate(
    filterId = base_name,
    .gate = boundary
  )
  
  # Store population names
  attr(gate_obj, "pop_names") <- unlist(pop_name) %>% rev()
  
  return(gate_obj)
}

#' Convert Boolean Gate
#' @keywords internal
#' @importFrom flowWorkspace booleanFilter
#' @importFrom magrittr %>%
convert_boolean_gate <- function(gate, pop_name) {
  
  # Boolean gates reference other populations
  specification <- gate$specification %||%
    gate$definition %||%
    gate$expression %||%
    gate$booleanDefinition %||%
    gate$gateDefinition
  
  if (is.null(specification) || specification == "") {
    warning("Boolean gate missing specification for: ", pop_name)
    return(NULL)
  }
  
  if (length(specification) > 1) {
    specification <- specification[1]
    warning("Boolean gate specification had multiple values for: ", pop_name, ". Using first value.")
  }
  
  # Parse boolean expression
  expr <- tryCatch(
    parse_boolean_expression(specification),
    error = function(e) {
      warning("Failed to parse boolean expression '", specification,
              "' for gate: ", pop_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(expr)) {
    return(NULL)
  }
  
  gate_obj <- tryCatch(
    flowWorkspace::booleanFilter(
      expr = expr,
      filterId = pop_name[[1]]
    ),
    error = function(e) {
      warning("Failed to create boolean filter for: ", pop_name,
              ": ", e$message)
      return(NULL)
    }
  )
  
  return(gate_obj)
}

#' Parse Boolean Expression
#' @keywords internal
parse_boolean_expression <- function(spec) {
  
  if (is.null(spec)) {
    stop("Boolean expression specification is NULL")
  }
  
  if (length(spec) > 1) {
    spec <- spec[1]
  }
  
  if (!is.character(spec)) {
    spec <- as.character(spec)
  }
  
  if (spec == "") {
    stop("Boolean expression specification is empty")
  }
  
  # Replace FlowJo operators with R operators
  expr_string <- spec
  expr_string <- gsub("&", " & ", expr_string)
  expr_string <- gsub("\\|", " | ", expr_string)
  expr_string <- gsub("!", "!", expr_string)
  expr_string <- gsub("\\s+", " ", expr_string)
  expr_string <- trimws(expr_string)
  
  if (expr_string == "") {
    stop("Boolean expression is empty after cleaning")
  }
  
  # Parse as expression
  expr <- tryCatch(
    parse(text = expr_string),
    error = function(e) {
      stop("Failed to parse boolean expression '", expr_string, "': ", e$message)
    }
  )
  
  return(expr)
}
