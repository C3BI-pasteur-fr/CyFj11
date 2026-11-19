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
#' @importFrom flowCore rectangleGate polygonGate ellipsoidGate compensation
#  logicleTransform arcsinhTransform logTransform linearTransform
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
extract_all_gates <- function(
    populationDefinitions,
    sample_uuids,
    channel.ignore.case = FALSE,
    extend_val = 0,
    extend_to = -4000,
    correct_faulty_gate = 0,
    use_transformed_coords = FALSE,
    schema_version = NULL
) {
    gates_list <- list()

    for (pop_uuid in names(populationDefinitions)) {
        pop_def <- populationDefinitions[[pop_uuid]]

        # Skip if no gate definition
        if (is.null(pop_def$definition)) next
        if (length(pop_def$definition$name) == 1 &&
            pop_def$definition$name == "Ungated") {
            next
        }

        gates_list <- gates_for_population_samples(
            gates_list, pop_uuid, pop_def, sample_uuids,
            channel.ignore.case, extend_val, extend_to,
            correct_faulty_gate, use_transformed_coords, schema_version
        )
    }

    return(gates_list)
}

#' Convert gates for one population across all samples
#'
#' @param gates_list Accumulated gates list (updated and returned)
#' @param pop_uuid Population UUID being processed
#' @param pop_def Population definition entry
#' @param sample_uuids Vector of sample UUIDs
#' @param channel.ignore.case Passed to convert_flowjo_gate
#' @param extend_val Passed to convert_flowjo_gate
#' @param extend_to Passed to convert_flowjo_gate
#' @param correct_faulty_gate Passed to convert_flowjo_gate
#' @param use_transformed_coords Passed to convert_flowjo_gate
#' @param schema_version Passed to convert_flowjo_gate
#' @return Updated gates list
#' @noRd
gates_for_population_samples <- function(gates_list, pop_uuid, pop_def,
    sample_uuids, channel.ignore.case, extend_val, extend_to,
    correct_faulty_gate, use_transformed_coords, schema_version = NULL) {
    gate_def <- pop_def$definition$gateDefinition
    desync_table <- pop_def$definition$desyncTable

    # Determine if gate has per-sample variations
    has_desync <- !is.null(desync_table) && length(desync_table) > 0

    for (sample_uuid in sample_uuids) {
        gate_to_use <- gate_def
        if (has_desync && sample_uuid %in% names(desync_table)) {
            gate_to_use <- desync_table[[sample_uuid]]
        }
        if (.pkgenv$verbose) {
            message("pop_def: ", pop_def, " ", sample_uuid)
        } # nocov
        if (is.null(gate_to_use)) next

        gate_obj <- convert_flowjo_gate(
            gate = gate_to_use,
            pop_name = pop_def$definition$name,
            pop_type = pop_def$definition$type,
            channel.ignore.case = channel.ignore.case,
            extend_val = extend_val,
            extend_to = extend_to,
            correct_faulty_gate = correct_faulty_gate,
            use_transformed_coords = use_transformed_coords,
            schema_version = schema_version
        )

        if (!is.null(gate_obj)) {
            gates_list[[paste0(pop_uuid, "_", sample_uuid)]] <- gate_obj
        } else {
            warning(
                "Failed to convert gate for population: ",
                pop_def$definition$name,
                " (", pop_uuid, "), sample: ", sample_uuid
            )
        }
    }

    gates_list
}

#' Convert FlowJo Gate to flowCore Gate Object
#'
#' @keywords internal
convert_flowjo_gate <- function(
    gate,
    pop_name,
    pop_type,
    channel.ignore.case = FALSE,
    extend_val = 0,
    extend_to = -4000,
    correct_faulty_gate = 0,
    use_transformed_coords = FALSE,
    schema_version = NULL
) {
    gate_type <- gate_infer_type(gate, pop_type)

    if (.pkgenv$verbose) message("Converting gate: ", gate_type, " - ",
        pop_name) # nocov

    tryCatch(
        {
            gate_dispatch_type(gate_type, gate, pop_name, extend_val,
                extend_to, correct_faulty_gate, use_transformed_coords,
                schema_version)
        },
        error = function(e) {
            warning(
                "Failed to convert gate ", gate_type, " for ", pop_name,
                ": ", e$message
            )
            NULL
        }
    )
}

#' Infer the FlowJo gate type from type field or structure
#'
#' @param gate Raw gate definition
#' @param pop_type Population type used when the gate has no type of its own
#' @return Gate type string used by gate_dispatch_type()
#' @noRd
gate_infer_type <- function(gate, pop_type) {
    gate_type <- gate$type %||% pop_type

    # Infer gate type from structure if needed
    if (is.null(gate_type) ||
        gate_type == "gate") {
        if (!is.null(gate$xVertices) && !is.null(gate$yVertices)) {
            gate_type <- "PolygonGate"
        } else if (!is.null(gate$xMin) || !is.null(gate$x$max) ||
            !is.null(gate$yMin) || !is.null(gate$y$max)) {
            gate_type <- "RectangleGate"
        } else if (!is.null(gate$centerX) || !is.null(gate$centerY)) {
            gate_type <- "EllipsoidGate"
        }
    }

    gate_type
}

#' Dispatch a typed gate to its flowCore converter
#'
#' @param gate_type Gate type string from gate_infer_type()
#' @param gate Raw gate definition
#' @param pop_name Population name for the flowCore gate
#' @param extend_val Passed through to the converters
#' @param extend_to Passed through to the converters
#' @param correct_faulty_gate Passed through to the converters
#' @param use_transformed_coords Passed through to the converters
#' @return flowCore gate object, or NULL for unsupported types
#' @noRd
gate_dispatch_type <- function(gate_type, gate, pop_name, extend_val,
    extend_to, correct_faulty_gate, use_transformed_coords,
    schema_version = NULL) {
    switch(gate_type,
        "RectangleGate" = ,
        "rectangle" = convert_rectangle_gate(
            gate, pop_name,
            extend_val, extend_to, correct_faulty_gate,
            use_transformed_coords, schema_version
        ),
        "PolygonGate" = ,
        "polygon" = convert_polygon_gate(
            gate, pop_name, extend_val,
            extend_to, correct_faulty_gate, use_transformed_coords,
            schema_version
        ),
        "EllipsoidGate" = ,
        "ellipse" = convert_ellipse_gate(
            gate, pop_name, extend_val,
            extend_to, correct_faulty_gate, use_transformed_coords,
            schema_version
        ),
        "RangeGate" = ,
        "range" = convert_range_gate(
            gate, pop_name, extend_val,
            extend_to, correct_faulty_gate, use_transformed_coords,
            schema_version
        ),
        "QuadrantGate" = ,
        "quad" = convert_quadrant_gate(
            gate, pop_name, extend_val,
            extend_to, correct_faulty_gate, use_transformed_coords,
            schema_version
        ),
        "BooleanGate" = convert_boolean_gate(gate, pop_name, schema_version),
        {
            warning(
                "Unsupported gate type: ", gate_type,
                " for population: ", pop_name
            )
            NULL
        }
    )
}


#' Transform coordinates from FlowJo display space to raw data space
#'
#' @description
#' FlowJo stores gate coordinates in display space \(0 to
#  gateResolution/vectorLength\).
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
display_to_raw <- function(
    display_coords, transform_spec,
    gate_resolution = NULL, correct_faulty_gate = 0,
    use_transformed_coords = FALSE
) {
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

    # Transform based on type
    trans_type <- transform_spec$transformType %||% "Linear"
    if (trans_type == "Linear") {
        display_to_raw_linear(
            display_coords, transform_spec,
            gate_resolution, correct_faulty_gate
        )
    } else if (trans_type == "Biex") {
        display_to_raw_biex(display_coords, transform_spec,
            use_transformed_coords)
    } else if (trans_type == "Log") {
        display_to_raw_log(
            display_coords, transform_spec,
            gate_resolution, use_transformed_coords
        )
    } else {
        warning(
            "Unsupported transform type: ", trans_type,
            ". Returning coordinates as-is."
        )
        display_coords
    }
}

#' Convert display coordinates through a Linear transform
#'
#' Scales from \[0, vectorLength\] to \[minRange, maxRange\].
#'
#' @param display_coords Numeric coordinates in display space
#' @param transform_spec Transform specification from the gate axis
#' @param gate_resolution Gate resolution (overrides vectorLength)
#' @param correct_faulty_gate Fallback maxRange value if maxRange=0
#' @return Numeric vector of coordinates in raw data space
#' @noRd
display_to_raw_linear <- function(display_coords, transform_spec,
    gate_resolution, correct_faulty_gate) {
    min_range <- transform_spec$minRange %||% 0
    max_range <- transform_spec$maxRange %||% 262144

    # Handle faulty gates with maxRange=0
    if (max_range == 0 && correct_faulty_gate != 0) {
        max_range <- correct_faulty_gate
    }

    if (max_range == 0) {
        stop(
            "Linear transform has maxRange=0. Set correct_faulty_gate",
            " parameter or fix workspace."
        )
    }

    # Get vector length (display space range)
    vector_length <- gate_resolution %||%
        transform_spec$vectorLength %||% 256

    # Linear scaling
    (display_coords / vector_length) * (max_range - min_range) + min_range
}

#' Convert display coordinates through a Biex transform
#'
#' @param display_coords Numeric coordinates in display space
#' @param transform_spec Transform specification from the gate axis
#' @param use_transformed_coords Keep coordinates in transformed space?
#' @return Numeric vector of coordinates in raw data space
#' @noRd
display_to_raw_biex <- function(display_coords, transform_spec,
    use_transformed_coords) {
    if (isTRUE(use_transformed_coords)) {
        # Data will be in transformed space (0-channelRange).
        # FlowJo display coords ARE the transformed coords. Return as-is.
        return(display_coords)
    }

    # Default: convert to raw space (no transform on data)
    trans_spec <- parse_transformation_info(transform_spec)
    trans_obj <- create_biexponential_transform(trans_spec)
    trans_obj$inverse(display_coords)
}

#' Convert display coordinates through a Log transform
#'
#' @param display_coords Numeric coordinates in display space
#' @param transform_spec Transform specification from the gate axis
#' @param gate_resolution Gate resolution (overrides vectorLength)
#' @param use_transformed_coords Keep coordinates in transformed space?
#' @return Numeric vector of coordinates in raw data space
#' @noRd
display_to_raw_log <- function(display_coords, transform_spec,
    gate_resolution, use_transformed_coords) {
    decades_offset <- transform_spec$decadesOffset %||% 1
    number_decades <- transform_spec$numberDecades %||% 4
    shift <- transform_spec$shift %||% 0
    vector_length <- gate_resolution %||%
        transform_spec$vectorLength %||% 256

    if (is.null(vector_length) || length(vector_length) == 0 ||
        vector_length == 0) {
        warning(
            "Log transform has invalid vectorLength (",
            vector_length, "). Using 256."
        )
        vector_length <- 256
    }

    if (isTRUE(use_transformed_coords)) {
        # Data will be in transformed (display) space, but the scale
        #   applied to the
        # data may differ from the gate's native gateResolution. Rescale
        #   the gate
        # display coordinates to the transformation's vectorLength so they
        #   line up.
        target_scale <- transform_spec$vectorLength %||% 256
        source_scale <- gate_resolution %||% target_scale
        scale_factor <- target_scale / source_scale
        return(display_coords * scale_factor)
    }

    # FlowJo Log transform: display coords are log-scaled.
    # Inverse: raw = 10^(display * numberDecades / vectorLength +
    #   decadesOffset - 1) - shift
    10^(display_coords * number_decades / vector_length +
        decades_offset - 1) - shift
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
        return(coords) # No extension requested
    }

    coords[!is.infinite(coords) & coords < extend_val] <- extend_to
    return(coords)
}


#' Convert Rectangle Gate
#' @keywords internal
#' @importFrom flowCore rectangleGate
convert_rectangle_gate <- function(
    gate, pop_name, extend_val, extend_to,
    correct_faulty_gate = 0, use_transformed_coords = FALSE,
    schema_version = NULL
) {
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
    x_raw <- display_to_raw(
        x_display, gate$xAxis$transform,
        gate_resolution, correct_faulty_gate, use_transformed_coords
    )
    x_raw <- apply_extension(x_raw, extend_val, extend_to)

    x_min <- min(x_raw)
    x_max <- max(x_raw)

    if (is.null(y_param)) {
        # 1D gate
        gate_obj <- flowCore::rectangleGate(
            filterId = pop_name,
            .gate = matrix(c(x_min, x_max),
                nrow = 2, ncol = 1,
                dimnames = list(c("min", "max"), x_param)
            )
        )
    } else {
        gate_obj <- gate_rect_2d(
            gate, pop_name, x_param, y_param, x_min, x_max,
            gate_resolution, extend_val, extend_to,
            correct_faulty_gate, use_transformed_coords
        )
    }

    return(gate_obj)
}

#' Build a 2D rectangle gate from already-converted X bounds
#'
#' @param gate Raw gate definition
#' @param pop_name Population name
#' @param x_param X channel name
#' @param y_param Y channel name
#' @param x_min Converted X minimum
#' @param x_max Converted X maximum
#' @param gate_resolution Gate resolution
#' @param extend_val Passed to display_to_raw/apply_extension
#' @param extend_to Passed to display_to_raw/apply_extension
#' @param correct_faulty_gate Passed to display_to_raw
#' @param use_transformed_coords Passed to display_to_raw
#' @return flowCore rectangleGate object
#' @noRd
gate_rect_2d <- function(gate, pop_name, x_param, y_param, x_min, x_max,
    gate_resolution, extend_val, extend_to, correct_faulty_gate,
    use_transformed_coords) {
    # 2D gate
    y_display <- unlist(gate$yVertices)
    y_raw <- display_to_raw(
        y_display, gate$yAxis$transform,
        gate_resolution, correct_faulty_gate, use_transformed_coords
    )
    y_raw <- apply_extension(y_raw, extend_val, extend_to)

    y_min <- min(y_raw)
    y_max <- max(y_raw)

    flowCore::rectangleGate(
        filterId = pop_name[[1]],
        .gate = matrix(c(x_min, y_min, x_max, y_max),
            nrow = 2, ncol = 2, byrow = TRUE,
            dimnames = list(c("min", "max"), c(x_param, y_param))
        )
    )
}

#' Convert Polygon Gate
#' @keywords internal
#' @importFrom flowCore polygonGate
convert_polygon_gate <- function(
    gate, pop_name, extend_val, extend_to,
    correct_faulty_gate = 0, use_transformed_coords = FALSE,
    schema_version = NULL
) {
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

    coords_raw <- gate_polygon_vertices(
        gate, pop_name, extend_val, extend_to,
        correct_faulty_gate, use_transformed_coords
    )

    # Create boundary matrix
    boundaries <- matrix(
        c(coords_raw$x, coords_raw$y),
        ncol = 2,
        dimnames = list(NULL, c(x_param, y_param))
    )

    # Create polygon gate
    flowCore::polygonGate(
        filterId = unlist(pop_name)[1],
        .gate = boundaries
    )
}

#' Validate polygon parameters and convert vertices to raw space
#'
#' Stops on invalid parameters or vertex counts; closes the polygon.
#'
#' @param gate Raw gate definition
#' @param pop_name Population name
#' @param extend_val Passed to display_to_raw/apply_extension
#' @param extend_to Passed to display_to_raw/apply_extension
#' @param correct_faulty_gate Passed to display_to_raw
#' @param use_transformed_coords Passed to display_to_raw
#' @return list(x, y) with converted, extended, closed coordinates
#' @noRd
gate_polygon_vertices <- function(gate, pop_name, extend_val, extend_to,
    correct_faulty_gate, use_transformed_coords) {
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
    x_raw <- display_to_raw(
        x_display, gate$xAxis$transform,
        gate_resolution, correct_faulty_gate, use_transformed_coords
    )
    y_raw <- display_to_raw(
        y_display, gate$yAxis$transform,
        gate_resolution, correct_faulty_gate, use_transformed_coords
    )

    # Apply extension
    x_raw <- apply_extension(x_raw, extend_val, extend_to)
    y_raw <- apply_extension(y_raw, extend_val, extend_to)

    # Ensure polygon is closed
    if (x_raw[1] != x_raw[length(x_raw)] || y_raw[1] != y_raw[length(y_raw)]) {
        x_raw <- c(x_raw, x_raw[1])
        y_raw <- c(y_raw, y_raw[1])
    }

    list(x = x_raw, y = y_raw)
}

#' Convert Ellipse Gate
#' @keywords internal
#' @importFrom flowCore ellipsoidGate
convert_ellipse_gate <- function(
    gate, pop_name, extend_val, extend_to,
    correct_faulty_gate = 0, use_transformed_coords = FALSE,
    schema_version = NULL
) {
    # Extract parameters
    x_param <- gate$xAxis$parameterSpec$name %||% gate$xParameter
    y_param <- gate$yAxis$parameterSpec$name %||% gate$yParameter

    # Get vertices (should be 2 for ellipse: major and minor axis endpoints)
    x_display <- unlist(gate$xVertices)
    y_display <- unlist(gate$yVertices)

    gate_resolution <- gate$gateResolution %||% gate$resolution

    if (length(x_display) != 2 || length(y_display) != 2) {
        stop("Ellipse gate has unexpected number of vertices for: ", pop_name)
    }

    ell <- gate_ellipse_geometry(
        gate, x_display, y_display, gate_resolution,
        correct_faulty_gate, use_transformed_coords, schema_version
    )

    # Set names for covariance matrix
    colnames(ell$cov_raw) <- c(x_param, y_param)
    rownames(ell$cov_raw) <- c(x_param, y_param)

    # Get distance parameter
    distance <- gate$distance %||% gate$radius %||% 1

    # Create ellipsoid gate
    gate_obj <- flowCore::ellipsoidGate(
        filterId = pop_name[[1]],
        .gate = ell$cov_raw,
        mean = c(ell$center_x_raw, ell$center_y_raw),
        distance = distance
    )

    flowCore::parameters(gate_obj) <- c(x_param, y_param)

    gate_obj
}

#' Compute ellipse geometry in raw data space
#'
#' Builds the display-space covariance, transforms the center, and scales the
#' covariance: Cov_raw = S * Cov_display * S^T.
#'
#' @param gate Raw gate definition
#' @param x_display X vertex pair in display space
#' @param y_display Y vertex pair in display space
#' @param gate_resolution Gate resolution
#' @param correct_faulty_gate Fallback maxRange value if maxRange=0
#' @param use_transformed_coords Passed to display_to_raw
#' @param schema_version Analysis schema version; selects the ellipse convention
#' @return list(cov_raw, center_x_raw, center_y_raw)
#' @noRd
gate_ellipse_geometry <- function(gate, x_display, y_display, gate_resolution,
    correct_faulty_gate, use_transformed_coords, schema_version = NULL) {
    # Calculate ellipse parameters in DISPLAY space first
    center_x_display <- mean(x_display)
    center_y_display <- mean(y_display)

    a_display <- abs(diff(x_display)) / 2 # semi-major axis
    b_display <- abs(diff(y_display)) / 2 # semi-minor axis

    # Get rotation angle. FlowJo stores rotationAngle in degrees.
    angle <- gate$rotationAngle %||% 0
    angle_rad <- angle * pi / 180
    cov_display <- gate_ellipse_cov_display(
        a_display, b_display, angle_rad, schema_version
    )

    # Transform center to raw data space
    center_x_raw <- display_to_raw(
        center_x_display,
        gate$xAxis$transform, gate_resolution, correct_faulty_gate,
        use_transformed_coords
    )
    center_y_raw <- display_to_raw(
        center_y_display,
        gate$yAxis$transform, gate_resolution, correct_faulty_gate,
        use_transformed_coords
    )

    # Calculate scale factors for covariance transformation
    # This depends on the transform type
    vector_length <- gate_resolution %||%
        gate$xAxis$transform$vectorLength %||% 256

    scale_x <- gate_ellipse_axis_scale(
        gate$xAxis$transform, vector_length, correct_faulty_gate
    )
    scale_y <- gate_ellipse_axis_scale(
        gate$yAxis$transform, vector_length, correct_faulty_gate
    )

    # Scale covariance matrix: Cov_raw = S * Cov_display * S^T
    scale_matrix <- diag(c(scale_x, scale_y))
    cov_raw <- scale_matrix %*% cov_display %*% t(scale_matrix)

    list(cov_raw = cov_raw, center_x_raw = center_x_raw,
        center_y_raw = center_y_raw)
}

#' Build the ellipse covariance matrix in display space
#'
#' From analysis schema 4.0.0 (FlowJo 11.3) on, the semi-axis \code{a} lies
#' along the direction (cos, sin) of the rotation angle and \code{b}
#' perpendicular to it, i.e. Sigma = R diag(a^2, b^2) R^T. This reproduces
#' FlowJo's own counts. Workspaces with an older or unknown schema version keep
#' the convention used by CyFj11 up to 0.99.113 (semi-axes swapped, rotation
#' reversed) so that earlier results stay reproducible.
#'
#' @param a_display Semi-major axis length
#' @param b_display Semi-minor axis length
#' @param angle_rad Rotation angle in radians
#' @param schema_version Analysis schema version (\code{numeric_version} or NULL)
#' @return 2x2 covariance matrix
#' @noRd
gate_ellipse_cov_display <- function(a_display, b_display, angle_rad,
    schema_version = NULL) {
    cos_a <- cos(angle_rad)
    sin_a <- sin(angle_rad)

    if (schema_at_least(schema_version, "4.0.0")) {
        diag_x <- a_display^2 * cos_a^2 + b_display^2 * sin_a^2
        diag_y <- a_display^2 * sin_a^2 + b_display^2 * cos_a^2
    } else {
        diag_x <- a_display^2 * sin_a^2 + b_display^2 * cos_a^2
        diag_y <- a_display^2 * cos_a^2 + b_display^2 * sin_a^2
    }
    off <- (a_display^2 - b_display^2) * sin_a * cos_a

    matrix(c(diag_x, off, off, diag_y), nrow = 2, ncol = 2)
}

#' Per-axis display-to-data scale factor for ellipse covariance
#'
#' For Linear transforms: scale from display to data. For other types
#' (e.g. Biex): scaling is approximately 1 near the center (simplified).
#'
#' @param axis_transform Axis transform spec
#' @param vector_length Display space range
#' @param correct_faulty_gate Fallback maxRange value if maxRange=0
#' @return Numeric scale factor
#' @noRd
gate_ellipse_axis_scale <- function(axis_transform, vector_length,
    correct_faulty_gate) {
    trans_type <- axis_transform$transformType %||% "Linear"

    if (trans_type != "Linear") {
        return(1)
    }

    max_range <- axis_transform$maxRange %||% 262144
    if (max_range == 0 && correct_faulty_gate != 0) {
        max_range <- correct_faulty_gate
    }
    max_range / vector_length
}

#' Convert Range Gate (1D)
#' @keywords internal
#' @importFrom flowCore rectangleGate
convert_range_gate <- function(
    gate, pop_name, extend_val, extend_to,
    correct_faulty_gate = 0, use_transformed_coords = FALSE,
    schema_version = NULL
) {
    param <- gate_range_param(gate, pop_name)
    x_display <- gate_range_vertices(gate, pop_name)

    # Get gate resolution
    gate_resolution <- gate$gateResolution %||% gate$resolution

    # Get transform
    transform_spec <- gate$xAxis$transform %||% gate$transform %||%
        gate$axis$transform

    # Transform to raw space
    x_raw <- display_to_raw(
        x_display, transform_spec, gate_resolution,
        correct_faulty_gate, use_transformed_coords
    )

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
        .gate = matrix(c(min_val, max_val),
            nrow = 2, ncol = 1,
            dimnames = list(c("min", "max"), param)
        )
    )
    # message("convert_range_gate: ", pop_name, "  ", min_val, "  ", max_val)

    return(gate_obj)
}

#' Resolve a range gate's parameter name
#'
#' @param gate Raw gate definition
#' @param pop_name Population name for error messages
#' @return Parameter name (stops when absent or empty)
#' @noRd
gate_range_param <- function(gate, pop_name) {
    # Extract parameter
    param <- gate$parameter %||%
        gate$xParameter %||%
        gate$xAxis$parameterSpec$name %||%
        gate$xAxis %||%
        NULL

    if (is.null(param) || (is.character(param) && nchar(param) == 0)) {
        stop("Range gate missing valid parameter for: ", pop_name)
    }

    param
}

#' Resolve a range gate's display vertices
#'
#' @param gate Raw gate definition
#' @param pop_name Population name for error messages
#' @return Numeric vector of display coordinates (stops when absent)
#' @noRd
gate_range_vertices <- function(gate, pop_name) {
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
    as.numeric(x_display)
}
#' Convert Quadrant Gate
#' @keywords internal
#' @importFrom flowCore quadGate
convert_quadrant_gate <- function(
    gate, pop_name, extend_val, extend_to,
    correct_faulty_gate = 0, use_transformed_coords = FALSE,
    schema_version = NULL
) {
    params <- gate_quadrant_params(gate, pop_name)

    div_display <- gate_quadrant_dividers(gate)

    # Get gate resolution
    gate_resolution <- gate$gateResolution %||% gate$resolution

    # Transform dividers to raw space
    x_div_raw <- display_to_raw(
        div_display$x, gate$xAxis$transform,
        gate_resolution, correct_faulty_gate, use_transformed_coords
    )
    y_div_raw <- display_to_raw(
        div_display$y, gate$yAxis$transform,
        gate_resolution, correct_faulty_gate, use_transformed_coords
    )

    # Apply extension
    x_div_raw <- apply_extension(x_div_raw, extend_val, extend_to)
    y_div_raw <- apply_extension(y_div_raw, extend_val, extend_to)

    # Validate population names
    if (length(pop_name) != 4) {
        stop(
            "Quadrant gate must define exactly 4 populations, got: ",
            length(pop_name)
        )
    }

    # Create boundary
    boundary <- c(x_div_raw, y_div_raw)
    names(boundary) <- c(params$x, params$y)

    # Create quad gate
    base_name <- paste0("quad_", params$x, "_", params$y)

    gate_obj <- flowCore::quadGate(
        filterId = base_name,
        .gate = boundary
    )

    # Store population names
    attr(gate_obj, "pop_names") <- unlist(pop_name) %>% rev()

    return(gate_obj)
}

#' Resolve a quadrant gate's parameter names
#'
#' @param gate Raw gate definition
#' @param pop_name Population name for error messages
#' @return list(x, y) parameter names (stops when either is missing/empty)
#' @noRd
gate_quadrant_params <- function(gate, pop_name) {
    # Extract parameters
    x_param <- gate$xAxis$parameterSpec$name %||% gate$xParameter %||%
        gate$xAxis
    y_param <- gate$yAxis$parameterSpec$name %||% gate$yParameter %||%
        gate$yAxis

    if (is.null(x_param) || is.null(y_param) ||
        (is.character(x_param) && nchar(x_param) == 0) ||
        (is.character(y_param) && nchar(y_param) == 0)) {
        stop(
            "Quadrant gate missing valid parameters for: ",
            paste(pop_name, collapse = ", ")
        )
    }

    list(x = x_param, y = y_param)
}

#' Extract a quadrant gate's divider positions in display space
#'
#' @param gate Raw gate definition
#' @return list(x, y) numeric divider positions
#' @noRd
gate_quadrant_dividers <- function(gate) {
    # Extract divider position
    x_div_display <- gate$xDivider %||% gate$divider$x %||%
        (if (!is.null(gate$xVertices)) {
            unlist(gate$xVertices)[[1]]
        } else NULL) %||%
        gate$x %||% 0

    y_div_display <- gate$yDivider %||% gate$divider$y %||%
        (if (!is.null(gate$yVertices)) {
            unlist(gate$yVertices)[[1]]
        } else NULL) %||%
        gate$y %||% 0

    list(x = as.numeric(x_div_display), y = as.numeric(y_div_display))
}

#' Convert Boolean Gate
#' @keywords internal
#' @importFrom flowWorkspace booleanFilter
#' @importFrom magrittr %>%
convert_boolean_gate <- function(gate, pop_name, schema_version = NULL) {
    specification <- gate_boolean_spec(gate, pop_name)
    if (is.null(specification)) {
        return(NULL)
    }

    # Parse boolean expression
    expr <- tryCatch(
        parse_boolean_expression(specification),
        error = function(e) {
            warning(
                "Failed to parse boolean expression '", specification,
                "' for gate: ", pop_name, ": ", e$message
            )
            return(NULL)
        }
    )

    if (is.null(expr)) {
        return(NULL)
    }

    gate_boolean_filter(expr, pop_name)
}

#' Extract a boolean gate's expression specification
#'
#' @param gate Raw gate definition
#' @param pop_name Population name for warnings
#' @return Specification string, or NULL when absent or empty (warns)
#' @noRd
gate_boolean_spec <- function(gate, pop_name) {
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
        warning(
            "Boolean gate specification had multiple values for: ",
            pop_name, ". Using first value."
        )
    }

    specification
}

#' Build a booleanFilter from a parsed expression
#'
#' @param expr Parsed boolean expression
#' @param pop_name Population name for the filterId and warnings
#' @return flowWorkspace booleanFilter, or NULL on failure (warns)
#' @noRd
gate_boolean_filter <- function(expr, pop_name) {
    tryCatch(
        flowWorkspace::booleanFilter(
            expr = expr,
            filterId = pop_name[[1]]
        ),
        error = function(e) {
            warning(
                "Failed to create boolean filter for: ", pop_name,
                ": ", e$message
            )
            return(NULL)
        }
    )
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
            stop(
                "Failed to parse boolean expression '", expr_string, "': ",
                e$message
            )
        }
    )

    return(expr)
}
