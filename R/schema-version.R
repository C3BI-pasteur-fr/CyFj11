#' @title FlowJo v11 analysis schema version helpers
#' @description
#' FlowJo v11 workspaces carry a \code{schemaVersion} field in the analysis JSON
#' (for example \code{"3.0.0"} or \code{"4.0.0"}). The version is stored on the
#' workspace object and handed to the conversion functions so that each of them
#' can decide whether version specific behaviour is needed.
#' @name schema-version
#' @keywords internal
NULL

# Schema versions the package has been tested against.
.tested_schema_versions <- c("3.0.0", "4.0.0")

#' Get the schema version of a FlowJo v11 workspace
#'
#' @param x A \code{flowjo11_workspace} (from \code{\link{read_flowjo11_workspace}}),
#'   a character string such as \code{"4.0.0"}, a \code{numeric_version}, or
#'   \code{NULL}.
#' @return A \code{numeric_version}, or \code{NULL} if the version is unknown.
#' @examples
#' ws <- read_flowjo11_workspace(system.file("extdata", "min_test.flowjo", package = "CyFj11"))
#' fj11_schema_version(ws)
#' fj11_schema_version("4.0.0") >= "4.0.0"
#' @export
fj11_schema_version <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "numeric_version")) return(x)
  if (is.list(x)) x <- x$schemaVersion
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))) return(NULL)
  tryCatch(numeric_version(as.character(x)), error = function(e) NULL)
}

#' Is the schema version at least a given version?
#'
#' Unknown versions (\code{NULL}) are treated as older than any version, so code
#' guarded by \code{schema_at_least()} keeps the legacy behaviour for them.
#'
#' @param schema_version A \code{numeric_version}, character, workspace or \code{NULL}.
#' @param min_version Character. Minimum version, e.g. \code{"4.0.0"}.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
schema_at_least <- function(schema_version, min_version) {
  v <- fj11_schema_version(schema_version)
  if (is.null(v)) return(FALSE)
  v >= numeric_version(min_version)
}

#' Warn if a workspace uses an untested schema version
#'
#' @param schema_version A \code{numeric_version} or \code{NULL}.
#' @return Invisibly, \code{TRUE} if the version is known and not newer than the
#'   newest tested version.
#' @keywords internal
#' @noRd
check_schema_version <- function(schema_version) {
  v <- fj11_schema_version(schema_version)
  if (is.null(v)) {
    warning("Workspace has no recognisable schemaVersion; assuming the oldest ",
            "supported format (", .tested_schema_versions[1], ").", call. = FALSE)
    return(invisible(FALSE))
  }
  tested <- numeric_version(.tested_schema_versions)
  if (v > max(tested)) {
    warning("Workspace schemaVersion ", as.character(v), " is newer than the newest ",
            "version CyFj11 was tested with (", as.character(max(tested)), "). ",
            "Results may be incomplete or wrong.", call. = FALSE)
    return(invisible(FALSE))
  }
  if (v < min(tested)) {
    warning("Workspace schemaVersion ", as.character(v), " is older than the oldest ",
            "version CyFj11 was tested with (", as.character(min(tested)), ").",
            call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
