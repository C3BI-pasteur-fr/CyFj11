#' @title Derived parameters stored inside a FlowJo v11 workspace
#' @description
#' From FlowJo 11.3 on, gates can be drawn on parameters that FlowJo algorithms
#' derive (event quality filter, FlowSOM, t-SNE, UMAP). Their per-event values are
#' not in the FCS files; they are stored as CSV files inside the \code{.flowjo}
#' archive, one per algorithm ("platform") and sample:
#' \code{analyses/analysis-<id>/platform-<uuid>/datasource-<sample uuid>.csv}.
#' The first column is \code{Event #} (1-based index into the FCS file), the
#' remaining columns are named like the gate parameters.
#' @name derived-parameters
#' @keywords internal
NULL

#' List derived-parameter CSV members of a workspace archive
#'
#' @param workspace_path Path to the \code{.flowjo} archive.
#' @return Character vector of archive member names (excluding segment summaries).
#' @keywords internal
#' @noRd
list_derived_parameter_files <- function(workspace_path) {
  if (is.null(workspace_path) || !file.exists(workspace_path)) return(character())
  members <- tryCatch(utils::unzip(workspace_path, list = TRUE)$Name,
                      error = function(e) character())
  files <- grep("(^|/)platform-[^/]+/datasource-[^/]+\\.csv$", members, value = TRUE)
  files[!grepl("-segmentSummaries\\.csv$", files)]
}

#' Read the derived parameters of one sample
#'
#' @param workspace_path Path to the \code{.flowjo} archive.
#' @param sample_uuid Data source (sample) UUID.
#' @param n_events Expected number of events (rows). Files with another row
#'   count are skipped with a warning.
#' @param members Optional result of \code{list_derived_parameter_files()}.
#' @return A numeric matrix with one row per event (ordered by \code{Event #})
#'   and one column per derived parameter, or \code{NULL} if there are none.
#' @keywords internal
#' @noRd
read_derived_parameters <- function(workspace_path, sample_uuid, n_events = NULL,
                                    members = list_derived_parameter_files(workspace_path)) {
  files <- members[grepl(paste0("/datasource-", sample_uuid, "\\.csv$"), members)]
  if (length(files) == 0L) return(NULL)

  cols <- lapply(files, function(f) {
    # read.csv opens and closes an unopened connection itself
    x <- tryCatch(utils::read.csv(unz(workspace_path, f), check.names = FALSE),
                  error = function(e) NULL)
    if (is.null(x) || !"Event #" %in% names(x)) return(NULL)
    if (!is.null(n_events) && nrow(x) != n_events) {
      warning("Derived parameters in '", f, "' have ", nrow(x), " rows but the FCS file has ",
              n_events, " events; skipping them.", call. = FALSE)
      return(NULL)
    }
    x <- x[order(x[["Event #"]]), setdiff(names(x), "Event #"), drop = FALSE]
    as.matrix(x)
  })
  cols <- Filter(Negate(is.null), cols)
  if (length(cols) == 0L) return(NULL)
  m <- do.call(cbind, cols)
  storage.mode(m) <- "double"   # cf_append_cols needs doubles (CSV may hold integers)
  # keep the first occurrence of a duplicated column name
  m[, !duplicated(colnames(m)), drop = FALSE]
}

#' Append derived parameters to the cytoframes of a cytoset
#'
#' Samples are matched by position to \code{sample_uuids}. Columns that already
#' exist in a cytoframe are left untouched.
#'
#' @param cytoset A \code{cytoset}.
#' @param workspace_path Path to the \code{.flowjo} archive.
#' @param sample_uuids Sample UUIDs, in the order of the cytoset.
#' @return The cytoset with derived parameters appended.
#' @keywords internal
#' @noRd
append_derived_parameters <- function(cytoset, workspace_path, sample_uuids) {
  members <- list_derived_parameter_files(workspace_path)
  if (length(members) == 0L) return(cytoset)

  sn <- flowWorkspace::sampleNames(cytoset)
  cfs <- lapply(seq_along(sn), function(i) cytoset[[i]])
  names(cfs) <- sn
  changed <- FALSE
  for (i in seq_len(min(length(sample_uuids), length(cfs)))) {
    cf <- cfs[[i]]
    m <- read_derived_parameters(workspace_path, sample_uuids[[i]],
                                 n_events = nrow(cf), members = members)
    if (is.null(m)) next
    keep <- !colnames(m) %in% flowCore::colnames(cf)
    if (!any(keep)) next
    if (.pkgenv$verbose) message("Appending derived parameters: ", # nocov
                                 paste(colnames(m)[keep], collapse = ", ")) # nocov
    cfs[[i]] <- flowWorkspace::cf_append_cols(cf, m[, keep, drop = FALSE])
    changed <- TRUE
  }
  if (!changed) return(cytoset)
  # A cytoset needs identical columns in all samples, so it is rebuilt
  tryCatch(flowWorkspace::cytoset(cfs),
           error = function(e) stop("Could not append derived parameters: ", conditionMessage(e),
                                    "\nAll samples must have the same derived parameters.",
                                    call. = FALSE))
}
