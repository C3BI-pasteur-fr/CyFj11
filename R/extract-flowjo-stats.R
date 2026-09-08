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

#' @title Extract Statistics from FlowJo Workspace (.wsp) XML Files
#' @name extract-flowjo-stats
#' @keywords internal
#'
#' @importFrom xml2 read_xml xml_ns_strip xml_find_all xml_attr xml_attrs xml_name xml_type xml_parent
#' @importFrom dplyr anti_join arrange bind_rows case_when coalesce count desc distinct enquo filter group_by if_else left_join mutate n select summarise ungroup
#' @importFrom readr read_csv write_csv
#' @importFrom tidyr pivot_wider
#' @importFrom rlang %||%
#' @importFrom tibble tibble
#' @importFrom utils capture.output
#' @importFrom magrittr %>%
NULL

utils::globalVariables(c(
  "better_name", "description", "file", "id", "missing_signatures",
  "population", "population_path", "sample_count", "sample_id",
  "sample_name", "stat_ancestor", "stat_name", "value_num", "value_raw"
))

#' Extract Statistics from FlowJo Workspace (.wsp) XML Files
#'
#' Parses FlowJo workspace XML files and extracts population counts (from the
#' \code{count} attribute of \code{<Population>}) and all \code{<Statistic>}
#' elements. The function implements a two-step workflow:
#'
#' \enumerate{
#'   \item Run \strong{without} \code{csv_file} to generate a template CSV
#'     containing every unique statistic signature found in the workspaces.
#'     The \code{better_name} column is pre-filled automatically; edit only
#'     rows that collide or need a more human-readable label.
#'   \item Run \strong{with} \code{csv_file} to validate the mapping and return
#'     a wide-format table: one row per sample, one column per renamed statistic.
#' }
#'
#' @param wsp_files Character vector of paths to \code{.wsp} files.
#' @param csv_file Optional path to a mapping CSV. Must contain the columns
#'   \code{population_path}, \code{stat_name}, and \code{stat_ancestor}.
#'   The \code{better_name} column is optional; if missing it is auto-generated
#'   from \code{population}, a short form of \code{stat_name}, and
#'   \code{stat_ancestor}.
#' @param output_csv Path where the template CSV is written when
#'   \code{csv_file} is \code{NULL}. Default: \code{"wsp_stat_mapping.csv"}.
#' @param write_csv Logical. If \code{TRUE} (default) and \code{csv_file} is
#'   \code{NULL}, write the template CSV to disk.
#' @param value_as_numeric Logical. If \code{TRUE} (default), coerce raw values
#'   to numeric where possible.
#'
#' @return
#' If \code{csv_file} is \code{NULL}, a named list:
#' \describe{
#'   \item{\code{long}}{A tibble with all extracted counts/statistics in long
#'     format.}
#'   \item{\code{mapping}}{A tibble of unique statistic signatures, ready to be
#'     written to \code{output_csv}.}
#' }
#'
#' If \code{csv_file} is provided, a wide tibble with one row per sample and
#' one column per renamed statistic.
#'
#' @details
#' The matching signature is the combination of \code{population_path},
#' \code{stat_name}, and \code{stat_ancestor}. \code{population_path} is the
#' full gating lineage (e.g., \code{"CD45+ / T cells / CD4+"}). The \code{id}
#' column is included for reference but is \emph{not} used for matching, so it
#' can be left blank in the CSV if desired.
#'
#' @examples
#' \dontrun{
#' # Step 1: discover all statistics and write a template CSV
#' result <- extract_flowjo_stats(wsp_files = c("panel1.wsp", "panel2.wsp"))
#'
#' # Step 2: edit wsp_stat_mapping.csv, then:
#' wide <- extract_flowjo_stats(
#'   wsp_files = c("panel1.wsp", "panel2.wsp"),
#'   csv_file = "wsp_stat_mapping.csv"
#' )
#' }
#'
#' @param preserve_slashes Logical. If TRUE, "/" characters in population names
#'   are preserved instead of being replaced with "_". Default FALSE.
#' @keywords internal
extract_flowjo_stats <- function(wsp_files,
                                 csv_file = NULL,
                                 output_csv = "wsp_stat_mapping.csv",
                                 write_csv = TRUE,
                                 value_as_numeric = TRUE,
                                 preserve_slashes = FALSE) {

  # ---- Input validation ----
  if (!is.character(wsp_files) || length(wsp_files) == 0) {
    stop("'wsp_files' must be a non-empty character vector.")
  }

  missing_files <- wsp_files[!file.exists(wsp_files)]
  if (length(missing_files) > 0) {
    stop("The following files do not exist:\n",
         paste(missing_files, collapse = "\n"))
  }

  # ---- Extract long-format table ----
  long_table <- bind_rows(lapply(wsp_files, extract_wsp_data_long))

  if (value_as_numeric) {
    long_table <- long_table %>%
      mutate(value_num = as_num_quiet(value_raw))
  }

  # Normalize empty ancestors to NA for matching
  long_table <- long_table %>%
    mutate(stat_ancestor = if_else(is.na(stat_ancestor) | stat_ancestor == "",
                                   NA_character_,
                                   as.character(stat_ancestor)))

  # ---- Mode 1: no CSV -> generate mapping template ----
  if (is.null(csv_file)) {

    mapping <- long_table %>%
      distinct(population_path, population, id, stat_name, stat_ancestor) %>%
      mutate(
        better_name = make_better_names(population, stat_name, stat_ancestor,
                                        preserve_slashes = preserve_slashes),
        description = case_when(
          stat_name == "PopulationCount" ~
            paste0("Event count for population '", population_path, "'"),
          is.na(stat_ancestor) ~
            paste0(stat_name, " for population '", population_path, "'"),
          TRUE ~
            paste0(stat_name, " for population '", population_path,
                   "' relative to '", stat_ancestor, "'")
        )
      ) %>%
      select(population_path, population, id, stat_name, stat_ancestor,
             better_name, description) %>%
      arrange(population_path, stat_name, stat_ancestor)

    if (write_csv) {
      write_csv(mapping, output_csv)
      message("Mapping template written to: ",
              normalizePath(output_csv, mustWork = FALSE))
    }

    return(list(long = long_table, mapping = mapping))
  }

  # ---- Mode 2: CSV provided -> validate and build wide table ----
  if (!file.exists(csv_file)) {
    stop("csv_file does not exist: ", csv_file)
  }

  mapping <- read_csv(csv_file, show_col_types = FALSE)

  required_cols <- c("population_path", "stat_name", "stat_ancestor")
  missing_cols <- required_cols[!required_cols %in% names(mapping)]
  if (length(missing_cols) > 0) {
    stop("csv_file is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Auto-generate better_name if the column is missing from the CSV
  if (!"better_name" %in% names(mapping)) {
    mapping <- mapping %>%
      mutate(
        population = coalesce(population, get_leaf_population(population_path)),
        better_name = make_better_names(population, stat_name, stat_ancestor,
                                        preserve_slashes = preserve_slashes)
      )
  }

  # Normalize mapping signatures
  mapping <- mapping %>%
    mutate(
      stat_ancestor = if_else(is.na(stat_ancestor) | stat_ancestor == "",
                              NA_character_,
                              as.character(stat_ancestor)),
      better_name = as.character(better_name)
    )

  # Signatures present in data vs. CSV
  data_sigs <- long_table %>%
    distinct(population_path, stat_name, stat_ancestor)

  csv_sigs <- mapping %>%
    distinct(population_path, stat_name, stat_ancestor)

  # Fail if the workspace contains stats not in the CSV
  missing_in_csv <- long_table %>%
    anti_join(csv_sigs, by = c("population_path", "stat_name", "stat_ancestor")) %>%
    distinct(file, population_path, stat_name, stat_ancestor)

  if (nrow(missing_in_csv) > 0) {
    affected_files <- sort(unique(missing_in_csv$file))
    summary_by_file <- missing_in_csv %>%
      count(file, name = "n_missing") %>%
      arrange(file)

    stop(
      "The following workspace signatures are missing from csv_file.\n",
      "This usually means the CSV was generated without all wsp_files.\n",
      "Affected file(s): ", paste(affected_files, collapse = ", "), "\n\n",
      "Missing by file:\n",
      paste(capture.output(format(as.data.frame(summary_by_file))),
            collapse = "\n"),
      "\n\nMissing signatures:\n",
      paste(capture.output(format(as.data.frame(missing_in_csv))),
            collapse = "\n"),
      call. = FALSE
    )
  }

  # Warn about extra mappings in CSV that are not in the data
  extra_in_csv <- csv_sigs %>%
    anti_join(data_sigs, by = c("population_path", "stat_name", "stat_ancestor"))

  if (nrow(extra_in_csv) > 0) {
    warning("csv_file contains mappings not present in the workspaces ",
            "(these will be ignored):\n",
            paste(capture.output(format(as.data.frame(extra_in_csv))),
                  collapse = "\n"))
  }

  # Check for missing better_name values
  missing_names <- mapping %>%
    filter(is.na(better_name) | better_name == "")

  if (nrow(missing_names) > 0) {
    stop("The following rows in csv_file have empty better_name values:\n",
         paste(capture.output(format(as.data.frame(missing_names))),
               collapse = "\n"))
  }

  # Check for duplicate better_name values and require manual resolution
  dup_count <- mapping %>%
    count(better_name) %>%
    filter(n > 1)

  if (nrow(dup_count) > 0) {
    stop("Duplicate better_name values found; please resolve manually:\n",
         paste(dup_count$better_name, collapse = "\n"))
  }

  # Join long table with mapping and pivot wide
  wide_table <- long_table %>%
    left_join(
      mapping %>% select(population_path, stat_name, stat_ancestor, better_name),
      by = c("population_path", "stat_name", "stat_ancestor")
    ) %>%
    mutate(
      better_name = if_else(is.na(better_name),
                            paste0(population_path, " | ", stat_name),
                            better_name)
    ) %>%
    select(file, sample_name, sample_id, sample_count, better_name, value_num) %>%
    pivot_wider(
      id_cols = c(file, sample_name, sample_id, sample_count),
      names_from = better_name,
      values_from = value_num
    )

  return(wide_table)
}


#' Check Workspace-to-CSV Coverage
#'
#' Reports which \code{.wsp} files have statistics that are missing from a
#' mapping CSV. Useful for detecting mismatches before calling
#' \code{extract_flowjo_stats()} in wide-table mode.
#'
#' @param wsp_files Character vector of paths to \code{.wsp} files.
#' @param csv_file Path to a mapping CSV with the columns
#'   \code{population_path}, \code{stat_name}, and \code{stat_ancestor}.
#'
#' @return A tibble with one row per \code{wsp_file} and columns:
#' \describe{
#'   \item{file}{The workspace file name.}
#'   \item{total_signatures}{Number of unique statistic signatures in the file.}
#'   \item{missing_signatures}{Number of signatures not present in the CSV.}
#'   \item{covered}{Logical; \code{TRUE} if all signatures are in the CSV.}
#' }
#'
#' @keywords internal
check_wsp_csv_coverage <- function(wsp_files, csv_file) {
  if (!is.character(wsp_files) || length(wsp_files) == 0) {
    stop("'wsp_files' must be a non-empty character vector.")
  }
  if (!is.character(csv_file) || length(csv_file) != 1) {
    stop("'csv_file' must be a single character string.")
  }
  if (!file.exists(csv_file)) {
    stop("csv_file does not exist: ", csv_file)
  }

  csv <- read_csv(csv_file, show_col_types = FALSE) %>%
    mutate(
      stat_ancestor = if_else(is.na(stat_ancestor) | stat_ancestor == "",
                              NA_character_,
                              as.character(stat_ancestor))
    ) %>%
    distinct(population_path, stat_name, stat_ancestor)

  long_table <- bind_rows(lapply(wsp_files, extract_wsp_data_long)) %>%
    mutate(
      stat_ancestor = if_else(is.na(stat_ancestor) | stat_ancestor == "",
                              NA_character_,
                              as.character(stat_ancestor))
    )

  file_totals <- long_table %>%
    distinct(file, population_path, stat_name, stat_ancestor) %>%
    group_by(file) %>%
    summarise(total_signatures = n(), .groups = "drop")

  missing_by_file <- long_table %>%
    distinct(file, population_path, stat_name, stat_ancestor) %>%
    anti_join(csv, by = c("population_path", "stat_name", "stat_ancestor")) %>%
    group_by(file) %>%
    summarise(missing_signatures = n(), .groups = "drop")

  coverage <- file_totals %>%
    left_join(missing_by_file, by = "file") %>%
    mutate(
      missing_signatures = coalesce(missing_signatures, 0L),
      covered = missing_signatures == 0L
    ) %>%
    arrange(desc(missing_signatures), file)

  coverage
}


# ============================================================================
# Helper functions
# ============================================================================

#' Build the full population path from a node up to (but not above) SampleNode
#' @noRd
get_population_path <- function(node) {
  path <- character()
  current <- node

  while (!is.null(current) && xml_type(current) == "element") {
    nm <- xml_name(current)
    if (nm == "SampleNode") break  # stop at the sample boundary
    if (nm == "Population") {
      path <- c(xml_attr(current, "name"), path)
    }
    current <- xml_parent(current)
  }

  if (length(path) == 0) return(NA_character_)
  paste(path, collapse = " / ")
}

#' Get the name of the immediate enclosing Population
#' @noRd
get_immediate_population <- function(node) {
  if (xml_type(node) == "element" && xml_name(node) == "Population") {
    return(xml_attr(node, "name"))
  }

  current <- xml_parent(node)
  while (!is.null(current) && xml_type(current) == "element") {
    nm <- xml_name(current)
    if (nm == "SampleNode") return(NA_character_)
    if (nm == "Population") return(xml_attr(current, "name"))
    current <- xml_parent(current)
  }
  return(NA_character_)
}

#' Extract all counts and statistics from one .wsp file into a long tibble
#' @noRd
extract_wsp_data_long <- function(file) {
  doc <- xml_ns_strip(read_xml(file))
  samples <- xml_find_all(doc, "//SampleNode")

  bind_rows(lapply(samples, function(sample) {
    sample_name  <- xml_attr(sample, "name")
    sample_id    <- xml_attr(sample, "sampleID")
    sample_count <- as_num_quiet(xml_attr(sample, "count"))

    # ---- Population-level counts (the count attribute) ----
    populations <- xml_find_all(sample, ".//Population")

    pop_counts <- bind_rows(lapply(populations, function(pop) {
      attrs <- xml_attrs(pop)
      tibble(
        file = basename(file),
        sample_name = sample_name,
        sample_id = sample_id,
        sample_count = sample_count,
        population_path = get_population_path(pop),
        population = attrs["name"] %||% NA_character_,
        id = attrs["id"] %||% NA_character_,
        field_type = "PopulationCount",
        stat_name = "PopulationCount",
        stat_ancestor = NA_character_,
        value_raw = attrs["count"]
      )
    }))

    # ---- Statistic elements ----
    stats <- xml_find_all(sample, ".//Statistic")

    stat_rows <- bind_rows(lapply(stats, function(stat) {
      attrs <- xml_attrs(stat)
      tibble(
        file = basename(file),
        sample_name = sample_name,
        sample_id = sample_id,
        sample_count = sample_count,
        population_path = get_population_path(stat),
        population = get_immediate_population(stat),
        id = attrs["id"] %||% NA_character_,
        field_type = "Statistic",
        stat_name = attrs["name"] %||% NA_character_,
        stat_ancestor = attrs["ancestor"] %||% NA_character_,
        value_raw = attrs["value"]
      )
    }))

    bind_rows(pop_counts, stat_rows)
  }))
}

#' Extract the leaf population name from a population path
#' @noRd
get_leaf_population <- function(population_path) {
  path <- as.character(population_path)
  leaf <- sub("^.*/\\s*", "", path)
  ifelse(leaf == "", NA_character_, leaf)
}

#' Shorten common FlowJo statistic names for compact column names
#' @noRd
shorten_stat_name <- function(stat_name) {
  nm <- as.character(stat_name)
  nm[is.na(nm)] <- ""

  # Named substitutions (long -> short)
  dict <- c(
    "PopulationCount" = "Count",
    "Frequency of Parent" = "FreqParent",
    "Freq. of Parent" = "FreqParent",
    "Frequency of Grandparent" = "FreqGrandparent",
    "Freq. of Grandparent" = "FreqGrandparent",
    "Frequency of Total" = "FreqTotal",
    "Freq. of Total" = "FreqTotal",
    "Mean" = "Mean",
    "Median" = "Median",
    "Geometric Mean" = "GeoMean",
    "Mode" = "Mode",
    "SD" = "SD",
    "CV" = "CV",
    "Min" = "Min",
    "Max" = "Max",
    "Count" = "Count"
  )

  out <- dict[nm]
  out <- ifelse(is.na(out), nm, out)

  # Remove punctuation / spacing to make it compact and name-safe
  out <- gsub("[[:space:]_]+", "", out)
  out <- gsub("[^A-Za-z0-9]", "", out)
  unname(out)
}

#' Sanitize a string so it is safe as a column name
#'
#' Keeps alphanumeric characters, plus/minus signs, and underscores. Replaces
#' everything else with underscores, collapses repeated separators, and makes
#' sure the name does not start with a digit.
#'
#' @param x Character vector to sanitize
#' @param preserve_slashes Logical. If TRUE, "/" characters are preserved
#'   instead of being replaced with "_". Default FALSE (replaces "/" with "_").
#' @return Character vector of sanitized names
#' @noRd
sanitize_name <- function(x, preserve_slashes = FALSE) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  # Handle slashes: either preserve or replace
  if (!preserve_slashes) {
    x <- gsub("/", "_", x)
  }

  # Keep only allowed characters (preserve +/- distinction and optionally /)
  if (preserve_slashes) {
    x <- gsub("[^A-Za-z0-9+/_-]", "_", x)
  } else {
    x <- gsub("[^A-Za-z0-9+_-]", "_", x)
  }

  # Collapse repeated separators
  x <- gsub("_+", "_", x)
  x <- gsub("(^_|_$)", "", x)
  # Ensure it does not start with a digit
  x <- ifelse(grepl("^[0-9]", x), paste0("X", x), x)
  # Empty strings become NA
  unname(ifelse(x == "", NA_character_, x))
}

#' Build better_name values from population, stat_name, and stat_ancestor
#'
#' The names are guaranteed to be unique. If the sanitized combination is not
#' unique, a minimal numeric suffix is appended.
#'
#' @param population Population name
#' @param stat_name Statistic name
#' @param stat_ancestor Ancestor population name (if any)
#' @param preserve_slashes Logical. If TRUE, "/" characters are preserved.
#'   Default FALSE (replaces "/" with "_").
#' @return Character vector of sanitized names
#' @noRd
make_better_names <- function(population, stat_name, stat_ancestor, preserve_slashes = FALSE) {
  pop <- sanitize_name(population, preserve_slashes = preserve_slashes)
  short <- shorten_stat_name(stat_name)
  anc <- sanitize_name(stat_ancestor, preserve_slashes = preserve_slashes)

  has_ancestor <- !is.na(anc) & anc != ""
  anc_part <- ifelse(has_ancestor, paste0("_rel_", anc), "")

  candidates <- paste0(pop, "_", short, anc_part)
  candidates <- sanitize_name(candidates, preserve_slashes = preserve_slashes)

  # Resolve collisions by appending a numeric suffix
  seen <- character()
  result <- character(length(candidates))

  for (i in seq_along(candidates)) {
    base <- candidates[i]
    if (is.na(base)) {
      result[i] <- NA_character_
      next
    }
    candidate <- base
    suffix <- 1L
    while (candidate %in% seen) {
      suffix <- suffix + 1L
      candidate <- paste0(base, "_", suffix)
    }
    result[i] <- candidate
    seen <- c(seen, candidate)
  }

  unname(result)
}
