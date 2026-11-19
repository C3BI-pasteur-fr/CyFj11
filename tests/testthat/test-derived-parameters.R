make_derived_zip <- function(n = 6L) {
  d <- tempfile("fjzip"); dir.create(d)
  pdir <- file.path(d, "analyses", "analysis-x", "platform-p1")
  dir.create(pdir, recursive = TRUE)
  # rows deliberately not in event order
  df <- data.frame(`Event #` = n:1, `QC of Experiment Data` = (n:1) %% 2,
                   check.names = FALSE)
  utils::write.csv(df, file.path(pdir, "datasource-S1.csv"), row.names = FALSE)
  utils::write.csv(df, file.path(pdir, "datasource-S1-segmentSummaries.csv"), row.names = FALSE)
  zipfile <- tempfile(fileext = ".flowjo")
  old <- setwd(d); on.exit(setwd(old))
  utils::zip(zipfile, files = "analyses", flags = "-r9Xq")
  zipfile
}

test_that("derived-parameter CSVs are found, ordered by Event # and sized-checked", {
  skip_if(Sys.which("zip") == "")
  zf <- make_derived_zip(6L)
  members <- CyFj11:::list_derived_parameter_files(zf)
  expect_length(members, 1L)                       # segmentSummaries excluded
  m <- CyFj11:::read_derived_parameters(zf, "S1", n_events = 6L)
  expect_equal(colnames(m), "QC of Experiment Data")
  expect_equal(as.vector(m), as.numeric((1:6) %% 2))   # reordered by Event #
  expect_null(CyFj11:::read_derived_parameters(zf, "other"))
  expect_warning(res <- CyFj11:::read_derived_parameters(zf, "S1", n_events = 5L),
                 "skipping")
  expect_null(res)
})

test_that("append_derived_parameters adds the columns to the cytoset", {
  skip_if(Sys.which("zip") == "")
  zf <- make_derived_zip(6L)
  mat <- matrix(as.numeric(1:12), ncol = 2, dimnames = list(NULL, c("FSC-A", "SSC-A")))
  cs <- flowWorkspace::cytoset(list(S1 = flowWorkspace::flowFrame_to_cytoframe(flowCore::flowFrame(mat))))
  out <- CyFj11:::append_derived_parameters(cs, zf, "S1")
  expect_true("QC of Experiment Data" %in% flowCore::colnames(out[[1]]))
  expect_equal(as.vector(flowCore::exprs(out[[1]])[, "QC of Experiment Data"]),
               as.numeric((1:6) %% 2))
  # a workspace without derived parameters leaves the cytoset alone
  expect_identical(CyFj11:::append_derived_parameters(cs, tempfile(), "S1"), cs)
})

test_that("FlowJo 11.3 reference workspace reproduces FlowJo counts", {
  dir <- Sys.getenv("CYFJ11_FJ113_DIR")
  skip_if(dir == "", "CYFJ11_FJ113_DIR not set")
  ws <- read_flowjo11_workspace(file.path(dir, "2026_11_18.flowjo"))
  gsl <- suppressWarnings(fj11_to_gatingset(ws, group_name = 1, path = file.path(dir, "fcs")))
  pop_name <- function(p)
    unlist(ws$populationDefinitions[[p$parents$populationDefinitions[[1]]]]$definition$name)[1]
  ref <- data.frame(
    ds = vapply(ws$populations, function(p) p$parents[["_dataSource"]][[1]], ""),
    pop = vapply(ws$populations, pop_name, ""),
    n = vapply(ws$populations, function(p) p$results$count, numeric(1)))
  fcs_of <- vapply(ws$dataSources, function(d) basename(d$definition$uri), "")
  for (g in gsl) {
    cnt <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(g))
    expect_equal(nrow(cnt), 14L)  # every gate must have been converted
    ds <- names(fcs_of)[fcs_of == flowWorkspace::sampleNames(g)]
    r <- ref[ref$ds == ds, ]
    for (k in seq_len(nrow(cnt))) {
      want <- r$n[r$pop == basename(cnt$Population[k])]
      expect_lte(abs(cnt$Count[k] - want), max(60, 0.01 * want))
    }
  }
})
