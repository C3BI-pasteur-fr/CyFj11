test_that("fj11_schema_version parses workspaces, strings and NULL", {
  expect_null(fj11_schema_version(NULL))
  expect_null(fj11_schema_version(""))
  expect_null(fj11_schema_version("not a version"))
  expect_s3_class(fj11_schema_version("4.0.0"), "numeric_version")
  expect_true(fj11_schema_version("4.0.0") >= "3.0.0")
  expect_equal(fj11_schema_version(list(schemaVersion = "3.0.0")),
               numeric_version("3.0.0"))
})

test_that("schema_at_least treats unknown versions as legacy", {
  expect_false(schema_at_least(NULL, "3.0.0"))
  expect_true(schema_at_least("4.0.0", "4.0.0"))
  expect_false(schema_at_least("3.0.0", "4.0.0"))
})

test_that("check_schema_version warns for untested versions only", {
  expect_silent(check_schema_version(numeric_version("3.0.0")))
  expect_silent(check_schema_version(numeric_version("4.0.0")))
  expect_warning(check_schema_version(numeric_version("9.0.0")), "newer")
  expect_warning(check_schema_version(numeric_version("1.0.0")), "older")
  expect_warning(check_schema_version(NULL), "schemaVersion")
})

test_that("read_flowjo11_workspace stores the schema version", {
  ws <- read_flowjo11_workspace(system.file("extdata", "min_test.flowjo", package = "CyFj11"))
  expect_s3_class(ws$schema_version, "numeric_version")
  expect_equal(fj11_schema_version(ws), ws$schema_version)
})
