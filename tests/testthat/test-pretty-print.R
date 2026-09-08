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

#' @title Tests for Pretty Print Function
#' @name test-pretty-print
#' @keywords internal
NULL



test_that("pretty_print_flowjo handles missing file", {
  # Test with non-existent file
  expect_error(pretty_print_flowjo("nonexistent.flowjo"),
               "File not found: nonexistent.flowjo")
})

test_that("pretty_print_flowjo works with example data", {
  # Skip if required packages are not available
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("ggcyto")
  
  # Set the path to your FlowJo workspace file
  workspace_path <- system.file("extdata", "test.data.flowjo", package = "CyFj11")
  
  # Skip if file doesn't exist
  skip_if(!file.exists(workspace_path), "Example workspace file not found")
  
  # Read the FlowJo v11 workspace
  ws <- read_flowjo11_workspace(workspace_path)
  
  # Test pretty printing (should not error); run in a temp working dir so
  # the output file does not pollute the package tree
  old_wd <- getwd()
  tmp_out <- tempfile("pretty_print_out_")
  dir.create(tmp_out)
  setwd(tmp_out)
  on.exit(setwd(old_wd), add = TRUE)
  expect_error(pretty_print_flowjo(workspace_path), NA)
  expect_true(file.exists(file.path(tmp_out, "test.data_pretty.json")))
  
  # Examine the workspace structure
  expect_true(is.list(ws))
  expect_true(all(c("groups", "dataSources", "populationDefinitions", "populations") %in% names(ws)))
})
