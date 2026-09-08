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

#' Pretty Print FlowJo Workspace
#'
#' This function takes a FlowJo v11 workspace file (.flowjo) and extracts the JSON content,
#' then pretty prints it to a text file with "_pretty.json" suffix.
#'
#' @param flowjo_file Path to the FlowJo workspace file (.flowjo)
#' @return Path to the created pretty printed JSON file
#' @export
#'
#' @examples
#' # Pretty print a FlowJo workspace to readable JSON
#' ws_path <- system.file("extdata", "min_test.flowjo", package = "CyFj11")
#' out_file <- pretty_print_flowjo(ws_path)
#' file.exists(out_file)  # Verify output was created
#' # Clean up
#' unlink(out_file)
pretty_print_flowjo <- function(flowjo_file) {
  # Check if the file exists
  if (!file.exists(flowjo_file)) {
    stop("File not found: ", flowjo_file)
  }
  
  # Create temporary directory for extraction
  temp_dir <- tempfile(pattern = "flowjo_extract_")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  
  # Extract the zip file
  system2("unzip", args = c(flowjo_file, "-d", temp_dir))
  
  # Find the analysis JSON file
  json_files <- list.files(file.path(temp_dir, "analyses"), pattern = "\\.json$", recursive = TRUE, full.names = TRUE)
  if (length(json_files) == 0) {
    stop("No JSON files found in the FlowJo workspace")
  }
  
  # Use the first JSON file (assuming there's only one main analysis file)
  json_file <- json_files[1]
  
  json_data <- jsonlite::read_json(json_file, simplifyVector = FALSE)
  
  # Create output file name next to the source file
  base_name <- sub("\\.flowjo$", "", basename(flowjo_file))
  output_file <- paste0(base_name, "_pretty.json")
  
  # Write the pretty printed JSON to file
  jsonlite::write_json(path = output_file, x = json_data, auto_unbox = TRUE, pretty = TRUE, always_decimal = TRUE)

  message("Pretty printed JSON saved to:", output_file, "\n")
  
  return(output_file)
}

# Example usage:
# pretty_print_flowjo("large.mi.flowjo")
