# Read locations and cohort membership without putting private information in Git.
# This helper performs no expression analysis or deconvolution.
pdac_config <- function() {
  filename <- Sys.getenv("PDAC_PRIVATE_CONFIG")
  if (!nzchar(filename)) stop("Set PDAC_PRIVATE_CONFIG; see docs/inputs.md")
  code_root <- Sys.getenv("PDAC_CODE_ROOT")
  if (!nzchar(code_root)) stop("Set PDAC_CODE_ROOT to this repository")
  repository <- normalizePath(code_root, mustWork = TRUE)
  config_path <- normalizePath(filename, mustWork = TRUE)
  if (identical(config_path, repository) || startsWith(config_path, paste0(repository, "/"))) {
    stop("The populated private configuration must be outside the repository")
  }
  jsonlite::fromJSON(filename, simplifyVector = TRUE)
}

pdac_path <- function(name) {
  value <- pdac_config()$paths[[name]]
  if (is.null(value) || !nzchar(value) || startsWith(value, "/path/to/")) {
    stop("Supply a real private location for paths.", name)
  }
  code_root <- Sys.getenv("PDAC_CODE_ROOT")
  if (!nzchar(code_root)) stop("Set PDAC_CODE_ROOT to this repository")
  result <- normalizePath(path.expand(value), mustWork = FALSE)
  repository <- normalizePath(code_root, mustWork = TRUE)
  if (identical(result, repository) || startsWith(result, paste0(repository, "/"))) {
    stop("Private inputs and outputs must stay outside the repository")
  }
  result
}

pdac_cohort <- function(name) {
  values <- as.character(pdac_config()$cohorts[[name]])
  if (!length(values) || anyDuplicated(values) || any(!grepl("^[0-9]+$", values))) {
    stop("Populate a unique ordered numeric patient list in cohorts.", name)
  }
  values
}
