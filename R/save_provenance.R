#'Save some information about where data objects came from
#' @importFrom digest digest
save_provenance = function(object){
  name = paste(as.character(substitute(object)), collapse = "")
  dir.create("provenance", showWarnings = FALSE)
  repo_object = git2r::repository(".")
  sink(paste0("provenance/", name, ".txt"))
  print(repo_object)
  print(git2r::status(git2r::repository(".")))
  cat("SHA1: ", digest::digest(object, "sha1"))
  sink(NULL)
}
