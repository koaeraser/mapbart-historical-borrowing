# Run a repo script verbatim with run_all.R-style header overrides.
REPO <- "/Users/yuanj/Dropbox (Personal)/YuanJi/research/Yunxuan-Zhang/MAP-BART/yunxuan-repo"
SCR  <- "/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon"
format_r_value <- function(x) { if (is.null(x)) return("NULL"); if (is.character(x) && length(x) == 1) return(paste0('"', x, '"')); paste(deparse(x), collapse = " ") }
apply_overrides <- function(text, overrides) {
  for (sk in c("seed_rct", "seed_method")) if (!is.null(overrides[[sk]])) {
    text <- gsub("set\\.seed\\(\\s*\\d+\\s*\\)", paste0("set.seed(", overrides[[sk]], ")"), text, perl = TRUE)
    overrides[[sk]] <- NULL; break }
  for (name in names(overrides)) { val <- overrides[[name]]; if (is.null(val)) next
    pattern <- paste0("(?m)^(\\s*)\\b", name, "\\b\\s*<-\\s*[^\n]*$")
    text <- gsub(pattern, paste0("\\1", name, " <- ", format_r_value(val)), text, perl = TRUE) }
  text
}
run_repo_script <- function(folder, file, overrides, stop_at = NULL) {
  txt <- paste(readLines(file.path(REPO, folder, file), warn = FALSE), collapse = "\n")
  txt <- sub("library(Rlab)", "rbern <- function(n, prob) rbinom(n, 1, prob)", txt, fixed = TRUE)
  txt <- sub("rm(list = ls())", "", txt, fixed = TRUE); txt <- sub("rm(list=ls())", "", txt, fixed = TRUE)
  if (!is.null(stop_at)) txt <- sub(stop_at, paste0("stop('__STOP__')\n", stop_at), txt, fixed = TRUE)
  txt <- apply_overrides(txt, c(list(mainDir = paste0(SCR, "/")), overrides))
  e <- new.env(parent = globalenv())
  tryCatch(eval(parse(text = txt), envir = e), error = function(err) if (!grepl("__STOP__", conditionMessage(err))) stop(err))
  invisible(e)
}
