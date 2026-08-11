suppressWarnings(suppressMessages(pkgload::load_all(".", quiet=TRUE, export_all=FALSE)))
options(warn=-1)
# Load every old definition into one env, mimicking the harness's own order so
# the winner of any collision is the same one the golden captured.
e <- new.env()
sf <- list.files("scripts","[.]R$",recursive=TRUE,full.names=TRUE)
ex <- grep("build_reference_databases[.]R$", sf, value=TRUE)
for (f in setdiff(sf,ex)) try(suppressWarnings(suppressMessages(sys.source(f, envir=e))), silent=TRUE)
for (f in ex) for (q in parse(f)) {
  ok <- is.call(q) && length(q)==3L && as.character(q[[1]]) %in% c("<-","=") &&
        is.call(q[[3]]) && identical(as.character(q[[3]][[1]]),"function")
  if (ok) try(eval(q, e), silent=TRUE)
}
frozen <- c("normalize_gsea_results","run_gsea","create_standard_volcano","format_pathway_name",
"gsea_running_sum_plot","list_to_term2gene","gsea_barplot","gsea_dotplot","load_reference_db",
"custom_minimal_theme_with_grid","gsea_dotplot_facet","create_MD_plot","empty_gsea_tibble",
"ensure_dir","run_gsea_analysis","save_gsea_log","plot_all_gsea_results","convert_human_to_mouse",
"parse_gmx","parse_mitoxplorer","filter_by_size","build_dge","list_reference_dbs","download_gatom_references")
fmt <- function(f) {
  fo <- formals(f)
  paste(mapply(function(n,v) {
    d <- paste(deparse(v), collapse="")
    if (identical(d,"")) n else paste0(n,"=",d)
  }, names(fo), fo), collapse=", ")
}
bad <- 0
for (nm in frozen) {
  o <- get0(nm, envir=e, mode="function"); n <- get0(nm, envir=asNamespace("bulkiRNA"), mode="function")
  if (is.null(o)) { cat(sprintf("%-32s OLD MISSING\n", nm)); next }
  if (is.null(n)) { cat(sprintf("%-32s NEW MISSING  <-- \n", nm)); bad<-bad+1; next }
  fo <- fmt(o); fn <- fmt(n)
  if (identical(names(formals(o)), names(formals(n))) && identical(fo,fn)) {
    cat(sprintf("%-32s identical\n", nm))
  } else if (identical(names(formals(o)), names(formals(n)))) {
    cat(sprintf("%-32s DEFAULTS DIFFER <--\n   old: %s\n   new: %s\n", nm, fo, fn)); bad<-bad+1
  } else {
    cat(sprintf("%-32s FORMALS DIFFER <--\n   old: %s\n   new: %s\n", nm, fo, fn)); bad<-bad+1
  }
}
cat("\nnames with a difference:", bad, "of", length(frozen), "\n")
