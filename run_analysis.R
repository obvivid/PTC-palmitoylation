#!/usr/bin/env Rscript
# Portable entry point. No changes are made outside the supplied project root.
if (.Platform$OS.type == "windows") suppressWarnings(try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE))
args <- commandArgs(trailingOnly=TRUE)
all_args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep("^--file=",all_args,value=TRUE)
if (!length(file_arg)) stop("Run this file with Rscript.")
SOURCE_ROOT <- dirname(normalizePath(sub("^--file=","",file_arg[1]),mustWork=TRUE))
PROJECT_ROOT <- SOURCE_ROOT
selected <- NULL
dry <- FALSE
for (a in args) {
  if (a == "--dry-run") dry <- TRUE
  else if (startsWith(a,"--root=")) PROJECT_ROOT <- normalizePath(sub("^--root=","",a),mustWork=TRUE)
  else if (startsWith(a,"--steps=")) selected <- strsplit(sub("^--steps=","",a),",")[[1]]
  else stop("Unknown argument; use --dry-run, --root=..., or --steps=01,02,...")
}
steps <- read.csv(file.path(SOURCE_ROOT,"config/steps.csv"),colClasses="character")
if (is.null(selected)) selected <- steps$step
if (!all(selected %in% steps$step)) stop("Unknown step identifier.")
plan <- steps[steps$step %in% selected,,drop=FALSE]
inputs <- read.csv(file.path(SOURCE_ROOT,"config/inputs.csv"),colClasses="character")
cat("Selected analysis steps:\n")
for (j in seq_len(nrow(plan))) cat(plan$step[j],plan$description[j],"\n")
need <- vapply(strsplit(inputs$steps,";",fixed=TRUE),function(v) any(v %in% selected),logical(1))
missing <- inputs$file[need & !file.exists(file.path(PROJECT_ROOT,inputs$file))]
if (length(missing)) cat("Missing public input files:\n",paste(missing,collapse="\n"),"\n")
if (dry) {
  cat("Dry run complete; no data downloaded or analyses executed.\n")
  quit(status=0)
}
if (length(missing)) stop("Supply public inputs listed above; see docs/DATA_SOURCES.md.")
OUTPUT_DIR <- file.path(PROJECT_ROOT,"outputs")
MATRIX_DIR <- file.path(OUTPUT_DIR,"fullmatrix")
RESULT_DIR <- file.path(OUTPUT_DIR,"analysis")
for (p in c(OUTPUT_DIR,MATRIX_DIR,RESULT_DIR)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
# Keep R-managed database caches inside the excluded project output directory.
cache <- file.path(OUTPUT_DIR,"cache")
dir.create(cache,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(R_USER_CACHE_DIR=cache)
packages <- read.delim(file.path(SOURCE_ROOT,"config/package_versions.tsv"),stringsAsFactors=FALSE)
required <- packages$package
absent <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(absent)) stop(paste("Missing packages:",paste(absent,collapse=", "),"; see install_dependencies.R"))
versions <- data.frame(package=required,version=vapply(required,function(p)as.character(packageVersion(p)),character(1)))
write.table(versions,file.path(OUTPUT_DIR,"package_versions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(c(paste("R",getRversion()),paste("Platform",R.version$platform)),file.path(OUTPUT_DIR,"runtime.txt"))
provenance <- inputs[need,c("file","source"),drop=FALSE]
provenance$md5 <- unname(tools::md5sum(file.path(PROJECT_ROOT,provenance$file)))
write.csv(provenance,file.path(OUTPUT_DIR,"input_checksums.csv"),row.names=FALSE)
common <- new.env(parent=globalenv())
for (n in c("SOURCE_ROOT","PROJECT_ROOT","OUTPUT_DIR","MATRIX_DIR","RESULT_DIR")) assign(n,get(n),envir=common)
sys.source(file.path(SOURCE_ROOT,"R/common.R"),envir=common,keep.source=FALSE)
for (j in seq_len(nrow(plan))) {
  cat("\nRunning",plan$step[j],plan$description[j],"\n")
  env <- new.env(parent=common)
  sys.source(file.path(SOURCE_ROOT,plan$script[j]),envir=env,keep.source=FALSE)
}
cat("\nSelected steps completed. Outputs remain excluded from the source distribution.\n")
