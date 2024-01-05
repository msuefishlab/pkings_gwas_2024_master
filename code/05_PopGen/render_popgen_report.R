#!/usr/bin/env Rscript


#run this script like this:
# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd

args = commandArgs(trailingOnly = TRUE)

root <- rprojroot::find_root(".git/index")
require(yaml)
require(tools)


rmarkdown::render(
  input = file.path(args[3]),
  output_file = basename(paste0(args[2], ".report.html")),
  output_dir = dirname(args[2]),
  params = list(
    popgen_prefix = args[1]
  )
)