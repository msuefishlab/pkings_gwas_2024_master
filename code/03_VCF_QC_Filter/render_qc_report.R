#!/usr/bin/env Rscript


#run this script like this:
# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd

args = commandArgs(trailingOnly=TRUE)

root<-rprojroot::find_root(".git/index")
require(yaml)
require(tools)


rmarkdown::render(file.path(args[3]),
                  params=list(
                  data_path=args[1],
                  output_file=file.path(paste0(args[2], ".report.html"))))
