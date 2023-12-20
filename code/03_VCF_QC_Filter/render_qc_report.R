#!/usr/bin/env Rscript


#run this script like this:
# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd

args = commandArgs(trailingOnly=TRUE)

root<-rprojroot::find_root(".git/index")
require(yaml)
require(tools)


rmarkdown::render(file.path(args[3]),
                  params=list(
                    data_path=file.path(root,'results/01_QC/',args[1])),
                    output_file=file.path(root,'results/01_QC',paste0(args[2],'.report.html')))
