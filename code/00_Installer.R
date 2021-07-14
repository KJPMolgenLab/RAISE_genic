### installer


package<-c("RColorBrewer", "dplyr", "httr", "parallel", "rJava",
           "compareGroups", "kableExtra", "tidyverse","dendextend",
           "knitr", "pheatmap", "glmpca", "plotly", "rsconnect", "xlsx",
           "lavaan", "semPlot", "ggplotify", "qqman", "webshot")


if(length(setdiff(package, rownames(installed.packages()))) > 0)	{
  install.packages(setdiff(package, rownames(installed.packages())),
                   INSTALL_opts = c('--no-lock'),dependencies=T,
                   ask=FALSE, repos="https://cloud.r-project.org/")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org/")

bpackage = c( "annotate","biomaRt", "DESeq2", "WGCNA", "edgeR","RRHO",
              "limma","Rsubread", "flashClust", "gprofiler2", "RCircos",
              "sva", "vsn", "EnhancedVolcano", "bumphunter",
              "TxDb.Hsapiens.UCSC.hg19.knownGene", "GenomicRanges", "Gviz",
              "rtracklayer")

if(length(setdiff(bpackage, rownames(installed.packages()))) > 0)	{
  BiocManager::install(setdiff(bpackage, rownames(installed.packages())))
}
