library(workflowr)


setwd("S:/KJP_Biolabor/00_KJP_Biolabor/Projects/RAISE_GENIC/")

wflow_start("./", existing = TRUE, git=T)
wflow_use_github(organization = "KJPMolgenLab",
                 repository = "RAISE_genic")

source("code/00_Installer.R")

wflow_build("analysis/index.Rmd")
wflow_build("analysis/07_Visualize_phenodata.Rmd")

system.time(wflow_publish(c("analysis/*","docs/*","code/*")))
system("git push -u origin master")

system("git add docs/* analysis/* code/*")
system("git commit -a -m \"adds codebook and phenotype analysis\"")
system("git push -u origin master")

