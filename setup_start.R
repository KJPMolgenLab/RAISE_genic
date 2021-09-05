library(workflowr)


setwd("S:/KJP_Biolabor/00_KJP_Biolabor/Projects/RAISE_GENIC/")
wflow_start("./", existing = TRUE, git=T)

wflow_use_github(organization = "KJPMolgenLab",
                 repository = "RAISE_genic")

source("code/00_Installer.R")

system("git add docs/* analysis/* code/*")
system("git commit -a -m \"adds QC and WGCNA\"")
system("git push -u origin master")
