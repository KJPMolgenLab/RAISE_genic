
library(data.table)

homedir="/scratch/fuchs/agmisc/chiocchetti/RAISEGenicData/VCF_Files/"
file="Merged_LGD_testset.vcf.gz"
gene_boudnaries =
    
    
setwd(homedir)

system(paste0("plink --vcf ", file, " --double-id --export A  --out ", file))

genotypes=fread(paste0(file,".raw"), header = T)
SNV_meta=fread(paste0(file,".bim"), header = T)

