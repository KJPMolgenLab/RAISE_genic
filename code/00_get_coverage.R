library(vcfR)

vcffileFFM="/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Frankfurt/DEUPUM_RAISE-GENIC_210526.vcf.gz"
vcffileKPO="/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Kuopio/AnVIL_CCDG_Broad_NP_Epilepsy_FINKPH_EPIL_CO_MORBIDI_MDS_WES_y2-only-raisegenic.vcf.gz"
vcffileBRS="/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Brussels/AnVIL_CCDG_Broad_NP_Epilepsy_BELULB_DS_EP_NPU_WES.vcf.gz"



vcf <- read.vcfR(vcffileFFM, verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
dp_median=apply(dp,2, function(x){median(x[x!=0],na.rm=T)})
dp_80=apply(dp,2, function(x){quantile(x[x!=0],0.2, na.rm=T)})
cat("FFM:DP_50 ",mean(dp_median), "\n")
cat("FFM:DP_80",mean(dp_80), "\n")

vcf <- read.vcfR(vcffileKPO, verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
dp_median=apply(dp,2, function(x){median(x[x!=0],na.rm=T)})
dp_80=apply(dp,2, function(x){quantile(x[x!=0],0.2, na.rm=T)})
cat("KPO:DP_50 ",mean(dp_median), "\n")
cat("KPO:DP_80",mean(dp_80), "\n")

vcf <- read.vcfR(vcffileBRS, verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
dp_median=apply(dp,2, function(x){median(x[x!=0],na.rm=T)})
dp_80=apply(dp,2, function(x){quantile(x[x!=0],0.2, na.rm=T)})
cat("BRS:DP_50 ",mean(dp_median), "\n")
cat("FFM:DP_80",mean(dp_80), "\n")



