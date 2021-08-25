###############################################
#						#
# 	01) Mapping And Alignment		#
#						#
###############################################


args <- commandArgs(trailingOnly = TRUE)
folder = args[1]

mode = "reanalyze"  # test only creates subset, "reanalyze" will calculate all BAM  files  again

fastqcfiles = folder
refhgfolder = "/scratch/fuchs/agmisc/chiocchetti/ReferenceGenomes/"
reffile = "hg38.fa.gz"
output=paste0(folder, "/Countmatrix/")


dir.create(output)

require("parallel")
require("Rsubread")

ncores = detectCores()


##Build Reference genome index

setwd(refhgfolder)
if(! file.exists("hg38.reads")){
    buildindex(basename=paste0(refhgfolder,"hg38"),reference=reffile)
}

setwd(output)
##Read in the fastq files
fastq.files<-list.files(path=fastqcfiles, pattern=".fq$",full.names=TRUE)

if (mode=="test")
    fastq.files=fastq.files[1:5]


if (mode == "reanalyze")
    align(index=paste0(refhgfolder,"hg38"),readfile1=fastq.files, nthreads = ncores)

if (mode=="normal"){
    chckfiles=paste0(fastq.files,".subread.BAM.summary")
    filestomake=setdiff(chckfiles, list.files(fastqcfiles,full.names = TRUE))
    filestomake=gsub(".subread.BAM.summary","",filestomake)
    align(index=paste0(refhgfolder,"hg38"),readfile1=filestomake, nthreads = ncores)
}


##Reading in Bam files
bam.files <- list.files(path =fastqcfiles, pattern = ".BAM$", full.names = TRUE)

#Getting the feature counts
fc <- featureCounts(bam.files, annot.inbuilt="hg38",countMultiMappingReads=TRUE,
                    nthreads = ncores)

#Deleting extra information from column names
Countdata <- data.frame(fc$counts)
colnames(Countdata)<-sub("_.*", "", colnames(Countdata))

save(Countdata, file=paste0(output,"Countmatrix.RData"))
