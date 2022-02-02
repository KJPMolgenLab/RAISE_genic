
home=getwd()
args = commandArgs(trailingOnly=TRUE)

##args=c("/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Merged_LGD_RAISE_Genic_WES_F_K_B_merged.vcf.gz.raw.part_aa", "/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Merged_LGD_RAISE_Genic_WES_F_K_B_merged.vcf.gz.bim", "part_aa")

print(args)

urltoLGDs=args[1]
filetarget.bim=args[2]
outputname=args[3]

output= "/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/output"
hg38refgenes="/scratch/fuchs/agchiocchetti/chiocchetti/ReferenceGenomes/RefSeqGeneshg38TxStartStop.txt"

source(paste0(home,"/code/custom_functions.R"))
load(paste0(output,"/WGCNA_adj_TOM.RData"))
load(paste0(output,"/dds_matrix.RData"))

require(M3C)
library(igraph)
library(DESeq2)
library(RCurl)
library(GenomicRanges)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(moments)

cat("number of cores detected", detectCores(), "\n")
registerDoParallel(cores=detectCores())

print(paste0("running on ", getDoParWorkers(), " cores"))


print("define functions")
get_network_params=function(modulename="", clusters, cormatrix, split=F, cores=NA, splitsize=500) {
    res=c()
    module=modulename
    index=which(clusters==module)
    modMatrixorig=cormatrix[index,index]
    modMatrixorigdist=1-cormatrix[index,index]
    
    diag(modMatrixorig) = NA
    diag(modMatrixorigdist)=NA
    modMatrix=modMatrixorig
    graph=graph_from_adjacency_matrix(modMatrix, weighted = T)
    graphdist=graph_from_adjacency_matrix(modMatrixorigdist, weighted = T)
                                        #calculate the parameters
    
                                        #print("... wConSum")
    res[["wConSum"]]=sum(modMatrixorig[lower.tri(modMatrixorig, diag = F)])
                                        #print("... wConMean")
    res[["wConMean"]]=mean(modMatrixorig[lower.tri(modMatrixorig, diag = F)])
                                        #print("... wAvgPath")
    if(length(V(graphdist))>splitsize & split){
        nc=min(c(getDoParWorkers(), cores), na.rm=T)
        cat("......", nc, " cores used\n")
        size=max(c(ceiling(length(V(graphdist))/nc), splitsize), na.rm=T)
        cat("......", size, " genes per batch\n")
        sets=split(V(graphdist), rep(1:ceiling(length(V(graphdist))/size), each=size)[1:length(V(graphdist))])
        cat("......", length(sets), " batches\n")
        dist <-  foreach(s=sets, .combine = 'c') %dopar% 
            igraph::distances(graphdist, v=s, mode = "all", algorithm = "dijkstra")
    } else {
        dist = igraph::distances(graphdist, mode = "all", algorithm = "dijkstra")
    }
    res[["wAwgPath"]] = mean(c(dist))
    res[["wKurtPath"]] = kurtosis(c(dist))
    res[["wSkewPath"]] = skewness(c(dist))
                                        #print("... wDistance")
    res[["wDistance"]]=farthest.nodes(graph)$distance
    
                                        #print("... wTransitivity")
                                        #make binary
                                        #print("... cutoff binarization")
    cutoff=1
    sum_of_noconn=ncol(modMatrixorig)
    
    plottab=data.frame(cutoff=cutoff, unconnected=sum_of_noconn)
    counter=0
    while(sum_of_noconn>0.05*ncol(modMatrixorig) & counter<1e5){
        cutoff=cutoff/1.1
        modMatrix = modMatrixorig
        modMatrix[modMatrix<cutoff]=NA
        diag(modMatrix)=NA
        modMatrix[! is.na(modMatrix)] = 1
        sum_of_noconn=sum(colSums(modMatrix, na.rm=T)==0)
        plottab=rbind(plottab,c(cutoff, sum_of_noconn))
        counter=counter+1
    }

    res[["uwBinCutoff"]]=cutoff
    
    graph=graph_from_adjacency_matrix(modMatrix)
                                        #print("... uwMaxcon")
    res[["uwMaxcon"]]=max(colSums(modMatrix, na.rm=T))
                                        #print("... uwConSum")
    res[["uwConSum"]]=sum(modMatrix[lower.tri(modMatrix, diag = F)], na.rm=T)
                                        #print("... uwDensity")
    res[["uwDensity"]]=graph.density(graph)
                                        #print("... uwDistance")
    res[["uwDistance"]]=farthest.nodes(graph)$distance
                                        #print("... uwMeanBetweeness")
    res[["uwMeanBetweeness"]]=mean(edge.betweenness(graph))
                                        #print("... uwTransitivity")
    res[["uwTransitivity"]]=transitivity(graph, type="global")
                                        #print("... uwEgoSize")
    res[["uwEgoSize"]]=mean(ego_size(graph))
    return(res)
}



print("loads LGD files")

LGDs=fread(urltoLGDs, data.table=F)
LGDs.bim=read.table(filetarget.bim, header=T)
colnames(LGDs.bim) = c("CHR", "rs", "cm", "pos", "alt","ref")
rownames(LGDs) = LGDs$IID
LGDS.gen = LGDs[,-(1:6)]

print("set NA to zero")
LGDS.gen[is.na(LGDS.gen)] <- 0


print("set colnames")

genes=read.table(hg38refgenes)
colnames(genes)=c("bin", "id", "CHR", "strand", "start", "stop", "hgnc", "id2")
genesGR=with(genes, GRanges(CHR, IRanges(start, stop),  ))
mcols(genesGR) = genes[,c("id", "hgnc")]

LGDs.bim.GR=with(LGDs.bim, GRanges(paste0("chr",CHR), IRanges(pos, pos)))
mcols(LGDs.bim.GR) = LGDs.bim[,c("rs", "alt", "ref")]


print("find hits")
hits=findOverlaps(LGDs.bim.GR, genesGR)
hits=as.data.frame(hits)
hits$hgnc=mcols(genesGR)$hgnc[hits$subjectHits]

hitgenes=tapply(hits$hgnc, hits$queryHits, function(x){paste0(unique(x),collapse=";")})
LGDs.bim$gene[as.numeric(names(hitgenes))] = hitgenes

gethitgenes=function(patID="FINKPH_EPI-BIO_AA0268A", genofile=LGDS.gen, bimfile=LGDs.bim){
  hgnchit = paste0(bimfile[,"gene"][which(genofile[patID,]>0)], collapse = ";") %>% gsub(";;", ";", .)
  hgnchit = unique(unlist(strsplit(hgnchit, ";")))
  return(hgnchit)
  }

reslist=list()
clustergene=mcols(ddsMat)$cluster
allgenes=mcols(ddsMat)$hgnc

## sort soit starts with the biggest 
modules=sort(table(clustergene), decreasing=T) %>% names

get_pat_parameters <- function(pat, verbose = T){
    cat(pat, "\n")
    modfeatures=list()
    hitcluster=clustergene
    hitcluster[allgenes %in% gethitgenes(pat)]=NA
    
    for (m in modules){
        if(verbose){cat("...",m, "\n")}
        tmp=get_network_params(modulename = m, 
                               clusters =hitcluster , 
                               cormatrix = adj,
                               splitsize=500, cores=10, split=F)
        modfeatures[[m]]=unlist(tmp)
    }
    return(modfeatures)
}



print("calculate features")

reslist=foreach(p=rownames(LGDs)[]) %dopar% get_pat_parameters(p)
save(reslist, file=paste0(output, "/ReslistWGCNA_pat", outputname, ".RData"))

NetwFeatures=sapply(dimnames(reslist)[[2]], function(x){unlist(reslist[,x])})
write.table(NetwFeatures, file=paste0(output, "/ReslistWGCNA_pat",outputname,".txt"))


if(FALSE){
    noiseridx=apply(NetwFeatures, 1, sd)==0

    NetwFeatures[noiseridx,] = matrix(data=rnorm(length(unlist(NetwFeatures[noiseridx,])), 0,0.1), 
                                      nrow=dim(NetwFeatures[noiseridx,])[1], ncol=dim(NetwFeatures[noiseridx,])[2])

    PCA=prcomp(t(scale(t(NetwFeatures))))
    res_umap=M3C::umap(NetwFeatures, colvec = rownames(PCA$rotation))

    pdf(paste0(output,"Patient_module_properties.pdf"))
    par(mfrow=c(1,2))
    plot(PCA$rotation[,"PC1"], PCA$rotation[,"PC2"], 
         main="PCA of WGCNA modules",
         xlab="PC1", ylab="PC2",
         col="black", pch=16)
    plot(res_umap$data$X1,res_umap$data$X2, 
         main="UMAP of WGCNA modules",
         xlab="UMAP 1", ylab="UMAP 2",
         col="black",pch=16)
    dev.off()
}


