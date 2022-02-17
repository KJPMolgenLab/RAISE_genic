
home=getwd()
args = commandArgs(trailingOnly=TRUE)

##args=c("/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Merged_LGD_RAISE_Genic_WES_F_K_B_merged.vcf.gz.raw.part_ab", "/scratch/fuchs/agchiocchetti/chiocchetti/RAISEGenicData/VCF_Files/Merged_LGD_RAISE_Genic_WES_F_K_B_merged.vcf.gz.bim", "part_ab")

args=c("S:/KJP_Biolabor/Projects/RAISE_GENIC/data/LGDSData.raw",
       "S:/KJP_Biolabor/Projects/RAISE_GENIC/data/LGDSData.bim", "Local_res")


print(args)

urltoLGDs=args[1]
filetarget.bim=args[2]
outputname=args[3]

output= paste0(home, "/output")
hg38refgenes="S:/KJP_Biolabor/RawData/99_Datasets/RefSeqGeneshg38TxStartStop.txt"

print("loading data")

labels=read.csv("S:/KJP_Biolabor/Projects/RAISE_GENIC/output/aed_responses.csv")

labels= labels %>% column_to_rownames("record_id")



rm(TOM)
gc()
gc()
gc()

print("loading pacakges")
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
library(xgboost)
source(paste0(home,"/code/custom_functions.R"))


# cat("number of cores detected", detectCores() , "\n")
#
# cl=makeCluster(detectCores()-2, type="PSOCK")
# registerDoParallel(cl)

# print(paste0("running on ", getDoParWorkers(), " cores"))

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
  gc()
  return(res)
}



gethitgenes=function(patID="FINKPH_EPI-BIO_AA0268A", genofile=LGDS.gen, bimfile=LGDs.bim){
  hgnchit = paste0(bimfile[,"gene"][which(genofile[patID,]>0)], collapse = ";") %>% gsub(";;", ";", .)
  hgnchit = unique(unlist(strsplit(hgnchit, ";")))
  return(hgnchit)
}




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


if(FALSE){
  print("calculate features")
  load(paste0(output,"/WGCNA_adj_TOM.RData"))
  load(paste0(output,"/dds_matrix.RData"))

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


  reslist=list()
  clustergene=mcols(ddsMat)$cluster
  allgenes=mcols(ddsMat)$hgnc

  ## sort soit starts with the biggest


  modules=sort(table(clustergene), decreasing=F) %>% names
  modules=modules[1:22]
  modules=c("black", "midnightblue", "grey60")

  reslist=lapply(rownames(LGDs), get_pat_parameters)
  names(reslist)=rownames(LGDs)[]

  print("DONE calculate features")

  save(reslist, file=paste0(output, "/ReslistWGCNA_pat", outputname, ".RData"))
} else {
  print("load features")

  load(paste0(output, "/ReslistWGCNA_pat", outputname, ".RData"))}


NetwFeatures_list=lapply(reslist, function(x){unlist(x)})
NetwFeatures=do.call(rbind, NetwFeatures_list)
labels_ordered = labels[rownames(NetwFeatures),]

write.table(NetwFeatures, file=paste0(output, "/ReslistWGCNA_pat",outputname,".txt"))

print("DONE write features")

gc()

# if no variance featurtes present add noise

NetwFeatures=t(NetwFeatures) %>%  as.data.frame()

noiseridx=which(apply(NetwFeatures, 1, sd)==0) %>% names()

NetwFeatures[noiseridx,] = matrix(data=unlist(NetwFeatures[noiseridx,])+rnorm(length(unlist(NetwFeatures[noiseridx,])), 0,0.1),
                                  nrow=dim(NetwFeatures[noiseridx,])[1], ncol=dim(NetwFeatures[noiseridx,])[2])

#NetwFeatures=t(NetwFeatures) %>%  as.data.frame()
PCA=prcomp(scale(NetwFeatures))
res_umap=M3C::umap(NetwFeatures, colvec = rownames(PCA$rotation))
res_tsne=M3C::tsne(NetwFeatures, labels = rownames(PCA$rotation))


pdf(paste0(output,"Patient_module_properties.pdf"))
par(mfrow=c(1,3))
plot(PCA$rotation[,"PC1"], PCA$rotation[,"PC2"],
     main="PCA of WGCNA modules",
     xlab="PC1", ylab="PC2",
     col=labels_ordered$Valproate_28, pch=16)
plot(res_umap$data$X1,res_umap$data$X2,
     main="UMAP of WGCNA modules",
     xlab="UMAP 1", ylab="UMAP 2",
     col=labels_ordered$Valproate_28,pch=16)
plot(res_tsne$data$X1,res_tsne$data$X2,
     main="t-sne of WGCNA modules",
     xlab="t-sne 1", ylab="t-sne 2",
     col=labels_ordered$Valproate_28,pch=16)

dev.off()



scaled_data = scale(t(NetwFeatures))



dataset = as.data.frame(scaled_data)
dataset$y=as.factor(labels_ordered$Levetiracetam_13 %>%  as.factor())

dataset_filt = dataset %>% filter( complete.cases(.))

trainset=dataset_filt[sample(rownames(dataset_filt), 0.7*nrow(dataset_filt)), ]
testset = dataset_filt[!rownames(trainset) %in% dataset_filt,]


# balance train set


maxlab=sum(table(trainset$y) %>% which.max() == trainset$y)

amplified=trainset[table(trainset$y) %>% which.max() == trainset$y,]

for( i in levels(trainset$y)){
  idx=trainset$y==i
  if(sum(idx)<maxlab){
    ridx=sample(rownames(testset)[idx],maxlab-sum(trainset$y==i), replace=T)
    if(FALSE){set=matrix(data=unlist(testset[ridx,])+
                           rnorm(length(unlist(testset[ridx,])),0,0.1),
                         nrow=dim(testset[ridx,])[1],
                         ncol=dim(testset[ridx,])[2]) %>%  as.data.frame()
    } else {
      set=testset[ridx,]
    }
    colnames(set)=colnames(testset)
    set$y=i
    rownames(set) = paste0("amp_",make.names(ridx, unique = T))
  }
  amplified=rbind(amplified,set)
}

amp_trainset=amplified

paramgrid=expand_grid(
  maxdepth=c(2:5),
  eta=seq(0,1,length.out=5)[2:5]
  )

for(i in 1:nrow(paramgrid))
{
  param <- list(max_depth = paramgrid$maxdepth[i], eta = paramgrid$eta[i],
                verbose = 0, nthread = 2,
                objective = "multi:softprob", num_class=4)

  rf = xgboost(data=as.matrix(amp_trainset %>% select(-"y")),
               label = as.numeric(amp_trainset$y)-1,
               nrounds =300,params=param)

  pred=predict(rf, as.matrix(as.matrix(testset %>% select(-"y"))))
  pred=matrix(data=pred,nrow = 413, ncol=4, byrow = T)

  pred_y = apply(pred,1, which.max)

  mat=as.matrix(table(pred_y, testset$y))

  accuracy <- sum(diag(mat)) / sum(mat)
  paramgrid[i,"accuracy"] = accuracy
}


param <- list(max_depth = paramgrid$maxdepth[which.max(paramgrid$accuracy)],
              eta = paramgrid$eta[which.max(paramgrid$accuracy)],
              verbose = 0, nthread = 2,
              objective = "multi:softprob", num_class=4)

rf = xgboost(data=as.matrix(amp_trainset %>% select(-"y")),
             label = as.numeric(amp_trainset$y)-1,
             nrounds =500,params=param)


pred=predict(rf, as.matrix(as.matrix(testset %>% select(-"y"))))
pred=matrix(data=pred,nrow = 413, ncol=4, byrow = T)

pred_y = apply(pred,1, which.max)

mat=as.matrix(table(pred_y, testset$y))

accuracy <- sum(diag(mat)) / sum(mat)
precision <- diag(mat) / rowSums(mat)
recall <- (diag(mat) / colSums(mat))
paramgrid[i,"accuracy"] = accuracy

cat("accuracy: ", accuracy, "\n",
    "precision: ", precision, "\n",
    "recall: ", recall, "\n")


########################

param <- list(max_depth = paramgrid$maxdepth[which.max(paramgrid$accuracy)],
              eta = paramgrid$eta[which.max(paramgrid$accuracy)],
              verbose = 0, nthread = 2,
              objective = "binary:logistic")

table(trainset$y==1)

rf = xgboost(data=as.matrix(trainset %>% select(-"y")),
             label = as.numeric(trainset$y==1),
             nrounds =500,params=param)


pred=predict(rf, as.matrix(as.matrix(testset %>% select(-"y"))))


mat=as.matrix(table(predicted=pred>0.9, observed=testset$y==1))

accuracy <- sum(diag(mat)) / sum(mat)
precision <- diag(mat) / rowSums(mat)
recall <- (diag(mat) / colSums(mat))
paramgrid[i,"accuracy"] = accuracy

cat("accuracy: ", accuracy, "\n",
    "precision: ", precision, "\n",
    "recall: ", recall, "\n")

library(pROC)

pROC_obj <- roc(testset$y==1,pred,
               smoothed = TRUE,
               # arguments for ci
               ci=TRUE, ci.alpha=0.9, stratified=FALSE,
               # arguments for plot
               plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
               print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")

predobj
