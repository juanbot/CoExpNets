

print.bootnet = function(net){
  cat("A bootstrapped network created with mode",net$mode,"\n",
      "Soft thresholding parameter (beta):",net$beta,"\n",
      "Adjacency summary\n")
  print(summary(net$adjacency))
  cat("Final number of modules",length(unique(net$moduleColors)),"\n")
  cat("All samples net final number of modules",
      length(unique(net$allsamplesnet$moduleColors)),"\n")

  rdiffs = NULL
  fdiffs = NULL
  adiffs = NULL
  clusters = NULL
  for(i in 1:length(net$subnets)){
    if(!is.null(net$subnets[[i]]$subcluster)){
      clusters = rbind(clusters,net$subnets[[i]]$subcluster)
    }
  }

  for(i in 2:nrow(clusters)){
    rdiffs = c(rdiffs,mclust::adjustedRandIndex(clusters[i,],clusters[i-1,]))
    fdiffs = c(fdiffs,mclust::adjustedRandIndex(clusters[i-1,],net$moduleColors))
    adiffs = c(adiffs,mclust::adjustedRandIndex(clusters[i,],net$allsamplesnet$moduleColors))

  }
  cat("Successive Rand simmilarities\n")
  print(rdiffs)
  cat("Rand differences with all samples net\n")
  print(adiffs)
  cat("Rand differences with final net\n")
  print(fdiffs)
}

#' plotMDS - Testing gene sepparability with multi-dimensional
#' scalling.
#'
#' @param rpkms.net A dataframe with the expression data with samples in columns.
#' Columns are assumed to use sample IDs as column names. These IDs should appear
#' as row names in the covariate file to properly index each value
#' @param path Folder to write the pdf plot generated
#' @param covvars The variables from the covs parameter to use in the plot
#' @param label A string for using at the plot, informative purposes
#' @param n.mds Number of genes to use, randomly sampled
#' @param covs A data frame with all covariates for the samples. The rows must
#' be named with the sample IDs.
#'
#' @return
#' @export
#'
#' @examples
plotMDS = function(rpkms.net,path,covs,covvars,label,n.mds=-1){

  intersect.s = intersect(rownames(covs),colnames(rpkms.net))
  covs = covs[intersect.s,]
  rpkms.net = rpkms.net[,intersect.s]
  stopifnot(identical(colnames(rpkms.net),rownames(covs)))

  #Each covariate has a name equal to a column name
  #This plot is useful to see whether the covariate has any
  #discriminatory power between their different values or, on
  #the contrary, all samples appear mixed across values
  for(covvar in covvars){
    cat(paste0("Working on MDS plot for ",covvar,"\n"))
    pdf(paste0(path,covvar,".pdf"),height=8,width=30)
    if(n.mds > 0)
      mask = sample(ncol(rpkms.net),n.mds)
    else
      mask = c(1:ncol(rpkms.net))
    colors = rainbow(length(levels(covs[,covvar])))
    plotMDS(rpkms.net[,mask],col=colors[as.numeric(covs[,covvar])],
            main=paste0("MDS using ",covvar," ",label))
    legend("topright",fill=colors,
           legend=levels(covs[,covvar]))
    dev.off()
  }

}

#' Title
#'
#' @param mode
#' @param expr.data
#' @param n.iterations
#' @param job.path
#' @param allsampsnet
#' @param each
#' @param tissue
#' @param b
#' @param ...
#' @param removeTOM
#' @param min.cluster.size
#' @param waitFor
#'
#' @return
#' @export
#'
#' @examples
getBootstrapNetworkCl = function(mode=c("leaveoneout","bootstrap"),
                                 expr.data,
                                 n.iterations=50,
                                 removeTOM=F,
                                 job.path,
                                 min.cluster.size=100,
                                 allsampsnet=F,
                                 each=1,
                                 tissue="Bootstrap",
                                 blockTOM=F,
                                 clParams=" -l nodes=1:nv ",
                                 clParamsPost=clParams,
                                 waitFor=24*3600,
                                 batchSize=1000,
                                 timePerBatch=1800,
                                 b=10,...){
  print(expr.data[1:5,1:5])
  if(typeof(expr.data) == "character")
    expr.data = readRDS(expr.data)

  library(sgefacilities)
  #Lets create indexes
  indexes = NULL
  if(mode == "leaveoneout"){
    lapply(1:nrow(expr.data),function(x){
      indexes <<- rbind(indexes,(1:nrow(expr.data))[-x])
    })
  }else if(mode == "bootstrap"){
    lapply(1:b,function(x){
      indexes <<- rbind(indexes,sample(1:nrow(expr.data),nrow(expr.data),replace=T))
    })
  }else stop(paste0("Mode",mode,"unknown\n"))

  ngenes = ncol(expr.data)
  count = 0
  allclusters = NULL
  allsubnets = NULL
  handlers = NULL
  maskd = NULL
  for(i in 1:nrow(indexes)){
    count = count + 1
    lexpr.data = expr.data[indexes[i,],]
    if(mode == "bootstrap"){
      maskd = duplicated(rownames(lexpr.data))
      while(sum(maskd)){

        rownames(lexpr.data)[maskd] = paste0(rownames(lexpr.data)[maskd],"_d")
        maskd = duplicated(rownames(lexpr.data))
      }
    }

    ltissue = paste0(tissue,"_b_",i)
    params = NULL
    params$its = n.iterations
    params$mode = mode
    params$outfolder = job.path
    params$datain = lexpr.data
    params$min.cluster.size=min.cluster.size
    params$save.tom =T
    params$tissue = ltissue
    params$maskd = maskd
    params$blockTOM = blockTOM
    params$fun = function(tissue,datain,mode,min.cluster.size,its,blockTOM,outfolder,save.tom,maskd){
      library(CoExpNets)
      genes = colnames(datain)
      sampids = rownames(datain)
      datain = log2(1 + as.matrix(datain))
      if(mode == "bootstrap"){
        maskdi = which(maskd)
        for(index in maskdi){
          datain[index,] = jitter(datain[index,])
        }
      }
      colnames(datain) = genes
      rownames(datain) = sampids

      net = getDownstreamNetwork(tissue=tissue,
                                 n.iterations=n.iterations,
                                 save.tom = save.tom,
                                 min.cluster.size=min.cluster.size,
                                 save.plots = F,
                                 excludeGrey = F,
                                 blockTOM=blockTOM,
                                 beta=-1,
                                 net.type = "signed",
                                 debug=F,
                                 expr.data=datain,
                                 job.path=outfolder)
      return(net)

    }
    job = launchJob(parameters = params,
                    clParams = clParams,
                    prefix=ltissue,
                    justPrepareForLaunch=T,
                    wd=job.path)
    handlers[[job$jobname]] = job
  }
  #Now we send the post job
  params = NULL
  ltissue = paste0(tissue,"_b_Final")
  params$handlers = handlers
  params$expr.data = expr.data
  params$removeTOM = removeTOM
  params$tissue = tissue
  params$allsampsnet = allsampsnet
  params$n.iterations = n.iterations
  params$min.cluster.size=min.cluster.size
  params$each = each
  params$indexes = indexes
  params$job.path = job.path
  params$mode = mode
  params$waitFor = waitFor
  params$blockTOM = blockTOM
  params$fun = postCluster
  singlehandler = launchJob(parameters = params,
                            clParams = clParamsPost,
                            justPrepareForLaunch=T,
                            prefix=ltissue,
                            wd=job.path)
  handf = paste0(job.path,"/",tissue,"_handlers.rds")

  #Lets launch first the sentinel
  sentinel = NULL
  sentinel[[singlehandler$jobname]] = singlehandler
  submitJobs(sentinel)
  #And now the rest of single networks
  submitJobs(handlers,batchSize=batchSize,timePerBatch=timePerBatch)

  #Now we store all handlers
  handlers[[singlehandler$jobname]] = singlehandler
  saveRDS(handlers,handf)
  return(handf)
}

postCluster = function(handlers,
                       expr.data,
                       tissue,
                       allsampsnet,
                       n.iterations,
                       min.cluster.size,
                       blockTOM=F,
                       removeTOM=F,
                       each=1,
                       mode,
                       waitFor=24*3600,
                       indexes,
                       job.path){

  time = Sys.time()
  ngenes = ncol(expr.data)
  TOM = matrix(nrow = ngenes,ncol=ngenes)
  TOM[] = 0
  allsubnets = NULL
  initHandlers = length(handlers)

  while(length(handlers) > 0 & (Sys.time() - time) < waitFor){

    cat("We will collect Jobs, still",length(handlers),"jobs to go\n")
    nets = waitForJobs(handlers=handlers,
                       timeLimit=20,
                       increment=20,
                       removeData=F,
                       removeLogs=F,
                       qstatworks=F,
                       wd=job.path)

    cat("Just came back from waitForJobs with",length(nets),"handlers\n")
    tomCount = 0
    for(jobname in names(nets)){
      net = nets[[jobname]]
      if(is.null(net)){
        cat("This work have not properly finished\n")
      }else{

        handlers[[jobname]] = NULL
        #print(net)
        net = net$result
        ftocheck = net$tom
        if(blockTOM)
          ftocheck = paste0(ftocheck,"_metadata.rds")

        if(file.exists(ftocheck)){

          #print(net)
          cat("Accumulating TOM",net$tom,"\n")
          if(blockTOM){
            cat("Reading TOM\n")
            TOM = TOM + CoExpNets::readTOM(net$tom)
            CoExpNets::removeTOM(net$tom)
            print("Done")
          }else{
            cat("Reading TOM\n")
            TOM = TOM + readTOM(net$tom)
            cat("Quantile normalization\n")
            tom = preprocessCore::normalize.quantiles(tom)
            print("Adding TOM")
            TOM = TOM + tom
            print("Done")
            rm("tom")
            if(removeTOM)
              file.remove(net$tom)
          }

          tomCount = tomCount + 1
          #file.remove(net$net)
          print("TOM updated\n")
          if((initHandlers - length(handlers)) %% each == 0 | length(handlers) == 0){
            dissTOM = 1 - TOM/tomCount

            geneTree = flashClust::flashClust(as.dist(dissTOM), method = "average")
            print("Now the genetree")
            n.mods = 0
            deep.split = 2
            while(n.mods < 10 & deep.split < 5){
              dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,
                                                          distM = dissTOM,
                                                          deepSplit = deep.split,
                                                          pamRespectsDendro = FALSE,
                                                          respectSmallClusters = F,
                                                          minClusterSize = min.cluster.size)
              n.mods = length(table(dynamicMods))
              deep.split = deep.split + 1
            }
            rm(dissTOM)
            # Convert numeric lables into colors
            #This will print the same, but using as label for modules the corresponding colors
            localnet = NULL
            localnet$moduleColors = CoExpNets::dropGreyModule(WGCNA::labels2colors(dynamicMods))
            print(table(localnet$moduleColors))

            outnet = CoExpNets::applyKMeans(tissue=tissue,
                                            n.iterations=n.iterations,
                                            net.file=localnet,
                                            expr.data=expr.data)

            net$subcluster = outnet$net$moduleColors

          }
          #net$indexes = indexes[tomCount,]
          allsubnets[[jobname]] = net
        }
      }
    }

  }

  tomCount = tomCount - 1
  finalnet = NULL
  finalnet$file = paste0(job.path,"/netBoot",tissue,".default.it.",n.iterations,".b.",nrow(indexes),".rds")

  if(length(allsubnets) > 0){
    finalnet$beta = as.integer(mean(unlist(lapply(allsubnets,function(x){return(x$beta)}))))
    finalnet$file = paste0(job.path,"/netBoot",tissue,".",finalnet$beta,".it.",n.iterations,".b.",nrow(indexes),".rds")

    TOM = TOM/tomCount
    adjacency = apply(TOM,2,sum)
    adjacency = adjacency/length(adjacency)
    names(adjacency) = colnames(expr.data)
    if(blockTOM){
      finalnet$tom = paste0(finalnet$file,".tom.")
      CoExpNets::saveTOM(TOM,outnet$net$moduleColors,finalnet$tom)
    }else{
      finalnet$tom = paste0(finalnet$file,".tom.rds")
      saveRDS(TOM,finalnet$tom)
    }


    finalnet$adjacency = adjacency
    finalnet$moduleColors = outnet$net$moduleColors
    finalnet$subnets = allsubnets
    finalnet$mode = mode
    rm("TOM")
    outnet = CoExpNets::applyKMeans(tissue=tissue,
                                    n.iterations=n.iterations,
                                    net.file=finalnet,
                                    expr.data=expr.data)
    finalnet$moduleColors = outnet$net$moduleColors
    finalnet$MEs = outnet$net$MEs
    finalnet$cgenes = outnet$cgenes
    finalnet$partitions = outnet$partitions
  }

  if(allsampsnet){
    finalnet$allsamplesnet = CoExpNets::getDownstreamNetwork(expr.data=expr.data,
                                                             tissue=paste0(tissue,"_allsamps_"),
                                                             job.path=job.path,
                                                             min.cluster.size = min.cluster.size,
                                                             save.plots=F,
                                                             excludeGrey=F,
                                                             net.type="signed",
                                                             debug=F,
                                                             n.iterations=n.iterations,
                                                             save.tom=F)$net
  }

  if(!is.null(finalnet)){
    attr(finalnet,"class") <- "bootnet"
    saveRDS(finalnet,finalnet$file)
    return(finalnet$file)
  }

  return(NULL)
}



#' Title
#'
#' @param mode
#' @param expr.data
#' @param n.iterations
#' @param job.path
#' @param allsampsnet
#' @param each
#' @param tissue
#' @param b
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
getBootstrapNetwork = function(mode=c("leaveoneout","bootstrap"),
                               expr.data,
                               n.iterations=50,
                               job.path,
                               allsampsnet=F,
                               excludeGrey=F,
                               each=1,
                               blockTOM=F,
                               min.cluster.size=100,
                               tissue="Bootstrap",
                               b=10,...){
  print(expr.data[1:5,1:5])
  if(typeof(expr.data) == "character")
    expr.data = readRDS(expr.data)

  #Lets create indexes
  indexes = NULL
  if(mode == "leaveoneout"){
    lapply(1:nrow(expr.data),function(x){
      indexes <<- rbind(indexes,(1:nrow(expr.data))[-x])
    })
  }else if(mode == "bootstrap"){
    lapply(1:b,function(x){
      indexes <<- rbind(indexes,sample(1:nrow(expr.data),nrow(expr.data),replace=T))
    })
  }else stop(paste0("Mode",mode,"unknown\n"))

  ngenes = ncol(expr.data)
  TOM = matrix(nrow = ngenes,ncol=ngenes)
  TOM[] = 0
  count = 0
  allclusters = NULL
  allsubnets = NULL

  for(i in 1:nrow(indexes)){
    count = count + 1
    lexpr.data = expr.data[indexes[i,],]
    if(mode == "bootstrap"){
      maskd = duplicated(rownames(lexpr.data))
      while(sum(maskd)){
        rownames(lexpr.data)[maskd] = paste0(rownames(lexpr.data)[maskd],"_d")
        maskd = duplicated(rownames(lexpr.data))
        #print(maskd)
      }

    }
    ltissue = paste0(tissue,"_b_",i)
    net = getDownstreamNetwork(expr.data=lexpr.data,
                               min.cluster.size=min.cluster.size,
                               tissue=ltissue,
                               job.path=job.path,
                               excludeGrey = excludeGrey,
                               blockTOM=blockTOM,
                               save.plots=F,
                               n.iterations=0,
                               save.tom=T,...)

    cat("Accumulating TOM",net$tom,"\n")
    cat("Reading TOM\n")
    if(!blockTOM){
      tom = readRDS(net$tom)
      file.remove(net$tom)
      cat("Quantile normalization\n")
      tom = preprocessCore::normalize.quantiles(tom)
    }else{
      tom = readTOM(net$tom)
      removeTOM(net$tom)
    }
    cat("Adding TOM\n")
    TOM = TOM + tom
    cat("Done\n")
    rm("tom")
    if(count %% each == 0 | count == nrow(indexes)){
      dissTOM = 1 - TOM/count
      geneTree = flashClust::flashClust(as.dist(dissTOM), method = "average")
      print("Now the genetree")
      n.mods = 0
      deep.split = 2
      while(n.mods < 10 & deep.split < 5){
        dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,
                                                    distM = dissTOM,
                                                    deepSplit = deep.split,
                                                    pamRespectsDendro = FALSE,
                                                    respectSmallClusters = F,
                                                    minClusterSize = min.cluster.size)
        n.mods = length(table(dynamicMods))
        deep.split = deep.split + 1
      }
      rm("dissTOM")

      # Convert numeric lables into colors
      #This will print the same, but using as label for modules the corresponding colors
      localnet = NULL
      localnet$moduleColors = CoExpNets::dropGreyModule(WGCNA::labels2colors(dynamicMods))
      outnet = applyKMeans(tissue=tissue,
                           n.iterations=n.iterations,
                           net.file=localnet,
                           expr.data=expr.data)

      net$subcluster = outnet$net$moduleColors

    }
    net$indexes = indexes[i,]
    allsubnets[[i]] = net

  }

  finalnet = NULL
  finalnet$beta = as.integer(mean(unlist(lapply(allsubnets,function(x){return(x$beta)}))))
  finalnet$file = paste0(job.path,"/netBoot",tissue,".",finalnet$beta,".it.",
                         n.iterations,".b.",nrow(indexes),".rds")
  finalnet$tom = paste0(finalnet$file,".tom.rds")
  TOM = TOM/count
  adjacency = apply(TOM,2,sum)
  adjacency = adjacency/length(adjacency)
  names(adjacency) = colnames(expr.data)
  saveRDS(TOM,finalnet$tom)
  finalnet$adjacency = adjacency
  finalnet$moduleColors = outnet$net$moduleColors
  finalnet$subnets = allsubnets

  outnet = applyKMeans(tissue=tissue,
                       n.iterations=n.iterations,
                       net.file=finalnet,
                       expr.data=expr.data)
  finalnet$moduleColors = outnet$net$moduleColors
  finalnet$MEs = outnet$net$MEs
  finalnet$mode = mode
  finalnet$cgenes = outnet$cgenes
  finalnet$partitions = outnet$partitions

  if(allsampsnet)
    finalnet$allsamplesnet = getDownstreamNetwork(expr.data=expr.data,
                                                  tissue=paste0(tissue,"_allsamps_"),
                                                  job.path=job.path,
                                                  save.plots=F,
                                                  min.cluster.size = min.cluster.size,
                                                  n.iterations=n.iterations,
                                                  save.tom=F)


  attr(finalnet,"class") <- "bootnet"
  saveRDS(finalnet,finalnet$file)
  return(finalnet$file)

}




#' getDownstreamNetwork - Create a network
#'
#' @param tissue A label to use to refer to the results in files and downstream results
#' @param n.iterations Number of iterations in the k-means algorithm to refine the
#' clustering process
#' @param expr.data A data frame with expression data, ready to get into WGCNA. Columns
#' are genes and rows are samples. Column names will be used as gene names in the
#' analysis product
#' @param beta The smoothing parameter of the WGCNA pipeline. If -1, then the method
#' will suggest one automatically by looking at the R2 between gene connectivity and
#' Scale Free Topology features
#' @param job.path This method will generate a number of files. It needs a folder to write
#' results
#' @param min.cluster.size Minimum number of genes for a group to be considered as cluster
#' for the tree cutting algorithm to convert from a dendogram to a cluster.
#' @param net.type Whether a signed ("signed") or unsigned ("unsigned") network type will be created
#' @param debug Set this to true if you want a quick run of the method to test how it works with
#' a small amount of your genes
#' @param excludeGrey If WGCNA detects grey genes, set it to TRUE if you want them removed
#' from the network before applying k-means
#'
#' @return A file name that can be used to access your network
#' @export
#'
#' @examples
getDownstreamNetwork = function(tissue="mytissue",
                                n.iterations=50,	#Number of iterations for k-means, 50 recommended
                                min.exchanged.genes=20,
                                expr.data, 	#We expect a file name pointing to a dataframe (RDS format) with
                                #genes in columns and samples in rows. Each gene name appears
                                #in the column name. Better to use gene symbols as names
                                beta=-1,	#If -1 the algorithm will seek for the best beta
                                job.path="~/tmp/",	#Where to store all results
                                min.cluster.size=30,		#Minimum number of genes to form a cluster
                                net.type="signed",			#Leave it like that (see WGCNA docs)
                                debug=F,
                                blockTOM=F,
                                save.tom=F,
                                save.plots=F,
                                excludeGrey=FALSE){

  final.net=NULL
  distance.type="cor"
  centroid.type="pca"
  cor.type="pearson"

  if(debug){
    if(typeof(expr.data) == "character")
      expr.data = readRDS(expr.data)
    expr.data = expr.data[,1:1500]
    n.iterations=5
  }

  net.and.tom = getAndPlotNetworkLong(expr.data=expr.data,
                                      beta=beta,
                                      tissue.name=tissue,
                                      min.cluster.size=min.cluster.size,
                                      save.plots=save.plots,
                                      excludeGrey=excludeGrey,
                                      additional.prefix=job.path,
                                      return.tom=T,
                                      cor.type=cor.type)

  if(is.null(final.net))
    final.net = paste0(job.path,"/","net",tissue,".",
                       net.and.tom$net$beta,".it.",n.iterations,".rds")

  outnet = CoExpNets::applyKMeans(tissue=tissue,
                                  n.iterations=n.iterations,
                                  net.file=net.and.tom$net,
                                  expr.data=expr.data,
                                  excludeGrey=excludeGrey,
                                  min.exchanged.genes = min.exchanged.genes)

  if(save.tom){
    if(blockTOM)
      saveTOM(tom=net.and.tom$tom,
              clusters=outnet$net$moduleColors,
              filepref=paste0(final.net,".tom."))
    else
      saveRDS(net.and.tom$tom,paste0(final.net,".tom.rds"))
  }

  foutnet = NULL
  foutnet$beta = net.and.tom$net$beta
  foutnet$type = net.and.tom$net$type
  foutnet$net = outnet$net
  foutnet$partitions = outnet$partitions
  foutnet$cgenes = outnet$cgenes

  if(save.tom){
    if(blockTOM)
      foutnet$tom = paste0(final.net,".tom.")
    else
      foutnet$tom = paste0(final.net,".tom.rds")
  }



  if(save.plots){
    cat("Generating mod sizes for",final.net,"\n")
    pdf(paste0(final.net,".mod_size.pdf"))
    plotModSizes(which.one="new",tissue=final.net)
    dev.off()
    pdf(paste0(final.net,".Eigengenes_clustering.pdf"))
    plotEGClustering(which.one="new",tissue=final.net)
    dev.off()
  }
  return(foutnet)
}

plotEGClustering = function(tissue,which.one){
  net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(net$MEs)
  # Cluster module eigengenes
  METree = flashClust::flashClust(as.dist(MEDiss), method = "average")
  MEDissThres = 0.1
  tb = table(net$moduleColors)
  names(tb) = paste0("ME",names(tb))
  tb <- tb[ METree$labels ]
  METree$labels <- paste0(names(tb), ":", tb)

  plot(METree, main = paste0("Eigengenes:",tissue," from ",which.one),
       xlab=paste0("Modules in ",tissue," from ",which.one))
}

getExchangedGenes <- function(old.partition,new.partition){
  return(old.partition[old.partition != new.partition])
}

#Arguments
#
#n.module				number of genes within the module
#n.module.and.specific	number of genes within module and with the condition (e.g. TFs, external list...)
#total.specific			number of genes in the global condition list
#total.net				number of genes in the network

testGeneSet = function(n.module,n.module.and.specific,total.specific,total.net,test="fisher",oldform=F){
  stopifnot(test == "fisher" | test == "chi")
  vector.data = matrix(ncol=2,nrow=2)

  vector.data[1,1] = n.module.and.specific
  vector.data[1,2] = n.module - n.module.and.specific
  vector.data[2,1] = total.specific - vector.data[1,1]
  if(oldform)
    vector.data[2,2] = total.net - n.module - total.specific
  else
    vector.data[2,2] = total.net - vector.data[1,2] - vector.data[2,1] - vector.data[1,1]
  #print(vector.data)
  if(test == "fisher")
    test.result = fisher.test(vector.data,alternative="greater")
  else if(test == "chi")
    test.result = chisq.test(vector.data)
  return(test.result)
}


plotModSizes = function(which.one,tissue){

  net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  tb2 <- table(net$moduleColors)[order(table(net$moduleColors))]
  barplot.title <- paste0("Module sizes for ",length(unique(net$moduleColors)),
                          " modules and ",tissue)
  barplot(tb2,col=names(tb2),
          main=barplot.title,
          ylab="Mod size",las=2)
}

generateBetaStudy <- function(expr.data,powers=c(1:30),title=NULL,plot.file=NULL,
                              net.type="signed",cor.type="pearson"){
  stopifnot(cor.type == "pearson" | cor.type == "spearman")

  if(cor.type == "pearson")
    sft = WGCNA::pickSoftThreshold(expr.data,powerVector=powers,verbose=5,moreNetworkConcepts=TRUE,
                                   networkType=net.type,corOptions=list(use='p'))
  else
    sft = WGCNA::pickSoftThreshold(expr.data,powerVector=powers,verbose=5,moreNetworkConcepts=TRUE,
                                   networkType=net.type,corFn=stats::cor,corOptions=list(method = "spearman"))

  if(!is.null(plot.file)){


    pdf(plot.file,width=18,height=8)

    old.par <- par()
    par(mfrow=c(1,2))
    cex1=0.9
    #Plotting the adjustment level

    plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab=paste0("Soft threshold, cor type ",cor.type," net type ",net.type),
         ylab="Scale Free Topology Model Fit.R^2",type="n",
         main=paste0("Scale independence ",title))

    text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,col="red",cex=cex1)

    abline(h=0.80,col="red",lwd=2)


    #Plotting mean network connectivity
    plot(sft$fitIndices[,1],sft$fitIndices[,5],xlab="Soft Threshold",
         ylab="Mean Connectivity", type="n", main=paste0("Mean connectivity (cor type ",cor.type,"net type  ",
                                                         net.type,") ",title))
    text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red",
         cex=cex1)

    par(old.par)
    dev.off()
    write.csv(sft$fitIndices,paste0(plot.file,".csv"))

  }


  sft
}

#' Title
#'
#' @param tissue
#' @param create.tom
#' @param distance.type
#' @param centroid.type
#' @param n.iterations
#' @param net.file
#' @param expr.data
#' @param beta
#' @param tom
#' @param job.path
#' @param plot.evolution
#' @param plot.go
#' @param debug
#' @param n.debug
#' @param norm.when.euc.mean
#' @param net.type
#' @param min.exchanged.genes
#' @param excludeGrey
#' @param final.net
#' @param cor.type
#'
#' @return
#' @export
#'
#' @examples
applyKMeans <- function(tissue,
                        net.file,
                        expr.data,
                        n.iterations=20,
                        debug=F,
                        n.debug=500,
                        net.type="signed",
                        min.exchanged.genes=20,
                        excludeGrey=F){

  if(typeof(expr.data) == "character")
    expr.data <- readRDS(expr.data)

  if(debug){
    cat("We are debugging, using only ",n.debug," genes")
    expr.data = expr.data[,1:n.debug]
  }

  #Step 2
  if(typeof(net.file) == "character"){
    net <- readRDS(net.file)
  }else
    net = net.file

  if(debug){
    net$moduleColors = net$moduleColors[1:n.debug]
  }


  #plot.evolution.file <- paste0(final.net,"_",distance.type,"_",centroid.type,"_evolution.kmeans.pdf")
  #partitions.file <- paste0(final.net,"_",distance.type,"_",centroid.type,"_partitions.rds")
  #evolution.file <- paste0(final.net,"_",distance.type,"_",centroid.type,"_evolution.kmeans.rds")
  #print(paste0("Evolution plot will be at ",plot.evolution.file))
  #print(paste0("Partitions data will be at ",partitions.file))
  #print(paste0("Evolution data will be at ",evolution.file))

  #ALGORITHM SPECIFICATION
  #Step 1. Let D be the expression data in which dij in D represents the expression value for
  #sample i and gene j, being s samples and g genes in total.
  #Step 2. Construct the partition by the WGCNA process, let P_D={m_1, m_2, ..., m_n} be
  #that partition where m_k is the k-th module.
  #Step 3. Get the eigengenes for each module within the partition, E={e_1, e_2, ..., e_n}
  #Step 4. Set up the k-means clustering
  #Step 4.1. Set k to n
  #Step 4.2. Set the centroids C to the eigengenes E, thus C to E
  #Step 5. Run the algorithm and monitor its evolution
  #Step 5.1 Set iterations to 0
  #Step 5.2 Create a new partition P', given C with n modules such that, for each gene, 1 <=
  #		j <= g, g_j belongs to the module c_t in C such that a distance meassure d(g_j,c_t) is
  #		minimum.
  #		Step 5.3 Calculate eigengenes of P', giving a new E'
  #		Step 5.4 Evaluate the progress. If progress done, set iterations to iterations + 1 and
  #		C to E' and go to step 5.2
  #Step 5.5 Finish
  #


  ##Step 1.

  #Gather the current partition we start from
  partition.in.colors <- net$moduleColors

  #Step 3
  eigengenes = WGCNA::moduleEigengenes(expr.data,partition.in.colors,
                                       excludeGrey=excludeGrey)

  #This variable is fixed and used as a reference to indicate the
  #modules used (they are colours but the position within the vector is
  #also relevant)
  centroid.labels <- substring(names(eigengenes$eigengenes),3)
  print("Module colors are")
  print(sort(centroid.labels))

  #Step 4
  k <- length(eigengenes$eigengenes)
  #Centroids must be a matrix with as much colums as centroids,
  #as much rows as samples
  centroids <- createCentroidMatrix(eigengenes$eigengenes)

  print("We have generated centroids")
  print(sort(centroids))


  #Step 5
  #For storage of all the partitions created during the iterations
  partitions <- list()
  #A partition will be a list of as much elements as genes and for the
  #i-th position it stores the index of the module the ith gene belongs
  #to, and the color can be found in "centroid.labels"
  new.partition <- match(partition.in.colors, centroid.labels)
  names(new.partition) <- centroid.labels[new.partition]
  partitions[[1]] <- new.partition

  #kmeans.evolution <- list()
  #kIMs <- getkIMs(tom.matrix,new.partition,length(centroid.labels))
  #kIMsWGCNA <- getkIMsFromWGCNA(expr.data,partition.in.colors,beta)
  #kMEs <- getkMEs(expr.data,new.partition,centroids)
  #kmeans.evolution[[1]] <- list(exchanged.genes=0,kMEs=kMEs,kIMs=kIMsWGCNA)

  #Launch the iterations
  #min.exchanged.genes = 20
  allgchanges = NULL
  exchanged.genes = min.exchanged.genes + 1
  iteration = 1

  while(exchanged.genes > min.exchanged.genes & iteration <= n.iterations){
    cat("k-means iteration:",iteration,"and",(n.iterations - iteration),"iterations left\n")

    new.partition <- apply(expr.data,MARGIN=2,
                           getBestModuleCor,
                           centroids=centroids,
                           signed=(net.type == "signed"),
                           cor.type="pearson")

    print(unlist(new.partition))
    print(unique(unlist(new.partition)))

    cat("We got",length(new.partition),"genes in partition\n")
    cat("We got",length(unique(new.partition)),"modules in partition\n")
    print(new.partition)
    partitions[[iteration + 1]] <- new.partition
    #Get the control values for the new partition
    exchanged.genes <- length(getExchangedGenes(partitions[[iteration]],
                                                partitions[[iteration + 1]]))
    cat(exchanged.genes,
        "genes moved to another module by k-means\n")
    new.partition.in.colors <- WGCNA::labels2colors(unlist(new.partition))
    names(new.partition.in.colors) = names(unlist(new.partition))
    print(table(new.partition.in.colors))
    cat("We got",length(unique(new.partition.in.colors)),"modules in partition\n")
    cat("We have",ncol(expr.data),"genes in expr.data\n")
    cat("We have",length(new.partition.in.colors),"genes in partition\n")

    centroids <- getNewCentroids(expr.data,new.partition.in.colors,centroid.labels)

    iteration = iteration + 1
    allgchanges = c(allgchanges,exchanged.genes)
  }
  cat("We finish with",iteration,"iterations\n")
  cat("Last number of gene changes where",exchanged.genes,"\n")
  #saveRDS(partitions,partitions.file)

  print("The algorithm finished correctly")
  net = NULL
  net$net = genNetFromPartition(expr.data.file=expr.data,
                                excludeGrey=excludeGrey,
                                partitions.file=partitions,
                                index=-1)
  net$partitions = partitions
  net$cgenes = allgchanges
  return(net)
}

getNewCentroids <- function(expr.data,partition.in.colors,centroid.labels){
  eg.vectors = WGCNA::moduleEigengenes(expr.data,
                                       partition.in.colors,excludeGrey=F)$eigengenes

  names(eg.vectors) <- substring(names(eg.vectors),3)
  eg.vectors <- eg.vectors[,centroid.labels]
  return(eg.vectors)
}


genNetFromPartition = function(expr.data.file,
                               partitions.file,
                               excludeGrey=F,
                               index=-1){

  if(typeof(partitions.file) == "character")
    parts = readRDS(partitions.file)
  else
    parts = partitions.file

  if(index < 0)
    index = length(parts)

  colors = unique(names(parts[[1]]))
  col.number = unique(parts[[1]])

  if(typeof(expr.data.file) == "character")
    expr.data = readRDS(expr.data.file)
  else
    expr.data = expr.data.file

  #if(typeof(net.file) == "character")
  #	net = readRDS(net.file)
  #else
  #	net = net.file
  new.net <- NULL
  if(index == 1){
    new.net$moduleLabels = parts[[index]]
    new.net$moduleColors = names(parts[[index]])
    names(new.net$moduleColors) = colnames(expr.data)
    names(new.net$moduleLabels) = colnames(expr.data)

  }else{
    new.net$moduleLabels = parts[[index]]
    new.net$moduleColors = colors[match(parts[[index]],col.number)]
    names(new.net$moduleColors) = colnames(expr.data)
    names(new.net$moduleLabels) = colnames(expr.data)

  }

  #If there are some grey genes as NA, add them again
  new.net$moduleColors[is.na(new.net$moduleColors)] = "grey"

  #¢if(sum(new.net$moduleColors == "grey") >= k.means.min.genes.to.consider.grey)
  #  new.net$MEs  = WGCNA::moduleEigengenes(expr.data,new.net$moduleColors,softPower=beta, excludeGrey=F)$eigengenes
  #else
  new.net$MEs  = WGCNA::moduleEigengenes(expr.data,
                                         new.net$moduleColors,
                                         excludeGrey=excludeGrey)$eigengenes
  return(new.net)
}

getBestModuleCor = function(gene,centroids,signed=TRUE,cor.type){

  return(which.max(corDistance(a=centroids,b=gene,signed=signed,cor.type=cor.type)))
}

createCentroidMatrix <- function(eigengenes){
  my.matrix <- NULL
  for(eigengene in eigengenes){
    my.matrix <- cbind(my.matrix,eigengene)
  }
  return(my.matrix)
}

#' Title
#'
#' @param colors
#' @param gcolor
#'
#' @return
#' @export
#'
#' @examples
dropGreyModule = function(colors,gcolor="grey"){
  availableColors = unique(colors)
  availableColors = availableColors[availableColors != gcolor]
  gmask = which(colors == gcolor)
  for(i in gmask)
    colors[i] = availableColors[sample(1:length(availableColors),1)]
  return(colors)

}

getAndPlotNetworkLong <- function(expr.data,beta,net.type="signed",
                                  tissue.name="SNIG",title=NULL,
                                  additional.prefix=NULL,
                                  min.cluster.size=100,
                                  save.plots=TRUE,
                                  return.tom=FALSE,
                                  excludeGrey=FALSE,
                                  max.k.cutoff = 150,
                                  r.sq.cutoff = 0.8,
                                  cor.type="pearson"){

  net <- NULL

  if(typeof(expr.data) == "character"){
    print(paste0("Reading expression from ",expr.data))
    expr.data = readRDS(expr.data)
  }

  print(paste0("We called getAndPlotNetworkLong with ",ncol(expr.data),
               " genes and ",nrow(expr.data)," samples and beta ",beta,
               " and correlation type ",cor.type,
               " and network type ", net.type," and min.cluster.size ",
               min.cluster.size, " for tissue ",tissue.name))

  print(paste0("Expression data is the following within getAndPlotNetworkLong"))
  print(expr.data[1:5,1:5])

  #We assume gene names are at columns
  gene.names <- colnames(expr.data)
  sample.names <- rownames(expr.data)
  #Lets delete unused memory
  print("Garbage collecting")
  gc()

  #If beta < 0 then we have to figure out by ourselves
  if(beta < 0){
    if(save.plots){
      b.study = generateBetaStudy(expr.data,title=paste0("Beta study for ",tissue.name),
                                  net.type=net.type,cor.type=cor.type,
                                  plot.file=paste0(additional.prefix,
                                                   "betastudy",gsub(" ","_",tissue.name),".",net.type,".pdf"))
    }else
      b.study = generateBetaStudy(expr.data,net.type=net.type,cor.type=cor.type)
    #beta = b.study$powerEstimate


    cat("Choosing beta from SFT.R.sq cut-off",r.sq.cutoff,"and max.k cut-off",max.k.cutoff,"\n")
    beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                    as.numeric(b.study$fitIndices$slope) < 0 &
                                    as.numeric(b.study$fitIndices$max.k) > max.k.cutoff,"Power"])

    if(beta == Inf){
      #OK, lets drop-off the maxk.cutoff
      beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                      as.numeric(b.study$fitIndices$slope) < 0,"Power"])
    }

    print(paste0("The estimated beta value is ",beta))
    print(paste0("The suggested beta was ",b.study$powerEstimate))

    if(beta == -Inf & is.na(b.study$powerEstimate)){
      stop("There is something wrong, our beta is",beta,"and suggested is",
           b.study$powerEstimate,"\n")
    }

    if(is.na(b.study$powerEstimate)){
      cat("We'll use",beta,"for beta\n")

    }else{
      if(beta == -Inf){
        cat("We'll use WGCNA's suggested beta\n")
        beta = b.study$powerEstimate
      }else if(beta - b.study$powerEstimate > 5){
        beta = trunc(0.5*(beta + b.study$powerEstimate))
        cat("We'll use average between WGCNA's suggested beta and ours.\n")
      }
    }

    if(beta == Inf){
      beta = 21
      cat("Warning, the final beta is",beta,"and SFT is compromised\n")
    }

    cat("The final beta value to use is:",beta,"\n")
  }

  additional.prefix = paste0(additional.prefix,tissue.name,".",beta)

  #Create the adjacency matrix of genes coming from the expression data, with the beta
  #passwd as argument
  print("Creating adjacency matrix")
  if(cor.type == "spearman")
    corOptions = "use = 'p', method = 'spearman'"
  else
    corOptions = "use = 'p'"

  adjacency = WGCNA::adjacency(expr.data, power = beta, type = net.type, corOptions = corOptions)
  print("Created")
  print(paste0("Adjacency is the following within getAndPlotNetworkLong"))
  print(adjacency[1:5,1:5])
  # Topological Overlap Matrix (TOM)
  # Turn adjacency into topological overlap
  print("Creating TOM")
  if(net.type == "signed")
    TOM = WGCNA::TOMsimilarity(adjacency)
  else if(net.type == "unsigned"){
    cat("As the network type is unsigned, the TOM type we'll create is signed")
    TOM = WGCNA::TOMsimilarity(adjacency,TOMType="signed")
  }else{
    stop(paste0("Nework type ",net.type," unrecognized when creating the network"))
  }
  #Now we can delete adjacency
  print("Deleting adjacency matrix")
  rm(adjacency)
  dissTOM = 1-TOM
  rm(TOM)
  print("Created TOM, dissTOM")
  print(dissTOM[1:5,1:5])



  # Clustering using TOM
  # Call the hierarchical clustering function that makes the clustering
  # dendrogram to grow until the leaves are genes

  geneTree = flashClust::flashClust(as.dist(dissTOM), method = "average")

  print("Now the genetree")
  print(geneTree)
  # Dynamic Tree Cut
  # We like large modules, so we set the minimum module size relatively high
  # Module identification using dynamic tree cut

  n.mods = 0
  deep.split = 2
  while(n.mods < 10 & deep.split < 5){
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                                deepSplit = deep.split,
                                                pamRespectsDendro = F,
                                                respectSmallClusters = F,
                                                minClusterSize = min.cluster.size)
    n.mods = length(table(dynamicMods))
    deep.split = deep.split + 1
  }

  dynamicColors = WGCNA::labels2colors(dynamicMods)
  print(table(dynamicColors))
  dynamicColors = CoExpNets::dropGreyModule(dynamicColors)
  print(tissue.name)
  print(table(dynamicColors))

  # Calculate eigengenes
  MEList = WGCNA::moduleEigengenes(expr.data,
                                   colors = dynamicColors,
                                   excludeGrey=excludeGrey)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs,use = "pairwise.complete.obs")
  # Cluster module eigengenes
  METree = flashClust::flashClust(as.dist(MEDiss), method = "average")

  MEDissThres = 0.1 #### MERGING THRESHOLD
  # Call an automatic merging function
  merge = WGCNA::mergeCloseModules(expr.data, dynamicColors, cutHeight = MEDissThres,
                                   verbose = 3, unassdColor="grey",getNewUnassdME = FALSE)
  # The merged module colors
  mergedColors = CoExpNets::dropGreyModule(merge$colors)
  # Eigengenes of the new merged modules
  mergedMEs = merge$newMEs
  print(table(mergedColors))

  if(save.plots){
    dendro.name = paste0(additional.prefix,"_dendro_colors.pdf")
    eigengenes.name = paste0(additional.prefix,"_Eigengenes_clustering.pdf")

    pdf(file=dendro.name)
    WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                               c("Dynamic Tree Cut", "Merged dynamic"),
                               dendroLabels = FALSE,
                               hang = 0.03, addGuide = TRUE,
                               guideHang = 0.05,main=title)
    dev.off()

    pdf(file=eigengenes.name)
    # Plot the result
    tb <- tb[ METree$labels ]
    METree$labels <- paste0(names(tb), ":", tb)

    plot(METree, main = paste0(title," ",tissue.name,
                               " Clustering of module eigengenes"),
         xlab = "", sub = "")
    # Plot the cut line into the dendrogram
    abline(h=MEDissThres, col = "red")
    dev.off()
  }
  #Creating the connectivity analysis
  #expr.data.sets <- NULL
  #expr.data.sets[[1]] <- expr.data
  #expr.data.sets.names <- vector(mode="character",length=1)
  #expr.data.sets.names[1] <- paste0("Tissue ",tissue.name)
  #pdf(file=paste0(additional.prefix,"_connectivity.pdf"))
  #connectivity(expr.data.sets,expr.data.sets.names,bethas=c(betha))
  #def.off()


  #Prepare for creating the network objecto to return
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", WGCNA::standardColors(400) )
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  tb2 <- table(moduleColors)[order(table(moduleColors))]
  print(tb2)

  if(save.plots){
    pdf(file=paste0(additional.prefix,"_mod_size.pdf"))
    barplot.title <- paste0("Mod. sizes for ",length(unique(moduleColors)),
                            " modules and ",title)
    barplot(tb2,col=names(tb2),
            main=barplot.title,
            xlab="Mod colors",ylab="Mod size")
    dev.off()
  }

  net$MEs <- MEs
  rownames(net$MEs) <- sample.names
  net$moduleLabels <- moduleLabels
  net$moduleColors <- moduleColors
  net$beta = beta
  net$type = net.type
  names(net$moduleColors) <- gene.names
  net$geneTree <- geneTree
  if(return.tom){
    return(list(net=net,tom=(1 - dissTOM)))
  }else
    return(net)
}

corDistance = function(a,b,signed=TRUE,cor.type="pearson"){
  if(cor.type=="pearson"){
    if(signed)
      #return(0.5 + 0.5*WGCNA::corFast(x=a,y=b)) #(Note they are equivalent)
      return(0.5 * (1 + stats::cor(a,b)))
    return(abs(stats::cor(a,b)))
  }else{
    if(signed)
      #return(0.5 + 0.5*WGCNA::corFast(a,b)) #(Note they are equivalent)
      return(0.5 * (1 + stats::cor(a,b,method=cor.type)))
    return(abs(stats::cor(a,b,method=cor.type)))
  }
}


getMM = function(net=NULL,
                 expr.data.file=NULL,
                 tissue,genes,
                 table.format=FALSE,
                 which.one="rnaseq",
                 silent=F,keep.grey=F,
                 alt.gene.index=NULL){

  net = getNetworkFromTissue(tissue,which.one)
  if(is.null(expr.data.file))
    expr.data = getExprDataFromTissue(tissue,which.one)
  else
    expr.data = readRDS(expr.data.file)

  if(is.null(expr.data)){
    cat("There is no expr.data file registered for category", which.one,"and tissue",tissue,"\n")
    return(expr.data)
  }


  if(which.one == "micro19K"){
    names(net$moduleColors) = colnames(expr.data)
  }


  mm = NULL

  if(is.null(genes)){
    #There is no correlation within grey module normally
    genes = names(net$moduleColors) #[net$moduleColors != "grey"]
  }

  if(table.format){
    out.table = data.frame(list(ensgene=character(0),name=character(0),
                                module=character(0),mm=numeric(0)),stringsAsFactors=FALSE)
    n.row = 1
    out.table[1:length(genes),1] = genes
    out.table[1:length(genes),2] = fromAny2GeneName(genes)
    out.table[1:length(genes),3] = net$moduleColors[match(genes,names(net$moduleColors))]
  }

  for(gene in genes){
    if(!is.null(alt.gene.index)){
      module = net$moduleColors[alt.gene.index[which(genes %in% gene)]]
      expr.data.gene.index = alt.gene.index[which(genes %in% gene)]
    }else{
      module = net$moduleColors[names(net$moduleColors) %in% gene]
      expr.data.gene.index = gene
    }
    if(length(module) == 0){
      if(!silent)
        cat("Gene ",gene," not in network\n")
      mm[[gene]] = -1
    }else{
      if(length(module) > 1)
        module = module[1]
      if(module == "grey" & !keep.grey)
        mm[[gene]] = 0
      else{
        if(length(module) > 1){
          if(!silent)
            print(paste0("Gene ",gene, " found in modules ",module))
          module = module[1]
          if(!silent)
            print(paste0("Keeping only ",module))
        }

        tryCatch(mm[[gene]] <- cor(net$MEs[paste0("ME",module)],expr.data[,expr.data.gene.index]),
                 error = function(e){
                   print(paste0("Error in module ",module," ",e))
                 })
      }
    }
  }
  if(table.format){
    out.table[1:length(genes),4] = unlist(mm[genes])
    return(out.table)
  }
  return(mm)

}


corWithCatTraits = function(tissue,which.one,covlist,covs=NULL){
  MEs = getNetworkEigengenes(tissue=tissue,which.one=which.one)
  if(is.null(covs))
    covs = getCovariates(tissue=tissue,which.one=which.one)
  covs = covs[,covlist]

  for(i in 1:ncol(covs)){
    if(typeof(covs[,i]) ==  "character")
      covs[,i] = as.factor(covs[,i])
  }
  factor.mask = NULL
  for(i in 1:ncol(covs)){
    factor.mask = c(factor.mask,is.factor(covs[,i]))
  }

  fcm = matrix(nrow=ncol(MEs),ncol=sum(factor.mask))
  index = 1
  for(i in which(factor.mask)){
    #cat("Factor",colnames(trait.data)[i],"\n")
    #cat("Levels",levels(trait.data[,i]))
    #print(trait.data[,i])
    for(j in 1:ncol(MEs)){
      #print(paste0(i,j)
      if(length(unique(covs[,i])) > 1){
        form = eg ~ cov
        data.in = data.frame(MEs[,j],covs[,i])
        colnames(data.in) = c("eg","cov")
        fcm[j,index] = anova(aov(form,data.in))$`Pr(>F)`[1]
      }else
        fcm[j,index] = 1

    }
    fcm[,index] = p.adjust(fcm[,index],method="BH")
    index = index + 1
  }

  moduleTraitCor = cor(MEs,covs[,!factor.mask],use="p")

  #Generate the p-values for significance of a given matrix of correlations, for all modules,
  #between traits data and eigengenes, both from samples
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nrow(MEs))
  moduleTraitPvalue = cbind(moduleTraitPvalue,fcm)
  colnames(moduleTraitPvalue) = c(colnames(covs)[!factor.mask],
                                  colnames(covs)[factor.mask])
  moduleTraitPvalue = -log10(moduleTraitPvalue)
  moduleTraitPvalue[moduleTraitPvalue > 10] = 10
  WGCNA::labeledHeatmap(Matrix=moduleTraitPvalue,
                        xLabels=colnames(moduleTraitPvalue),
                        yLabels=gsub("ME","",names(MEs)),
                        ySymbols=names(MEs),
                        colorLabels=FALSE,
                        colors=rev(heat.colors(50)),
                        cex.text=0.5,
                        zlim = c(0,10),
                        main="Module-trait relationships")
}


corWithNumTraits = function(tissue,which.one,covlist,covs=NULL){
  MEs = getNetworkEigengenes(tissue=tissue,which.one=which.one)
  if(is.null(covs))
    covs = getCovariates(tissue=tissue,which.one=which.one)
  covs = covs[,covlist]
  moduleTraitCor = cor(MEs, covs, use = "p")
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(MEs))
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                        xLabels = covlist,
                        yLabels = names(MEs),
                        ySymbols = names(MEs),
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = 0.5,
                        zlim = c(-1,1),
                        main = paste0("Module-trait relationships"))
}


trasposeDataFrame = function(file.in,first.c.is.name=F){
  if(typeof(file.in) == "character")
    data.in = readRDS(file.in)
  else{
    data.in = file.in
    rm(file.in)
  }

  if(first.c.is.name){
    data.t = as.data.frame(cbind(apply(data.in[,-1],MARGIN=1,function(x){ return(as.numeric(x))})))
    colnames(data.t) = data.in[,1]
    rownames(data.t) = colnames(data.in[-1])
  }else{
    data.t = as.data.frame(cbind(apply(data.in,MARGIN=1,function(x){ return(as.numeric(x))})))
    colnames(data.t) = rownames(data.in)
    rownames(data.t) = colnames(data.in)
  }
  return(data.t)
}

saveTOM = function(tom, clusters, filepref){
  modules = unique(clusters)
  metadata = NULL
  lapply(modules,function(x){
    mask = clusters %in% x
    tomname = paste0(filepref,"_",x,".rds")
    saveRDS(tom[mask,mask],tomname)
    metadata[[x]] <<- list(mask=mask,tomname=tomname)
  })
  saveRDS(metadata,paste0(filepref,"_metadata.rds"))
}

readTOM = function(filepref){
  mdfile = paste0(filepref,"_metadata.rds")
  if(!file.exists(mdfile))
    stop(paste0("Error: metadata file",mdfile," not found"))

  metadata = readRDS(mdfile)
  modules = names(metadata)
  size = length(metadata[[modules[1]]]$mask)
  tom = matrix(ncol=size,nrow=size)
  tom[] = 0
  for(module in modules){
    cat("Reading module",module,"\n")
    smalltom = readRDS(metadata[[module]]$tomname)
    mask = metadata[[module]]$mask
    tom[mask,mask] = smalltom
  }
  return(tom)
}

removeTOM = function(filepref){
  mdfile = paste0(filepref,"_metadata.rds")
  if(!file.exists(mdfile))
    stop(paste0("Error: metadata file",mdfile," not found"))

  metadata = readRDS(mdfile)
  modules = names(metadata)
  for(module in modules)
    file.remove(metadata[[module]]$tomname)
  file.remove(paste0(filepref,"_metadata.rds"))
}

