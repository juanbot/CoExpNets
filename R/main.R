

#' Title
#'
#' @param net
#' @param folder
#'
#' @return
#' @export
#'
#' @examples
milkNet = function(net,folder){
  netName = NULL
  if(typeof(net) == "character"){
    netName = gsub(".rds","",basename(net))
    stopifnot(file.exists(net))
    net = readRDS(net)
  }

  if(is.null(netName))
    netName = gsub(".rds","",basename(net$file))
  #Print the clusters
  clusters = NULL
  if(is.null(names(net$moduleColors)))
    names(net$moduleColors) = names(net$adjacency)
  clusters = list(genes=names(net$moduleColors),modules=net$moduleColors)
  #              stringsAsFactors=F)
  dir.create(folder)
  write.table(clusters,paste0(folder,netName,"_clusters.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  if(!is.null(net$adjacency)){
    adj = list(genes=names(net$moduleColors),adjacency=net$adjacency)
    write.table(adj,paste0(folder,netName,"_adjacency.tsv"),quote=F,row.names=F,col.names=T,sep="\t")

  }

  if(!is.null(net$go))
    write.table(net$go,paste0(folder,netName,"_funcAnnotation.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  else
    cat("No functional annotation to generate for",netName,"\n")
  if(!is.null(net$ct)){
    ct = net$ct
    ct = cbind(rownames(net$ct),ct)
    colnames(ct)[1] = "markerset"
    write.table(ct,paste0(folder,netName,"_cellTypes.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  }

  else
    cat("No cell annotation to generate for",netName,"\n")

  basicFacts = c(netName,net$mode,net$beta,
                 length(net$moduleColors),
                 length(unique(net$moduleColors)),
                 net$file)
  print(basicFacts)

  write.table(list(info=basicFacts,
                   value=c("name","mode","beta","ngenes","nmodules","netfile")),
              paste0(folder,netName,"_basicInfo.tsv"),
              quote=F,row.names=F,col.names=T,sep="\t")

  if(!is.null(net$subnets)){
    rdiffs = NULL
    fdiffs = NULL
    adiffs = NULL
    clusters = NULL

    indexes = names(net$subnets)
    if(is.null(indexes))
      indexes = 1:length(net$subnets)
    for(i in indexes){
      #cat(i)
      if(!is.null(net$subnets[[i]]$subcluster)){
        #cat("Other here\n")
        clusters = rbind(clusters,net$subnets[[i]]$subcluster)
      }else{
        #cat("here",i,"\n")
        #str(net$subnets[[i]])
        cluster = unlist(net$subnets[[i]]$partitions[length(net$subnets[[i]]$partitions)])
        #print(str(cluster))
        discGenes = net$subnets[[i]]$discGenes
        if(length(discGenes) > 0){
          newcluster = NULL
          for(gene in 1:length(genes)){
            if(genes[gene] %in% discGenes)
              newcluster = c(newcluster,"null")
            else{
              module = cluster[names(cluster) == genes[gene]]
              newcluster = c(newcluster,module)
            }

          }
          cat("#")
          #newcluster = cluster[match(names(cluster),genes)]
          ##newcluster = vector(mode="character",length = length(genes))
          #cluster[match()]
          #print(discGenes)
          #print(length(genes))
          #maskgenes = genes[!(genes %in% discGenes)]
          #maskgenes = match(maskgenes,genes)
          #print(length(maskgenes))
          #print(length(cluster))
          #newcluster = vector(mode="character",length = length(genes))
          #newcluster[mask] = cluster
          #newcluster[!mask] = "null"
        }else
          cat("-")
        clusters = rbind(clusters,cluster)
      }


    }

    if(!is.null(net$allsamplesnet)){
      if(is.null(net$allsamplesnet$moduleColors)){
        net$allsamplesnet$moduleColors = net$allsamplesnet$partitions[[length(net$allsamplesnet$partitions)]]
      }
    }
    for(i in 2:nrow(clusters)){
      rdiffs = c(rdiffs,mclust::adjustedRandIndex(clusters[i,],clusters[i-1,]))
      fdiffs = c(fdiffs,mclust::adjustedRandIndex(clusters[i-1,],net$moduleColors))
      if(!is.null(net$allsamplesnet)){
        adiffs = c(adiffs,mclust::adjustedRandIndex(clusters[i,],net$allsamplesnet$moduleColors))

      }
    }
    if(is.null(adiffs))
      adiffs = rep(0,length(rdiffs))
    print(rdiffs)
    print(adiffs)
    print(fdiffs)

    bootStats = list(rdiffs=rdiffs,adiffs=adiffs,fdiffs=fdiffs)
    write.table(bootStats,paste0(folder,netName,"_bootStrapStats.tsv"),
                quote=F,row.names=F,col.names=T,sep="\t")

  }

}

print.bootnet = function(net){
  cat("A bootstrapped network created with mode",net$mode,"\n",
      "Soft thresholding parameter (beta):",net$beta,"\n",
      "Adjacency summary\n")
  genes = names(net$moduleColors)
  if(length(genes) == 0)
    genes = names(net$adjacency)

  print(summary(net$adjacency))
  cat("Final number of modules",length(unique(net$moduleColors)),"\n")
  if(!is.null(net$allsamplesnet))
    cat("All samples net final number of modules",
        length(unique(net$allsamplesnet$moduleColors)),"\n")
  else
    cat("No all samples net available\n")

  rdiffs = NULL
  fdiffs = NULL
  adiffs = NULL
  clusters = NULL

  indexes = names(net$subnets)
  if(is.null(indexes))
    indexes = 1:length(net$subnets)
  for(i in indexes){
    #cat(i)
    if(!is.null(net$subnets[[i]]$subcluster)){
      #cat("Other here\n")
      clusters = rbind(clusters,net$subnets[[i]]$subcluster)
    }else{
      #cat("here",i,"\n")
      #str(net$subnets[[i]])
      cluster = unlist(net$subnets[[i]]$partitions[length(net$subnets[[i]]$partitions)])
      #print(str(cluster))
      discGenes = net$subnets[[i]]$discGenes
      if(length(discGenes) > 0){
        newcluster = NULL
        for(gene in 1:length(genes)){
          if(genes[gene] %in% discGenes)
            newcluster = c(newcluster,"null")
          else{
            module = cluster[names(cluster) == genes[gene]]
            newcluster = c(newcluster,module)
          }

        }
        cat("#")
        #newcluster = cluster[match(names(cluster),genes)]
        ##newcluster = vector(mode="character",length = length(genes))
        #cluster[match()]
        #print(discGenes)
        #print(length(genes))
        #maskgenes = genes[!(genes %in% discGenes)]
        #maskgenes = match(maskgenes,genes)
        #print(length(maskgenes))
        #print(length(cluster))
        #newcluster = vector(mode="character",length = length(genes))
        #newcluster[mask] = cluster
        #newcluster[!mask] = "null"
      }else
        cat("-")
      clusters = rbind(clusters,cluster)
    }


  }

  for(i in 2:nrow(clusters)){
    rdiffs = c(rdiffs,mclust::adjustedRandIndex(clusters[i,],clusters[i-1,]))
    fdiffs = c(fdiffs,mclust::adjustedRandIndex(clusters[i-1,],net$moduleColors))
    if(!is.null(net$allsamplesnet))
      adiffs = c(adiffs,mclust::adjustedRandIndex(clusters[i,],net$allsamplesnet$moduleColors))

  }
  cat("Successive Rand simmilarities\n")
  print(rdiffs)
  cat("Rand differences with all samples net\n")
  print(adiffs)
  cat("Rand differences with final net\n")
  print(fdiffs)
  return(list(rdiffs=rdiffs,adiffs=adiffs,fdiffs=fdiffs))
}

plot.bootnet = function(net){
  cat("A bootstrapped network created with mode",net$mode,"\n",
      "Soft thresholding parameter (beta):",net$beta,"\n")
  plotModSizes(which.one = "new",tissue=net)
  rdiffs = NULL
  fdiffs = NULL
  adiffs = NULL
  clusters = NULL

  indexes = names(net$subnets)
  if(is.null(indexes))
    indexes = 1:length(net$subnets)
  for(i in indexes){
    if(!is.null(net$subnets[[i]]$subcluster)){
      clusters = rbind(clusters,net$subnets[[i]]$subcluster)
    }else{
      cluster = unlist(net$subnets[[i]]$partitions[length(net$subnets[[i]]$partitions)])
      discGenes = net$subnets[[i]]$discGenes
      if(length(discGenes) > 0){
        newcluster = NULL
        for(gene in 1:length(genes)){
          if(genes[gene] %in% discGenes)
            newcluster = c(newcluster,"null")
          else{
            module = cluster[names(cluster) == genes[gene]]
            newcluster = c(newcluster,module)
          }
        }
      }
      clusters = rbind(clusters,cluster)
    }
  }

  for(i in 2:nrow(clusters)){
    rdiffs = c(rdiffs,mclust::adjustedRandIndex(clusters[i,],clusters[i-1,]))
    fdiffs = c(fdiffs,mclust::adjustedRandIndex(clusters[i-1,],net$moduleColors))
    if(!is.null(net$allsamplesnet))
      adiffs = c(adiffs,mclust::adjustedRandIndex(clusters[i,],net$allsamplesnet$moduleColors))

  }
  print(rdiffs)

  plot(rdiffs,type="lp",col="black",main="Bootstrap evolution",
       xlab="Number of subnets",ylab="Rand simmilarities")
  if(!is.null(adiffs))
    lines(adiffs,col="blue")
  lines(fdiffs,col="red")
  legend("topleft",
         legend = c("Succesive Rand simmilarity",
                    "Simmilarity with final net",
                    "Simmilarity with cannonical net"),
         fill=c("black","red","blue"))

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
    limma::plotMDS(rpkms.net[,mask], #col=colors[as.numeric(covs[,covvar])],
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
#' @param removeTOMF
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
                                 removeTOMF=F,
                                 job.path,
                                 min.cluster.size=100,
                                 allsampsnet=F,
                                 each=1,
                                 tissue="Bootstrap",
                                 blockTOM=F,
                                 #clParams=" -l nodes=1:nv ",
                                 clParams="-l h_rt=72:0:0 -l tmem=16.9G,h_vmem=16.9G",
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
      #datain = log2(1 + as.matrix(datain))
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
                                 fullAnnotation = F,
                                 expr.data=datain,
                                 job.path=outfolder)
      return(readRDS(net))

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
  params$removeTOMF = removeTOMF
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
                       removeTOMF=F,
                       each=100,
                       mode,
                       waitFor=24*3600,
                       indexes,
                       job.path){

  time = Sys.time()
  ngenes = ncol(expr.data)
  genes = colnames(expr.data)

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
            #localTOM = CoExpNets::readTOM(net$tom)
            mdfile = paste0(net$tom,"_metadata.rds")
            if(!file.exists(mdfile))
              stop(paste0("Error: metadata file",mdfile," not found"))

            metadata = readRDS(mdfile)
            modules = names(metadata)
            for(module in modules){
              cat("Reading module",module,"\n")
              smalltom = readRDS(metadata[[module]]$tomname)
              mask = metadata[[module]]$mask
              TOM[mask,mask] = TOM[mask,mask] + smalltom
            }

            CoExpNets::removeTOM(net$tom)
          }else{
            localTOM = readRDS(net$tom)
            if(removeTOMF)
              file.remove(net$tom)
            cat("Quantile normalization\n")
            localTOM = preprocessCore::normalize.quantiles(localTOM)
            print("Done")
            if(length(net$discGenes) > 0){
              mask = !(genes %in% net$discGenes)
              TOM[mask,mask] = TOM[mask,mask] + localTOM
            }else
              TOM = TOM + localTOM
            rm("localTOM")
            print("Done")
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

            validgenes = unlist(apply(expr.data,2,function(x){
              return(var(x) != 0)
            }))

            outnet = CoExpNets::applyKMeans(tissue=tissue,
                                            n.iterations=n.iterations,
                                            net.file=localnet,
                                            expr.data=expr.data)

            net$subcluster = outnet$moduleColors

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
    finalnet$file = paste0(job.path,"/netBoot",tissue,".",
                           finalnet$beta,".it.",n.iterations,
                           ".b.",nrow(indexes),".rds")

    TOM = TOM/tomCount
    adjacency = apply(TOM,2,sum)
    adjacency = adjacency/length(adjacency)
    names(adjacency) = colnames(expr.data)
    if(blockTOM){
      finalnet$tom = paste0(finalnet$file,".tom.")
      CoExpNets::saveTOM(TOM,outnet$moduleColors,finalnet$tom)
    }else{
      finalnet$tom = paste0(finalnet$file,".tom.rds")
      saveRDS(TOM,finalnet$tom)
    }


    finalnet$adjacency = adjacency
    finalnet$moduleColors = outnet$moduleColors
    finalnet$subnets = allsubnets
    finalnet$mode = mode
    rm("TOM")
    outnet = CoExpNets::applyKMeans(tissue=tissue,
                                    n.iterations=n.iterations,
                                    net.file=finalnet,
                                    expr.data=expr.data)
    finalnet$moduleColors = outnet$moduleColors
    finalnet$MEs = outnet$MEs
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
                                                             save.tom=F)
  }

  if(!is.null(finalnet)){
    attr(finalnet,"class") <- "bootnet"
    saveRDS(finalnet,finalnet$file)
    return(finalnet$file)
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
                               annotateFinalNet=F,
                               each=1,
                               blockTOM=F,
                               removeTOMF=F,
                               min.cluster.size=100,
                               tissue="Bootstrap",
                               b=10,...){
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
  genes = colnames(expr.data)
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
    net = readRDS(net)
    cat("Accumulating TOM",net$tom,"\n")
    if(blockTOM){
      cat("Reading TOM\n")
      #localTOM = CoExpNets::readTOM(net$tom)
      localTOM = readTOM(net$tom)
      CoExpNets::removeTOM(net$tom)
    }else{
      localTOM = readRDS(net$tom)
      if(removeTOMF)
        file.remove(net$tom)
      cat("Quantile normalization\n")
      localTOM = preprocessCore::normalize.quantiles(localTOM)
      print("Done")

    }
    if(length(net$discGenes) > 0){
      print(net$discGenes)
      print(genes)
      mask = !(genes %in% net$discGenes)
      #print(sum(mask))
      #print(str(TOM))
      #print(str(localTOM))
      TOM[mask,mask] = TOM[mask,mask] + localTOM
    }else{
      cat("No discarged genes\n")
      #print(str(TOM))
      #print(str(localTOM))
      TOM = TOM + localTOM
    }


    rm("localTOM")
    print("Done")

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

      net$subcluster = outnet$moduleColors

    }else
      cat("No entering here\n")
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
  finalnet$moduleColors = outnet$moduleColors
  finalnet$subnets = allsubnets

  outnet = applyKMeans(tissue=tissue,
                       n.iterations=n.iterations,
                       net.file=finalnet,
                       silent=F,
                       expr.data=expr.data)

  finalnet$moduleColors = outnet$moduleColors
  names(finalnet$moduleColors) = names(finalnet$adjacency)
  finalnet$MEs = outnet$MEs
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
  if(annotateFinalNet){
    go = CoExpNets::getGProfilerOnNet(net.file=finalnet$file,
                                      out.file=paste0(finalnet$file,"_gprof.csv"))

    write.csv(CoExpNets::genAnnotationCellType(net.in=finalnet$file,
                                               return.processed = F),
              paste0(finalnet$file,"_celltype.csv"))
  }

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
                                excludeGrey=FALSE,
                                fullAnnotation=T,
                                silent=T){

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

  validgenes = unlist(apply(expr.data,2,function(x){
    return(var(x) != 0)
  }))
  if(sum(validgenes) < ncol(expr.data))
    cat("There are genes with 0 variance:",
        paste0(colnames(expr.data)[!validgenes],collapse=","),"\n")
  discGenes = colnames(expr.data)[!validgenes]
  expr.data = expr.data[,validgenes]

  net.and.tom = getAndPlotNetworkLong(expr.data=expr.data,
                                      beta=beta,
                                      tissue.name=tissue,
                                      min.cluster.size=min.cluster.size,
                                      save.plots=save.plots,
                                      excludeGrey=excludeGrey,
                                      additional.prefix=job.path,
                                      return.tom=T,
                                      cor.type=cor.type,
                                      silent=silent)

  if(is.null(final.net))
    final.net = paste0(job.path,"/","net",tissue,".",
                       net.and.tom$net$beta,".it.",n.iterations,".rds")
  outnet = applyKMeans(tissue=tissue,
                       n.iterations=n.iterations,
                       net.file=net.and.tom$net,
                       expr.data=expr.data,
                       excludeGrey=excludeGrey,
                       min.exchanged.genes = min.exchanged.genes,
                       silent=silent)


  if(save.tom){
    if(blockTOM)
      saveTOM(tom=net.and.tom$tom,
              clusters=outnet$moduleColors,
              filepref=paste0(final.net,".tom."))
    else
      saveRDS(net.and.tom$tom,paste0(final.net,".tom.rds"))
  }

  outnet$beta = net.and.tom$net$beta
  outnet$file = final.net
  outnet$adjacency = net.and.tom$net$adjacency
  names(outnet$moduleColors) = colnames(expr.data)
  outnet$discGenes = discGenes

  if(save.tom){
    if(blockTOM)
      outnet$tom = paste0(final.net,".tom.")
    else
      outnet$tom = paste0(final.net,".tom.rds")
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
  saveRDS(outnet,final.net)
  if(fullAnnotation){
    go = CoExpNets::getGProfilerOnNet(net.file=final.net,
                                      out.file=paste0(final.net,"_gprof.csv"))

    write.csv(CoExpNets::genAnnotationCellType(net.in=final.net,
                                               return.processed = F),
              paste0(final.net,"_celltype.csv"))
  }
  return(final.net)
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
getSCDownstreamNetwork = function(tissue="mytissue",
                                  n.iterations=50,	#Number of iterations for k-means, 50 recommended
                                  min.exchanged.genes=20,
                                  alg="TOM",
                                  beta=-1,
                                  expr.data, 	#We expect a file name pointing to a dataframe (RDS format) with
                                  #genes in columns and samples in rows. Each gene name appears
                                  #in the column name. Better to use gene symbols as names
                                  job.path="~/tmp/",	#Where to store all results
                                  min.cluster.size=30,		#Minimum number of genes to form a cluster
                                  net.type="signed",			#Leave it like that (see WGCNA docs)
                                  debug=F,
                                  save.plots=F,
                                  excludeGrey=FALSE,
                                  fullAnnotation=T,
                                  silent=T){

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

  validgenes = unlist(apply(expr.data,2,function(x){
    return(var(x) != 0)
  }))
  if(sum(validgenes) < ncol(expr.data))
    cat("There are genes with 0 variance:",
        paste0(colnames(expr.data)[!validgenes],collapse=","),"\n")
  discGenes = colnames(expr.data)[!validgenes]
  expr.data = expr.data[,validgenes]

  if(alg == "singlecell")
    net = getAndPlotSCNetworkLong(expr.data=expr.data,
                                  tissue.name=tissue,
                                  min.cluster.size=min.cluster.size,
                                  save.plots=save.plots,
                                  excludeGrey=excludeGrey,
                                  additional.prefix=job.path,
                                  cor.type=cor.type,
                                  silent=silent)
  else if(alg == "notom")
    net = getAndPlotSCNoTOMNetworkLong(expr.data=expr.data,
                                       tissue.name=tissue,
                                       beta=beta,
                                       min.cluster.size=min.cluster.size,
                                       save.plots=save.plots,
                                       excludeGrey=excludeGrey,
                                       additional.prefix=job.path,
                                       cor.type=cor.type,
                                       silent=silent)
  else
    stop("Algorithm type for network creation not recognised")

  if(is.null(net))
    return(net)
  if(is.null(final.net))
    final.net = paste0(job.path,"/","net",tissue,".",
                       net$beta,".it.",n.iterations,".rds")

  outnet = applyKMeans(tissue=tissue,
                       n.iterations=n.iterations,
                       net.file=net,
                       expr.data=expr.data,
                       excludeGrey=excludeGrey,
                       min.exchanged.genes = min.exchanged.genes,
                       silent=silent)



  outnet$beta = NA
  outnet$file = final.net
  outnet$adjacency = NA
  names(outnet$moduleColors) = colnames(expr.data)
  outnet$discGenes = discGenes
  saveRDS(outnet,final.net)

  if(save.plots){
    cat("Generating mod sizes for",final.net,"\n")
    pdf(paste0(final.net,".mod_size.pdf"))
    plotModSizes(which.one="new",tissue=final.net)
    dev.off()
    pdf(paste0(final.net,".Eigengenes_clustering.pdf"))
    plotEGClustering(which.one="new",tissue=final.net)
    dev.off()
  }
  if(fullAnnotation){
    go = CoExpNets::getGProfilerOnNet(net.file=final.net,
                                      out.file=paste0(final.net,"_gprof.csv"))

    write.csv(CoExpNets::genAnnotationCellType(net.in=final.net,
                                               return.processed = F),
              paste0(final.net,"_celltype.csv"))
  }
  return(final.net)
}

plotEGClustering = function(tissue,which.one){
  net = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one=which.one)
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

#' Title
#'
#' @param n.module
#' @param n.module.and.specific
#' @param total.specific
#' @param total.net
#' @param test
#' @param oldform
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param which.one
#' @param tissue
#'
#' @return
#' @export
#'
#' @examples
plotModSizes = function(which.one,tissue){

  net = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one=which.one)
  tb2 <- table(net$moduleColors)[order(table(net$moduleColors))]
  if(typeof(tissue) == "character")
    barplot.title <- paste0("Module sizes for ",length(unique(net$moduleColors)),
                            " modules and ",tissue)
  else
    barplot.title <- paste0("Module sizes for ",length(unique(net$moduleColors)),
                            " modules")
  barplot(tb2,col=names(tb2),
          main=barplot.title,
          ylab="Mod size",las=2)
}

#' Title
#'
#' @param expr.data
#' @param powers
#' @param title
#' @param plot.file
#' @param net.type
#' @param cor.type
#'
#' @return
#' @export
#'
#' @examples
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

#' Title K-means refinement process after getting TOM and clusters from `WGCNA`
#' This method applies `n.iterations` k-means iterations to the hierarchical clustering
#' generated partition from `WGCNA`. It uses the eigengenes as centroids and the distance
#' between gene pairs is calculated on the basis of whether the network is signed or not.
#' For details on the approach see paper
#' <https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6>
#' #Step 1. Let D be the expression data in which dij in D represents the expression value for
#' sample i and gene j, being s samples and g genes in total.
#' Step 2. Construct the partition by the WGCNA process, let P_D={m_1, m_2, ..., m_n} be
#' that partition where m_k is the k-th module.
#' Step 3. Get the eigengenes for each module within the partition, E={e_1, e_2, ..., e_n}
#' Step 4. Set up the k-means clustering
#' Step 4.1. Set k to n
#' Step 4.2. Set the centroids C to the eigengenes E, thus C to E
#' Step 5. Run the algorithm and monitor its evolution
#' Step 5.1 Set iterations to 0
#' Step 5.2 Create a new partition P', given C with n modules such that, for each gene, 1 <=
#' j <= g, g_j belongs to the module c_t in C such that a distance meassure d(g_j,c_t) is
#' minimum.
#' Step 5.3 Calculate eigengenes of P', giving a new E'
#' Step 5.4 Evaluate the progress. If progress done, set iterations to iterations + 1 and
#' C to E' and go to step 5.2
#' Step 5.5 Finish
#
#'
#' @param tissue A tissue name
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
                        excludeGrey=F,
                        silent=T){

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

  ##Step 1.

  #Gather the current partition we start from
  #partition.in.colors <- net$moduleColors
  clusters = net$moduleColors
  geneNames = names(clusters)
  #Step 3
  eigengenes = WGCNA::moduleEigengenes(expr.data,net$moduleColors,
                                       excludeGrey=excludeGrey)

  #This variable is fixed and used as a reference to indicate the
  #modules used (they are colours but the position within the vector is
  #also relevant)
  #centroid.labels <- substring(names(eigengenes$eigengenes),3)
  if(!silent){
    print("Module colors are")
    print(sort(unique(net$moduleColors)))
  }

  #Step 4
  k = length(eigengenes$eigengenes)
  if(!silent)
    cat("Working with",k,"modules/centroids\n")
  #Centroids must be a matrix with as much colums as centroids,
  #as much rows as samples
  centroids = createCentroidMatrix(eigengenes$eigengenes)

  #print("We have generated centroids")
  #print(sort(centroids))

  #Step 5
  #For storage of all the partitions created during the iterations
  partitions = list()
  #A partition will be a list of as much elements as genes and for the
  #i-th position it stores the index of the module the ith gene belongs
  #to, and the color can be found in "centroid.labels"
  #new.partition <- match(partition.in.colors, centroid.labels)
  #names(new.partition) <- centroid.labels[new.partition]
  partitions[[1]] = clusters

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

  new.clusters = NULL
  while(exchanged.genes > min.exchanged.genes & iteration <= n.iterations){
    if(!silent)
      cat("k-means iteration:",iteration,"and",(n.iterations - iteration),"iterations left\n")
    #print(corDistance)
    #print(centroids)
    #print(sum(is.na(centroids)))
    #print(sum(is.na(expr.data)))
    new.clusters = NULL
    for(i in 1:ncol(expr.data)){
      newc = getBestModuleCor(gene=expr.data[,i],centroids,signed=(net.type == "signed"),
                              cor.type="pearson")
      if(length(newc) == 0 | is.null(newc)){
        if(!silent)
          cat("Gene",colnames(expr.data)[i],"\n")
        if(!silent)
          print(expr.data[,i])
        if(!silent)
          print(newc)
        stop("Something wrong in the data")
      }
      new.clusters = c(new.clusters,newc)

    }
    #new.clusters <- unlist(apply(expr.data,2,
    #                             getBestModuleCor,
    #                             centroids=centroids,
    #                             signed=(net.type == "signed"),
    #                             cor.type="pearson"))
    if(!silent)
      print(sum(is.na(expr.data)))

    if(!silent)
      print(table(new.clusters))
    new.clusters = colnames(centroids)[new.clusters]
    names(new.clusters) = geneNames
    #print(table(new.clusters))

    if(!silent)
      cat("We got",length(new.clusters),"genes in partition\n")
    if(!silent)
      cat("We got",length(unique(new.clusters)),"modules in partition\n")
    #print(unique(new.clusters))
    partitions[[iteration + 1]] <- new.clusters
    #Get the control values for the new partition
    exchanged.genes <- length(getExchangedGenes(partitions[[iteration]],
                                                partitions[[iteration + 1]]))
    if(!silent)
      cat(exchanged.genes,
          "genes moved to another module by k-means\n")
    if(!silent)
      cat("We have",ncol(expr.data),"genes in expr.data\n")
    centroids <- getNewCentroids(expr.data,new.clusters)
    if(!silent)
      cat("We got",ncol(centroids),"new centroids\n")


    iteration = iteration + 1
    allgchanges = c(allgchanges,exchanged.genes)
  }

  if(!silent)
    cat("We finish with",(iteration-1),"iterations\n")
  if(iteration > 1){
    if(!silent)
      cat("Last number of gene changes where",exchanged.genes,"\n")
    #saveRDS(partitions,partitions.file)

    net = NULL
    net$moduleColors = WGCNA::labels2colors(new.clusters)
    names(net$moduleColors) = geneNames
    names(net)
    net$MEs = WGCNA::moduleEigengenes(expr.data,
                                      net$moduleColors,
                                      excludeGrey=excludeGrey)$eigengenes

    net$partitions = partitions
    net$cgenes = allgchanges
  }

  print("The k-means algorithm finished correctly")
  return(net)
}



getNewCentroids <- function(expr.data,partition.in.colors){
  eg.vectors = WGCNA::moduleEigengenes(expr.data,
                                       partition.in.colors,
                                       excludeGrey=F)$eigengenes

  names(eg.vectors) <- substring(names(eg.vectors),3)
  #eg.vectors <- eg.vectors[,centroid.labels]
  return(eg.vectors)
}



getBestModuleCor = function(gene,centroids,signed=TRUE,cor.type){
  return(which.max(corDistance(a=centroids,b=gene,signed=signed,cor.type=cor.type)))
}

createCentroidMatrix <- function(eigengenes){
  my.matrix <- NULL
  for(eigengene in eigengenes){
    my.matrix <- cbind(my.matrix,eigengene)
  }
  colnames(my.matrix) = substring(names(eigengenes),3)
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

#' Title Create a WGCNA network + TOM + some plots
#' If you are not interested in the k-means or before getting into that you want to
#' have a look at the basic WGCNA network use this method
#'
#' @param expr.data Aan object or a file, genes at columns, samples at rows.
#' The network net$moduleColors vector will be generated as genes are in the data.frame.
#' @param beta If beta is < 0 then it is obtained automatically. You can set your own.
#' @param net.type Either "signed" or "unsigned" for easier or trickier MM values interpretation
#' @param tissue.name The label you want the code to use for naming the network
#' @param title For possible plots
#' @param additional.prefix Something you want to add to the naming of files generated to avoid
#' potential clashes
#' @param min.cluster.size Minimum number of genes for a group of them to be a cluster
#' @param save.plots You want plots?
#' @param return.tom Set this to TRUE if you want the TOM returned with the rest of stuff (it may be really big)
#' @param excludeGrey Discard grey genes from the result network
#' @param max.k.cutoff Connectivity parameter
#' @param r.sq.cutoff Only beta values above this value for the adjusted linear regression model are
#' considered
#' @param cor.type Values in "pearson", "kendall", "spearman"
#'
#' @return If all goes well, a network with the following elements (as a list)
#' * MEs will be a matrix with the eigengenes (columns) and the values for each sample (rows)
#' * moduleLabels A vector of labels corresponding to into which cluster each gene is.
#' labels appear in the same order as genes were disposed in the expression data matrix
#' * moduleColors A vector of colors corresponding to into which cluster each gene is. Useful to
#' refer to each module with a color and to signal them at plots. The vector is named with the
#' genes and they appear in the same order as genes were disposed in the expression data matrix
#' beta Is the beta value used
#' type The type of network as in "signed" or "unsigned"
#' geneTree The dendrogram from which the partition was generated
#' @export
#'
#' @examples
getAndPlotNetworkLong <- function(expr.data,beta,
                                  net.type="signed",
                                  tissue.name="MyTissue",
                                  title=NULL,
                                  additional.prefix=NULL,
                                  min.cluster.size=100,
                                  save.plots=TRUE,
                                  return.tom=FALSE,
                                  excludeGrey=FALSE,
                                  max.k.cutoff = 150,
                                  r.sq.cutoff = 0.8,
                                  cor.type="pearson",
                                  silent=T){

  net <- NULL

  if(typeof(expr.data) == "character"){
    if(!silent)
      print(paste0("Reading expression from ",expr.data))
    expr.data = readRDS(expr.data)
  }

  if(!silent)
    print(paste0("We called getAndPlotNetworkLong with ",ncol(expr.data),
                 " genes and ",nrow(expr.data)," samples and beta ",beta,
                 " and correlation type ",cor.type,
                 " and network type ", net.type," and min.cluster.size ",
                 min.cluster.size, " for tissue ",tissue.name))

  if(!silent)
    print(paste0("Expression data is the following within getAndPlotNetworkLong"))
  if(!silent)
    print(expr.data[1:5,1:5])

  #We assume gene names are at columns
  gene.names <- colnames(expr.data)
  sample.names <- rownames(expr.data)
  #Lets delete unused memory
  if(!silent)
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


    if(!silent)
      cat("Choosing beta from SFT.R.sq cut-off",r.sq.cutoff,"and max.k cut-off",max.k.cutoff,"\n")
    beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                    as.numeric(b.study$fitIndices$slope) < 0 &
                                    as.numeric(b.study$fitIndices$max.k) > max.k.cutoff,"Power"])

    if(beta == Inf){
      #OK, lets drop-off the maxk.cutoff
      beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                      as.numeric(b.study$fitIndices$slope) < 0,"Power"])
    }

    if(!silent)
      print(paste0("The estimated beta value is ",beta))
    if(!silent)
      print(paste0("The suggested beta was ",b.study$powerEstimate))

    if(beta == -Inf & is.na(b.study$powerEstimate)){
      stop("There is something wrong, our beta is",beta,"and suggested is",
           b.study$powerEstimate,"\n")
    }

    if(is.na(b.study$powerEstimate)){
      if(!silent)
        cat("We'll use",beta,"for beta\n")

    }else{
      if(beta == -Inf){
        if(!silent)
          cat("We'll use WGCNA's suggested beta\n")
        beta = b.study$powerEstimate
      }else if(beta - b.study$powerEstimate > 5){
        beta = trunc(0.5*(beta + b.study$powerEstimate))
        if(!silent)
          cat("We'll use average between WGCNA's suggested beta and ours.\n")
      }
    }

    if(beta == Inf){
      beta = 21
      if(!silent)
        cat("Warning, the final beta is",beta,"and SFT is compromised\n")
    }

    cat("The final beta value to use is:",beta,"\n")
  }

  additional.prefix = paste0(additional.prefix,tissue.name,".",beta)

  #Create the adjacency matrix of genes coming from the expression data, with the beta
  #passwd as argument
  if(!silent)
    print("Creating adjacency matrix")
  if(cor.type == "spearman")
    corOptions = "use = 'p', method = 'spearman'"
  else
    corOptions = "use = 'p'"

  adjacency = WGCNA::adjacency(expr.data, power = beta, type = net.type, corOptions = corOptions)
  if(!silent)
    print("Created")
  if(!silent)
    print(paste0("Adjacency is the following within getAndPlotNetworkLong"))
  if(!silent)
    print(adjacency[1:5,1:5])
  # Topological Overlap Matrix (TOM)
  # Turn adjacency into topological overlap
  if(!silent)
    print("Creating TOM")
  if(net.type == "signed")
    TOM = WGCNA::TOMsimilarity(adjacency)
  else if(net.type == "unsigned"){
    if(!silent)
      cat("As the network type is unsigned, the TOM type we'll create is signed")
    TOM = WGCNA::TOMsimilarity(adjacency,TOMType="signed")
  }else{
    stop(paste0("Nework type ",net.type," unrecognized when creating the network"))
  }
  #Now we can delete adjacency
  if(!silent)
    print("Deleting adjacency matrix")
  rm(adjacency)
  adjacency = apply(TOM,2,sum)
  adjacency = adjacency/length(adjacency)
  names(adjacency) = colnames(expr.data)

  dissTOM = 1-TOM
  rm(TOM)
  if(!silent)
    print("Created TOM, dissTOM")
  if(!silent)
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

  if(n.mods <= 3){
    print(table(dynamicMods))
    stop(paste0("There are only ",n.mods," modules in the network. Makes no sense to continue"))
  }
  dynamicColors = WGCNA::labels2colors(dynamicMods)
  if(!silent)
    print(table(dynamicColors))
  dynamicColors = CoExpNets::dropGreyModule(dynamicColors)
  if(!silent)
    print(tissue.name)
  if(!silent)
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
  if(!silent)
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


  #Prepare for creating the network objecto to return
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", WGCNA::standardColors(400) )
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  tb2 <- table(moduleColors)[order(table(moduleColors))]
  if(!silent)
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
  net$adjacency = adjacency
  net$beta = beta
  net$type = net.type
  names(net$moduleColors) <- gene.names
  net$geneTree <- geneTree
  if(return.tom){
    return(list(net=net,tom=(1 - dissTOM)))
  }else
    return(net)
}

getAndPlotSCNoTOMNetworkLong <- function(expr.data,beta,
                                         net.type="signed",
                                         tissue.name="MyTissue",
                                         title=NULL,
                                         additional.prefix=NULL,
                                         min.cluster.size=100,
                                         save.plots=TRUE,
                                         return.tom=FALSE,
                                         excludeGrey=FALSE,
                                         max.k.cutoff = 150,
                                         r.sq.cutoff = 0.8,
                                         cor.type="pearson",
                                         silent=T){

  net <- NULL

  if(typeof(expr.data) == "character"){
    if(!silent)
      print(paste0("Reading expression from ",expr.data))
    expr.data = readRDS(expr.data)
  }

  if(!silent)
    print(paste0("We called getAndPlotNetworkLong with ",ncol(expr.data),
                 " genes and ",nrow(expr.data)," samples and beta ",beta,
                 " and correlation type ",cor.type,
                 " and network type ", net.type," and min.cluster.size ",
                 min.cluster.size, " for tissue ",tissue.name))

  if(!silent)
    print(paste0("Expression data is the following within getAndPlotNetworkLong"))
  if(!silent)
    print(expr.data[1:5,1:5])

  #We assume gene names are at columns
  gene.names <- colnames(expr.data)
  sample.names <- rownames(expr.data)
  #Lets delete unused memory
  if(!silent)
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


    if(!silent)
      cat("Choosing beta from SFT.R.sq cut-off",r.sq.cutoff,"and max.k cut-off",max.k.cutoff,"\n")
    beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                    as.numeric(b.study$fitIndices$slope) < 0 &
                                    as.numeric(b.study$fitIndices$max.k) > max.k.cutoff,"Power"])

    if(beta == Inf){
      #OK, lets drop-off the maxk.cutoff
      beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                      as.numeric(b.study$fitIndices$slope) < 0,"Power"])
    }

    if(!silent)
      print(paste0("The estimated beta value is ",beta))
    if(!silent)
      print(paste0("The suggested beta was ",b.study$powerEstimate))

    if(beta == -Inf & is.na(b.study$powerEstimate)){
      stop("There is something wrong, our beta is",beta,"and suggested is",
           b.study$powerEstimate,"\n")
    }

    if(is.na(b.study$powerEstimate)){
      if(!silent)
        cat("We'll use",beta,"for beta\n")

    }else{
      if(beta == -Inf){
        if(!silent)
          cat("We'll use WGCNA's suggested beta\n")
        beta = b.study$powerEstimate
      }else if(beta - b.study$powerEstimate > 5){
        beta = trunc(0.5*(beta + b.study$powerEstimate))
        if(!silent)
          cat("We'll use average between WGCNA's suggested beta and ours.\n")
      }
    }

    if(beta == Inf){
      beta = 21
      if(!silent)
        cat("Warning, the final beta is",beta,"and SFT is compromised\n")
    }

    cat("The final beta value to use is:",beta,"\n")
  }
  additional.prefix = paste0(additional.prefix,tissue.name,".",beta)

  #Create the adjacency matrix of genes coming from the expression data, with the beta
  #passwd as argument
  if(!silent)
    print("Creating adjacency matrix")
  if(cor.type == "spearman")
    corOptions = "use = 'p', method = 'spearman'"
  else
    corOptions = "use = 'p'"

  adjacency = WGCNA::adjacency(expr.data, power = beta,
                               type = net.type, corOptions = corOptions)

  dissTOM = 1-adjacency
  if(!silent)
    print("Created TOM, dissTOM")
  if(!silent)
    print(dissTOM[1:5,1:5])

  geneTree = flashClust::flashClust(as.dist(dissTOM), method = "average")

  print("Now the genetree")
  print(geneTree)
  # Dynamic Tree Cut
  # We like large modules, so we set the minimum module size relatively high
  # Module identification using dynamic tree cut

  n.mods = 0
  deep.split = 2
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,
                                              method="tree",
                                              distM = dissTOM,
                                              deepSplit = TRUE,
                                              minClusterSize = min.cluster.size)
  n.mods = length(table(dynamicMods))
  dynamicColors = WGCNA::labels2colors(dynamicMods)
  if(n.mods == 1 | (n.mods == 2 & "grey" %in% dynamicColors)){
    print(table(dynamicMods))
    return(NULL)
  }

  if(!silent)
    print(table(dynamicColors))
  dynamicColors = CoExpNets::dropGreyModule(dynamicColors)
  if(!silent)
    print(tissue.name)
  if(!silent)
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

  if(save.plots){
    dendro.name = paste0(additional.prefix,"_dendro_colors.pdf")
    eigengenes.name = paste0(additional.prefix,"_Eigengenes_clustering.pdf")

    pdf(file=dendro.name)
    WGCNA::plotDendroAndColors(geneTree, dynamicColors,
                               "Dynamic Tree Cut",
                               dendroLabels = FALSE,
                               hang = 0.03, addGuide = TRUE,
                               guideHang = 0.05,main=title)
    dev.off()

    pdf(file=eigengenes.name)
    # Plot the result
    tb = dynamicColors
    tb <- tb[ METree$labels ]
    METree$labels <- paste0(names(tb), ":", tb)

    plot(METree, main = paste0(title," ",tissue.name,
                               " Clustering of module eigengenes"),
         xlab = "", sub = "")
    dev.off()
  }


  #Prepare for creating the network objecto to return
  # Rename to moduleColors
  moduleColors = dynamicColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", WGCNA::standardColors(400) )
  moduleLabels = match(moduleColors, colorOrder)-1
  tb2 <- table(moduleColors)[order(table(moduleColors))]
  if(!silent)
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
  adjacency = apply(adjacency,2,sum)
  adjacency = adjacency/length(adjacency)
  names(adjacency) = colnames(expr.data)

  net$MEs <- MEs
  rownames(net$MEs) <- sample.names
  net$moduleLabels <- moduleLabels
  net$moduleColors <- moduleColors
  net$adjacency = adjacency
  net$beta = beta
  net$type = net.type
  names(net$moduleColors) <- gene.names
  net$geneTree <- geneTree
  return(net)
}

#' Title Create a WGCNA network + TOM + some plots
#' If you are not interested in the k-means or before getting into that you want to
#' have a look at the basic WGCNA network use this method
#'
#' @param expr.data Aan object or a file, genes at columns, samples at rows.
#' The network net$moduleColors vector will be generated as genes are in the data.frame.
#' @param beta If beta is < 0 then it is obtained automatically. You can set your own.
#' @param net.type Either "signed" or "unsigned" for easier or trickier MM values interpretation
#' @param tissue.name The label you want the code to use for naming the network
#' @param title For possible plots
#' @param additional.prefix Something you want to add to the naming of files generated to avoid
#' potential clashes
#' @param min.cluster.size Minimum number of genes for a group of them to be a cluster
#' @param save.plots You want plots?
#' @param return.tom Set this to TRUE if you want the TOM returned with the rest of stuff (it may be really big)
#' @param excludeGrey Discard grey genes from the result network
#' @param max.k.cutoff Connectivity parameter
#' @param r.sq.cutoff Only beta values above this value for the adjusted linear regression model are
#' considered
#' @param cor.type Values in "pearson", "kendall", "spearman"
#'
#' @return If all goes well, a network with the following elements (as a list)
#' * MEs will be a matrix with the eigengenes (columns) and the values for each sample (rows)
#' * moduleLabels A vector of labels corresponding to into which cluster each gene is.
#' labels appear in the same order as genes were disposed in the expression data matrix
#' * moduleColors A vector of colors corresponding to into which cluster each gene is. Useful to
#' refer to each module with a color and to signal them at plots. The vector is named with the
#' genes and they appear in the same order as genes were disposed in the expression data matrix
#' beta Is the beta value used
#' type The type of network as in "signed" or "unsigned"
#' geneTree The dendrogram from which the partition was generated
#' @export
#'
#' @examples
getAndPlotSCNetworkLong <- function(expr.data,beta,
                                    net.type="scell",
                                    tissue.name="MyTissue",
                                    title=NULL,
                                    additional.prefix=NULL,
                                    min.cluster.size=100,
                                    save.plots=TRUE,
                                    excludeGrey=FALSE,
                                    max.k.cutoff = 150,
                                    cor.type="pearson",
                                    silent=T){

  net <- NULL

  if(typeof(expr.data) == "character"){
    if(!silent)
      print(paste0("Reading expression from ",expr.data))
    expr.data = readRDS(expr.data)
  }

  if(!silent)
    print(paste0("We called getAndPlotSCNetworkLong with ",ncol(expr.data),
                 " genes and ",nrow(expr.data)," samples and network type scell and min.cluster.size ",
                 min.cluster.size, " for tissue ",tissue.name))

  if(!silent)
    print(paste0("Expression data is the following within getAndPlotSCNetworkLong"))
  if(!silent)
    print(expr.data[1:5,1:5])

  #We assume gene names are at columns
  gene.names <- colnames(expr.data)
  sample.names <- rownames(expr.data)
  #Lets delete unused memory
  if(!silent)
    print("Garbage collecting")
  gc()

  cat("Creating now clustering with expression\n")
  expr.data.dist = dist(t(expr.data))
  geneTree = flashClust::flashClust(expr.data.dist, method = "average")

  if(!silent)
    print("Created")

  print("Now the genetree")
  print(geneTree)
  # Dynamic Tree Cut
  # We like large modules, so we set the minimum module size relatively high
  # Module identification using dynamic tree cut

  n.mods = 0
  deep.split = 2
  while(n.mods < 10 & deep.split < 5){
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,
                                                distM = as.matrix(expr.data.dist),
                                                deepSplit = deep.split,
                                                pamRespectsDendro = F,
                                                respectSmallClusters = F,
                                                minClusterSize = min.cluster.size)
    n.mods = length(table(dynamicMods))
    deep.split = deep.split + 1
  }

  if(n.mods <= 3){
    print(table(dynamicMods))
    stop(paste0("There are only ",n.mods," modules in the network. Makes no sense to continue"))
  }
  dynamicColors = WGCNA::labels2colors(dynamicMods)
  if(!silent)
    print(table(dynamicColors))
  dynamicColors = CoExpNets::dropGreyModule(dynamicColors)
  if(!silent)
    print(tissue.name)
  if(!silent)
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



  #Prepare for creating the network objecto to return
  # Rename to moduleColors
  moduleColors = dynamicColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", WGCNA::standardColors(400) )
  moduleLabels = match(moduleColors, colorOrder)-1
  tb2 <- table(moduleColors)[order(table(moduleColors))]
  if(!silent)
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
  net$adjacency = NA
  net$beta = NA
  net$type = net.type
  names(net$moduleColors) <- gene.names
  net$geneTree <- geneTree
  return(net)
}

corDistance = function(a,b,signed=TRUE,cor.type="pearson"){
  if(cor.type=="pearson"){
    if(signed)
      #return(0.5 + 0.5*WGCNA::corFast(x=a,y=b)) #(Note they are equivalent)
      return(0.5 * (1 + stats::cor(a,b,use="pairwise.complete.obs")))
    return(abs(stats::cor(a,b)))
  }else{
    if(signed)
      #return(0.5 + 0.5*WGCNA::corFast(a,b)) #(Note they are equivalent)
      return(0.5 * (1 + stats::cor(a,b,method=cor.type)))
    return(abs(stats::cor(a,b,method=cor.type,use="pairwise.complete.obs")))
  }
}


#' Title Obtain the Module Membership of genes
#' This function works by getting the Pearson correlation between the
#' eigengene of the module and the gene(s).
#'
#' @param net The network can be either and object or a file
#' @param expr.data.file Again an object or a file, genes at columns, samples at rows.
#' It is expected that the genes are in the same order as they are expressed in the
#' network net$moduleColors vector. Use this parameter together with the net one and
#' which.one set to "new" for a network which is not in the DDBB.
#' @param tissue A tissue as it can be find in the network DDBB
#' @param genes The gene names. If set to NULL, the MM of all genes in the network is obtained
#' @param which.one The category the network belongs to. If not in the DDBB just ignore it
#' @param silent Set to true if you want no log
#' @param keep.grey Use grey genes too
#' @param alt.gene.index You can pass genes in a different order, use this index order for that
#'
#' @returnk Depending on the value of table.format, a data.frame with genes and MM values or a list
#' with genes as keys and the MMs as values
#' @export
#'
#' @examples
getMM = function(net=NULL,
                 expr.data.file=NULL,
                 tissue,
                 genes,
                 which.one="rnaseq",
                 silent=F,
                 keep.grey=F,
                 identicalNames=T, #When gene ids in expr data and net are identical
                 alt.gene.index=NULL,
                 dupAware=T){

  if(is.null(net))
    net = getNetworkFromTissue(tissue,which.one)
  else{
    if(typeof(net) == "character")
      net = readRDS(net)
  }
  if(is.null(expr.data.file))
    expr.data = getExprDataFromTissue(tissue=tissue,which.one=which.one,only.file = F)
  else{
    if(typeof(expr.data.file) == "character")
      expr.data = readRDS(expr.data.file)
    else
      expr.data = expr.data.file
  }

  if(!identicalNames){
    colnames(expr.data) = fromAny2Ensembl(colnames(expr.data))
    names(net$moduleColors) = fromAny2Ensembl(names(net$moduleColors))
    if(is.null(genes))
      genes = names(net$moduleColors)
    ens.genes = fromAny2Ensembl(genes)
  }else{
    if(is.null(genes))
      genes = names(net$moduleColors)
    ens.genes = genes
  }


  if(is.null(expr.data)){
    cat("There is no expr.data file registered for category", which.one,"and tissue",tissue,"\n")
    return(expr.data)
  }

  mm = NULL

  if(is.null(genes)){
    #There is no correlation within grey module normally
    genes = names(net$moduleColors) #[net$moduleColors != "grey"]
  }

  if(dupAware){
    newgenes = NULL
    mods = NULL
    for(gene in ens.genes){
      localmods = net$moduleColors[names(net$moduleColors) == gene]
      if(length(localmods)){
        newgenes = c(newgenes,rep(gene,length(localmods)))
        mods = c(mods,localmods)
      }
    }
    genes = newgenes
    modules = mods
    finalgenes = NULL
    finalmodules = NULL
    for(module in unique(modules)){
      genesinmodule = genes[modules == module]
      finalgenes = c(finalgenes,genesinmodule[!duplicated(genesinmodule)])
      finalmodules = c(finalmodules,rep(module,sum(!duplicated(genesinmodule))))
    }
    genes = finalgenes
    modules = finalmodules

  }else
    modules = net$moduleColors[match(genes,names(net$moduleColors))]

  if(is.null(genes))
    return(NULL)

  out.table = data.frame(list(ensgene=character(0),name=character(0),
                              module=character(0),mm=numeric(0)),
                         stringsAsFactors=FALSE)
  n.row = 1
  out.table[1:length(genes),1] = genes
  out.table[1:length(genes),2] = fromAny2GeneName(genes)
  out.table[1:length(genes),3] = modules



  for(i in 1:length(genes)){
    gene = genes[i]
    module = modules[i]
    expr.data.gene.index = gene

    if(length(module) == 0){
      if(!silent)
        cat("Gene ",gene," not in network\n")
      out.table[i,4] = -1
    }else{
      if(module == "grey" & !keep.grey)
        out.table[i,4] = 0
      else{

        tryCatch(out.table[i,4] <- cor(net$MEs[paste0("ME",module)],expr.data[,expr.data.gene.index]),
                 error = function(e){
                   print(paste0("Error in module ",module," ",e))
                 })
      }
    }
  }
  return(out.table)
}



corWithCatTraits = function(tissue,which.one,covlist,covs=NULL,retPVals=F){
  if(is.null(covs))
    covs = getCovariates(tissue=tissue,which.one=which.one)

  if(!is.null(covlist))
    covs = covs[,covlist,drop=F]

  for(i in 1:ncol(covs)){
    if(typeof(covs[,i]) ==  "character")
      covs[,i] = as.factor(covs[,i])
  }
  factor.mask = unlist(lapply(covs,is.factor))
  cat("We will work with",sum(factor.mask),"factors\n")
  stopifnot(sum(factor.mask) > 0)

  MEs = getNetworkEigengenes(tissue=tissue,which.one=which.one)


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

  if(sum(!factor.mask) > 0){
    moduleTraitCor = cor(MEs,covs[,!factor.mask,drop=FALSE],use="p")
    #Generate the p-values for significance of a given matrix of correlations, for all modules,
    #between traits data and eigengenes, both from samples
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nrow(MEs))
    moduleTraitPvalue = cbind(moduleTraitPvalue,fcm)
    colnames(moduleTraitPvalue) = c(colnames(covs)[!factor.mask],
                                    colnames(covs)[factor.mask])
  }else{
    moduleTraitPvalue = fcm
    colnames(moduleTraitPvalue) = colnames(covs)[factor.mask]
  }
  if(retPVals)
    toReturn = moduleTraitPvalue
  else
    toReturn = -log10(moduleTraitPvalue)
  rownames(toReturn) = gsub("ME","",names(MEs))
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
  return(toReturn)
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

#' Title Managing big TOM matrices
#' This function allows saving a bit TOM matrix aquared matrix, keeping only those
#' cells referring to genes clustering together in the same module
#'
#' @param tom The TOM matrix itself
#' @param clusters A vector with the labels (in the same order as in the expression data)
#' showing how genes cluster together
#' @param filepref This full path including file prefix will be used for the TOM files to
#' be saved
#'
#' @return No value returned
#' @export
#'
#' @examples
saveTOM = function(tom, clusters, filepref){
  print("Hello saveTOM!!!!")
  print(clusters)
  modules = unique(clusters)
  metadata = NULL
  lapply(modules,function(x){
    print(paste0("saving",x))
    mask = clusters %in% x
    tomname = paste0(filepref,"_",x,".rds")
    saveRDS(tom[mask,mask],tomname)
    metadata[[x]] <<- list(mask=mask,tomname=tomname)
  })
  saveRDS(metadata,paste0(filepref,"_metadata.rds"))
}

#' Title Managing big TOM matrices
#'
#' @param filepref The prefix we used to previously save the matrix
#'
#' @return A squared TOM matrix disposed in the same order as genes in the
#' expression data matrix
#' @export
#'
#' @examples
readTOM = function(filepref){
  cat("Hello!!!!!")
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

incrementTOM = function(filepref){
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

#' Title Managing big TOM matrices
#'
#' @param filepref
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param prince
#' @param label
#' @param smallest
#' @param note
#' @param notecol
#' @param notecex
#' @param breaks
#' @param col
#' @param margins
#' @param key
#' @param cexRow
#' @param cexCol
#' @param xlab
#' @param colsep
#' @param rowsep
#' @param sepcolor
#' @param sepwidth
#' @param Rsquared
#' @param breaksRsquared
#' @param main
#'
#' @return
#' @export
#'
#' @examples
princePlot = function (prince, label = colnames(prince$o), smallest = -20,
                       note = F, notecol = "black", notecex = 1,
                       breaks = seq(-20, 0, length.out = 100),
                       col = heat.colors(99), margins = c(5,
                                                          7), key = T, cexRow = 1, cexCol = 1, xlab = "Principal Components (Variation)",
                       colsep = NULL, rowsep = NULL, sepcolor = "black",
                       sepwidth = c(0.05, 0.05),
                       Rsquared = F, breaksRsquared = seq(0, 1, length.out = 100),main)
{
  if (class(prince) != "prince") {
    stop("prince is not an object generated by the prince function")
  }
  if (smallest > 0) {
    stop("smallest has to be less than 0")
  }
  require(gplots)
  linp10 <- log10(prince$linp)
  linp10 <- replace(linp10, linp10 <= smallest, smallest)
  tonote <- signif(prince$linp, 1)
  prop <- round(prince$prop, 0)
  if (Rsquared == T) {
    linp10 <- prince$rsquared
    breaks <- breaksRsquared
    tonote <- round(prince$rsquared, 2)
    col <- col[length(col):1]
  }
  heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none",
            trace = "none", symbreaks = F, symkey = F, breaks = breaks,
            key = key, col = col, cexRow = cexRow, cexCol = cexCol,
            colsep = colsep, rowsep = rowsep, sepcolor = sepcolor,
            sepwidth = sepwidth, main = main, labCol = paste(1:ncol(linp10),
                                                             " (", prop, ")", sep = ""), margins = margins, labRow = label,
            xlab = xlab, cellnote = if (note == T) {
              tonote
            }
            else {
              matrix(ncol = ncol(prince$linp), nrow = nrow(prince$linp))
            }, notecol = notecol, notecex = notecex)
}

adjacencyGeneration = function(beta,which.one,package=NULL,from.name=T){

  if(is.null(package))
    package = which.one

  #Adjacency for microarray
  ts = getAvailableNetworks(which.one)
  #ts = ts[(1+which(ts == "Testis")):length(ts)]
  for(tissue in ts){
    cat("Working on tissue's adjacency",tissue,"\n")
    netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one=which.one,only.file=T)
    if(from.name){
      if(which.one == "CoExpROSMAP" | which.one == "gtexv6"){
        beta = as.numeric(stringr::str_split(basename(CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                                                      which.one=which.one,
                                                                                      only.file=T)),
                                             "\\.")[[1]][2])
      }else if(which.one == "10UKBEC"){
        beta = 12
      }else if(which.one == "nabec"){
        beta = 8
      }else
        stop(paste0("Unknown category", which.one))
    }
    cat("Using beta value",beta,"\n")
    adjfile = paste0(netFile,".adj.rds")
    expr.data = CoExpNets::getExprDataFromTissue(tissue=tissue,which.one=which.one)
    adjacency = WGCNA::adjacency(expr.data, power = beta, type = "signed" )
    TOM = WGCNA::TOMsimilarity(adjacency)
    colnames(TOM) = colnames(expr.data)
    localadj = apply(TOM,2,sum)
    names(localadj) = colnames(TOM)
    cat("Saving at",adjfile,"\n")
    saveRDS(localadj,adjfile)
  }
}

mmGeneration = function(which.one,package=NULL){

  if(is.null(package))
    package = which.one

  #AMM for microarray
  ts = getAvailableNetworks(which.one)
  #ts = ts[(1+which(ts == "Testis")):length(ts)]
  for(tissue in ts){
    cat("Working on tissue's MM for",tissue,"\n")
    netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one=which.one,only.file=T)
    mm = CoExpNets::getMM(tissue=tissue,which.one=which.one,genes=NULL)
    mmfile = paste0(netFile,".mm.rds")
    cat("Saving at",mmfile,"\n")
    saveRDS(mm,mmfile)
  }
}

adjacencyStudy = function(categories = c("CoExpROSMAP",
                                         "10UKBEC","gtexv6"),
                          common.genes=F){
  batchid = NULL
  nsamples = NULL
  betas = NULL
  if(common.genes){
    tcount = 0
    tnames = NULL
    allgenes = NULL
    for(category in categories){
      cat("On category",category,"now\n")
      ts = getAvailableNetworks(category=category)
      if(length(ts))
        for(tissue in ts){

          cat("Working category",category," and on tissue's adjacency",tissue,"\n")
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          adjfile = paste0(netFile,".adj.rds")
          net = readRDS(netFile)

          if(file.exists(adjfile)){
            genes = fromAny2GeneName(names(net$moduleColors))

            if(tcount == 0)
              allgenes = genes
            else
              allgenes = intersect(allgenes,genes)

            print("########")
            cat("Genes are",length(allgenes),"\n")
            print("########")
            tcount = tcount + 1
            nsamples = c(nsamples,nrow(net$MEs))
            batchid = c(batchid,category)
            tnames = c(tnames,paste0(category,"_",tissue))

            if(category == "nabec" | category == "CoExpROSMAP" | category == "gtexv6" | category == "CoExpGTExV7"){
              beta = as.numeric(stringr::str_split(basename(CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                                                            which.one=category,
                                                                                            only.file=T)),
                                                   "\\.")[[1]][2])
            }else if(which.one == "nabec"){
              beta = 8
            }else if(category == "10UKBEC"){
              beta = 12
            }
            betas = c(betas,beta)
          }
        }
    }
    adjs = matrix(nrow=tcount,ncol=length(allgenes))
    rownames(adjs) = tnames
    colnames(adjs) = allgenes
    adjs[] = 0
    tcount = 0
    for(category in categories){
      ts = getAvailableNetworks(category=category)

      if(length(ts))
        for(tissue in ts){
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          adjfile = paste0(netFile,".adj.rds")
          if(file.exists(adjfile)){
            tcount = tcount + 1
            adjgenes = readRDS(adjfile)
            if(category == "CoExpROSMAP")
              names(adjgenes) = unlist(lapply(names(adjgenes),function(x){
                stringr::str_split(x,"\\.")[[1]][1]}))
            names(adjgenes) = fromAny2GeneName(names(adjgenes))
            mask = names(adjgenes) %in% colnames(adjs)
            adjgenes = adjgenes[mask]
            adjs[tcount,] = adjgenes[match(colnames(adjs),names(adjgenes))]
          }
        }
    }



  }else{
    tcount = 0
    tnames = NULL
    allgenes = NULL

    for(category in categories){
      cat("On category",category,"now\n")
      ts = getAvailableNetworks(category=category)
      if(length(ts))
        for(tissue in ts){

          cat("Working category",category," and on tissue's adjacency",tissue,"\n")
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          adjfile = paste0(netFile,".adj.rds")
          if(file.exists(adjfile)){
            allgenes = c(allgenes,
                         fromAny2GeneName(names(CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                                                which.one=category)$moduleColors)))
            tcount = tcount + 1
            tnames = c(tnames,paste0(category,"_",tissue))
            batchid = c(batchid,category)
            if(category == "nabec" | category == "CoExpROSMAP" | category == "gtexv6" | category == "CoExpGTExV7"){
              beta = as.numeric(stringr::str_split(basename(CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                                                            which.one=category,
                                                                                            only.file=T)),
                                                   "\\.")[[1]][2])
            }else if(which.one == "nabec"){
              beta = 8
            }else if(category == "10UKBEC"){
              beta = 12
            }
            betas = c(betas,beta)

          }

        }
    }
    allgenes = unique(allgenes)
    adjs = matrix(nrow=tcount,ncol=length(allgenes))
    rownames(adjs) = tnames
    colnames(adjs) = allgenes
    adjs[] = 0

    tcount = 0
    for(category in categories){
      ts = getAvailableNetworks(category=category)

      if(length(ts))
        for(tissue in ts){
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          adjfile = paste0(netFile,".adj.rds")
          if(file.exists(adjfile)){
            tcount = tcount + 1

            adjgenes = readRDS(adjfile)
            names(adjgenes) = fromAny2GeneName(names(adjgenes))
            mask = colnames(adjs) %in% names(adjgenes)
            adjgenes = adjgenes[names(adjgenes) %in% colnames(adjs)]
            adjs[tcount,mask] = adjgenes[match(colnames(adjs)[mask],names(adjgenes))]
          }
        }
    }
  }

  #adjn = normalize.quantiles(adjs)
  resids = adjs
  tissues = NULL
  for(i in 1:length(categories)){
    localbatch = batchid
    localbatch[localbatch != categories[i]] = "theother"
    localbatch = as.numeric(as.factor(localbatch))
    tissues = cbind(tissues,localbatch)
  }
  tissues = cbind(tissues,betas)
  tissues = cbind(tissues,nsamples)
  tissues = scale(tissues)
  tissues = as.data.frame(tissues)
  colnames(tissues) = c(categories,"beta","nsamples")

  resids = apply(resids, 2, function(y){
    lm( y ~ . , data=tissues)$residuals
  })

  stsne = Rtsne(resids,dims=2,perplexity=15,check_duplicates = F)
  plot(stsne$Y)
  text(x=stsne$Y[,1],y=stsne$Y[,2],labels=gsub("gtexv6_|10UKBEC|CoExpROSMAP","",
                                               rownames(adjs)),cex=0.5)
  rownames(resids) = rownames(adjs)
  return(list(resids = resids,categories = batchid))
}

prettyPlot = function(mms,what="PCA",zoom=NULL){
  labels = rownames(mms$resids)
  colors = labels
  colors[grep("gtexv6",labels)] = "GTExV6"
  colors[grep("CoExpROSMAP_",labels)] = "ROSMAP"
  colors[grep("10UKBEC_",labels)] = "10UKBEC"
  colors[grep("nabec",labels)] = "NABEC"

  labels = gsub("gtexv6_|10UKBEC_|CoExpROSMAP_","",labels)
  labels[labels == "notad"] = "Not_AD"
  labels[labels == "probad"] = "Probably_AD"
  labels[labels == "ad"] = "AD"
  labels[labels == "allsamples"] = "ROSMAP"
  labels[labels == "nabec_Cortex"] = "NABEC_Cortex"

  if(what == "tSNE"){
    stsne = Rtsne(mms$resids,dims=2,perplexity=15,check_duplicates = F)

    if(!is.null(zoom)){
      mask = stsne$Y[,1] > zoom[1] & stsne$Y[,1] < zoom[2] &
        stsne$Y[,2] > zoom[3] & stsne$Y[,2] < zoom[4]
      stsne$Y = stsne$Y[mask,]
      labels = labels[mask]
      colors = colors[mask]

    }


    mydata = data.frame(stsne$Y)
    colnames(mydata) = paste0("tSNEaxis",1:ncol(mydata))
    ggplot(data=mydata) + aes(x=tSNEaxis1,y=tSNEaxis2) +
      geom_text(aes(x=tSNEaxis1,y=tSNEaxis2,label=labels,color=colors),size=5) +
      theme_minimal() + theme(legend.position="bottom")
    return(NULL)
  }else if(what == "PCA"){
    if(is.null(zoom)){
      pca = prcomp(scale(mms$resids))
      mydata = data.frame(pca$x)
      colnames(mydata) = paste0("PCA",1:ncol(mydata))
      ggplot(data=mydata) + aes(x=PCA1,y=PCA2) +
        geom_text(aes(x=PCA1,y=PCA2,label=labels,color=colors),size=3) +
        theme_minimal() + theme(legend.position="bottom")

    }else{
      pca = prcomp(scale(mms$resids))
      mydata = data.frame(pca$x)
      colnames(mydata) = paste0("PCA",1:ncol(mydata))
      mask = mydata$PCA1 > zoom[1] & mydata$PCA1 < zoom[2] &
        mydata$PCA2 > zoom[3] & mydata$PCA2 < zoom[4]
      mydata = mydata[mask,]
      labels = labels[mask]
      colors = colors[mask]


      ggplot(data=mydata) + aes(x=PCA1,y=PCA2) +
        geom_text(aes(x=PCA1,y=PCA2,label=labels,color=colors),size=3) +
        theme_minimal() + theme(legend.position="bottom")

    }
  }else if(what == "UMAP"){
    if(is.null(zoom)){
      pca = umap(scale(mms$resids))
      plot(x=pca$layout[,1],y=pca$layout[,2])
      text(x=pca$layout[,1],y=pca$layout[,2],label=labels,color=colors)

    }else{
      pca = prcomp(scale(mms$resids))
      mydata = data.frame(pca$x)
      colnames(mydata) = paste0("PCA",1:ncol(mydata))
      mask = mydata$PCA1 > zoom[1] & mydata$PCA1 < zoom[2] &
        mydata$PCA2 > zoom[3] & mydata$PCA2 < zoom[4]
      mydata = mydata[mask,]
      labels = labels[mask]
      colors = colors[mask]


      ggplot(data=mydata) + aes(x=PCA1,y=PCA2) +
        geom_text(aes(x=PCA1,y=PCA2,label=labels,color=colors),size=3) +
        theme_minimal() + theme(legend.position="bottom")

    }
  }



}

MMStudy = function(categories = c("CoExpROSMAP","10UKBEC","gtexv6","nabec"),
                   common.genes=F){
  batchid = NULL
  betas = NULL
  if(common.genes){
    tcount = 0
    tnames = NULL
    allgenes = NULL
    for(category in categories){
      cat("On category",category,"now\n")
      ts = getAvailableNetworks(category=category)
      if(length(ts))
        for(tissue in ts){

          cat("Working category",category," and on tissue's adjacency",tissue,"\n")
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          mmfile = paste0(netFile,".mm.rds")
          if(file.exists(mmfile)){
            genes = fromAny2GeneName(names(CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                                           which.one=category)$moduleColors))
            if(tcount == 0)
              allgenes = genes
            else
              allgenes = intersect(allgenes,genes)

            print("########")
            cat("Genes are",length(allgenes),"\n")
            print("########")
            tcount = tcount + 1
            batchid = c(batchid,category)
            tnames = c(tnames,paste0(category,"_",tissue))
          }
        }
    }
    mms = matrix(nrow=tcount,ncol=length(allgenes))
    rownames(mms) = tnames
    colnames(mms) = allgenes
    mms[] = 0
    tcount = 0
    for(category in categories){
      ts = getAvailableNetworks(category=category)
      if(length(ts))
        for(tissue in ts){
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          mmfile = paste0(netFile,".mm.rds")
          if(file.exists(mmfile)){
            tcount = tcount + 1
            mmgenes = readRDS(mmfile)
            if(category == "CoExpROSMAP")
              names(mmgenes) = unlist(lapply(names(mmgenes),function(x){
                stringr::str_split(x,"\\.")[[1]][1]}))
            names(mmgenes) = fromAny2GeneName(names(mmgenes))
            mask = names(mmgenes) %in% colnames(mms)
            mmgenes = mmgenes[mask]
            mms[tcount,] = mmgenes[match(colnames(mms),names(mmgenes))]
          }
        }
    }



  }else{
    tcount = 0
    tnames = NULL
    allgenes = NULL

    for(category in categories){
      cat("On category",category,"now\n")
      ts = getAvailableNetworks(category=category)
      ts = ts[ts != "Testis"]
      if(length(ts))
        for(tissue in ts){

          cat("Working category",category," and on tissue's MM",tissue,"\n")
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          mmfile = paste0(netFile,".mm.rds")
          if(file.exists(mmfile)){
            allgenes = c(allgenes,
                         fromAny2GeneName(names(CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                                                which.one=category)$moduleColors)))
            tcount = tcount + 1
            tnames = c(tnames,paste0(category,"_",tissue))
            batchid = c(batchid,category)

          }

        }
    }
    allgenes = unique(allgenes)
    mms = matrix(nrow=tcount,ncol=length(allgenes))
    rownames(mms) = tnames
    colnames(mms) = allgenes
    mms[] = 0

    tcount = 0
    for(category in categories){
      ts = getAvailableNetworks(category=category)
      ts = ts[ts != "Testis"]

      if(length(ts))
        for(tissue in ts){
          netFile = CoExpNets::getNetworkFromTissue(tissue=tissue,
                                                    which.one=category,only.file=T)
          mmfile = paste0(netFile,".mm.rds")
          if(file.exists(mmfile)){
            tcount = tcount + 1

            mmgenes = readRDS(mmfile)
            names(mmgenes) = fromAny2GeneName(names(mmgenes))
            mask = colnames(mms) %in% names(mmgenes)
            mmgenes = mmgenes[names(mmgenes) %in% colnames(mms)]
            mms[tcount,mask] = mmgenes[match(colnames(mms)[mask],names(mmgenes))]
          }
        }
    }
    dn = dimnames(mms)
    mms = normalize.quantiles(mms)
    dimnames(mms) = dn
  }

  resids = mms
  tissues = NULL
  for(i in 1:length(categories)){
    localbatch = batchid
    localbatch[localbatch != categories[i]] = "theother"
    localbatch = as.numeric(as.factor(localbatch))
    tissues = cbind(tissues,localbatch)
  }
  tissues = as.data.frame(tissues)
  colnames(tissues) = categories

  resids = apply(resids, 2, function(y){
    lm( y ~ . , data=tissues)$residuals
  })

  stsne = Rtsne(resids,dims=2,perplexity=15,check_duplicates = F)
  plot(stsne$Y,main="MM based patterns across CoExp networks")
  text(x=stsne$Y[,1],y=stsne$Y[,2],labels=gsub("CoExpGTExV7|gtexv6_|10UKBEC|CoExpROSMAP","",
                                               rownames(mms)),cex=0.5)

  rownames(resids) = rownames(mms)
  return(list(resids = resids,categories = batchid))


}

testCoExpNetworks = function(){
  if(!exists("coexp.nets"))
    return(NULL)

  categories = unique(coexp.nets$which.one)
  for(category in categories){
    tissues = getAvailableNetworks(category=category)
    for(tissue in tissues){
      lapply(getModulesFromTissue(tissue=tissue,which.one=category),
             function(x){ reportOnModule(tissue=tissue,which.one=category,
                                         module=x)})
    }
  }

  cat("Success!!\n")
}


bottomUpNetwork = function(tissue="SNIG",
                           which.one="micro19K",
                           threshold=0,
                           folder="~/tmp/",
                           seed.genes,
                           disease.genes=NULL,
                           loaded=F,
                           permutations=100,
                           include.all.edges=F,width=20,height=12,times=7){


  net = CoExpNets::getNetworkFromTissue(tissue=tissue,which.one=which.one)
  expr.data = CoExpNets::getExprDataFromTissue(tissue=tissue,which.one=which.one)

  cat("Disease genes will be",paste0(disease.genes,collapse=", "),"\n")
  cat("Candidate genes will be",paste0(seed.genes,collapse=", "),"\n")

  if(threshold == 0){
    threshold = times*(length(disease.genes) + length(seed.genes))
    cat("We'll use",threshold,"genes within the plot (",times,"times disease + seed genes)\n")
  }else{
    cat("We'll use",threshold,"genes within the plot \n")
  }

  biotypes = fromEnsembl2Function(fromAny2Ensembl(names(net$moduleColors)))
  names(biotypes) = names(net$moduleColors)
  edge.cols <- c("fromNode", "toNode", "weight", "direction", "fromAltName", "toAltName")
  node.cols <- c("nodeName", "altName", "biotype","module","adjacency","mm","type")

  #Let's prepare the naming of files
  folder = paste0(folder,tissue,".",which.one,".")
  if(include.all.edges)
    folder = paste0(folder,"complete.")
  else
    folder = paste0(folder,"partial.")

  if(!loaded){
    tom.matrix <<- getTOMFromTissue(tissue,which.one,"signed")
    colnames(tom.matrix) <<- names(net$moduleColors)
    rownames(tom.matrix) <<- colnames(tom.matrix)
    adjacencies <<- apply(tom.matrix,2,sum)
    names(adjacencies) <<- colnames(tom.matrix)
    mms <<- getMM(net=net,genes=NULL,expr.data.file=expr.data,tissue=tissue,
                  which.one=which.one,in.cluster=F)
  }else{
    adjacencies = apply(tom.matrix,2,sum)
    names(adjacencies) = colnames(tom.matrix)
    initmms = CoExpNets::getMM(genes=NULL,tissue=tissue,
                               which.one=which.one,identicalNames = T)
    mms = initmms$mm
    names(mms) = initmms$ensgene
  }

  knn.data = NULL
  apriori.genes = unique(c(seed.genes,disease.genes))
  gene.modules = unique(net$moduleColors[names(net$moduleColors) %in% apriori.genes])
  cat("Genes are in modules",gene.modules,"\n")
  gene.names <- names(net$moduleColors)[net$moduleColors %in% gene.modules]
  tom.matrix.small <- tom.matrix[gene.names,gene.names]


  the.order = order(apply(tom.matrix.small,2,sum),decreasing=T)
  tom.matrix.small = tom.matrix.small[the.order,the.order]

  all.seed.genes = apriori.genes
  apriori.genes = apriori.genes[apriori.genes %in% colnames(tom.matrix.small)]
  cat("Genes out of the plot will be",
      paste0(all.seed.genes[!(all.seed.genes %in% apriori.genes)],collapse=", "),"\n")

  if(threshold < length(apriori.genes)){
    stop("We need to include",length(apriori.genes),"but the size limit is lower:",threshold,"\n")
  }
  rounds = 1
  genes.included = ""

  #We'll init knn.data
  for(i in 1:length(apriori.genes))
    knn.data[[apriori.genes[i]]] = character(0)

  #We'll use a length on the rank to stop
  rank.limit = 100
  rank.length = 1
  while(rank.length < rank.limit | length(genes.included) < threshold){
    for(apriori.index in 1:length(apriori.genes)){
      seed.gene = apriori.genes[apriori.index]
      new.gene = colnames(tom.matrix.small)[order(tom.matrix.small[seed.gene,],decreasing=T)][rank.length]
      if(length(genes.included) < threshold & !(new.gene %in% genes.included)){
        genes.included = c(genes.included,new.gene)
        cat("New gene",new.gene,"added to the plot\n")
      }
      knn.data[[seed.gene]] = c(knn.data[[seed.gene]],new.gene)
    }
    rank.length = rank.length + 1
  }


  if(permutations > 0){
    #We need to generate estimates of significance for the selected genes
    #Number of random seed genes to choose
    random.n = length(apriori.genes)
    #p.values for each gene included above
    estimates = vector(mode="numeric",length=length(genes.included) - length(apriori.genes))
    estimates[] = 0
    names(estimates) = genes.included[!(genes.included %in% apriori.genes)]

    #Permutate!!
    for(permutation in 1:permutations){
      rounds = 1
      r.genes.included = NULL
      r.apriori.genes = sample((1:ncol(tom.matrix.small))[-which(colnames(tom.matrix.small) %in% apriori.genes)],
                               random.n)
      while(length(r.genes.included) < threshold){
        for(apriori.index in 1:length(r.apriori.genes)){
          seed.gene = r.apriori.genes[apriori.index]
          new.gene = colnames(tom.matrix.small)[order(tom.matrix.small[seed.gene,],decreasing=T)][rounds]
          if(!(new.gene %in% r.genes.included))
            r.genes.included = c(r.genes.included,new.gene)
        }
        rounds = rounds + 1
      }
      #Now we update estimates
      estimates[names(estimates) %in% r.genes.included] = estimates[names(estimates) %in% r.genes.included] + 1
    }
    estimates = estimates/permutations
    names(estimates) = fromAny2GeneName(names(estimates))
  }

  #Now we prepare the tom for generating the network
  tom.matrix.small = tom.matrix.small[rownames(tom.matrix.small) %in% genes.included,
                                      colnames(tom.matrix.small) %in% genes.included]
  #Generate the cluster to further verify the simmilarities
  ready.for.cluster = tom.matrix.small
  colnames(ready.for.cluster) = fromAny2GeneName(colnames(tom.matrix.small))
  rownames(ready.for.cluster) = colnames(ready.for.cluster)


  edge.data <- list()
  node.data <- list()

  n.edges <- 1
  n.nodes <- 1
  alt.names <- fromAny2GeneName(rownames(tom.matrix.small))
  nodes.stored <- vector(mode="logical",length=nrow(tom.matrix.small))

  for(i in 1:nrow(tom.matrix.small)){
    for(j in i:ncol(tom.matrix.small)){
      if(i != j){
        gene.i = rownames(tom.matrix.small)[i]
        gene.j = rownames(tom.matrix.small)[j]
        #Depending on the mode
        if(include.all.edges | gene.i %in% apriori.genes | gene.j %in% apriori.genes){

          edge.data[[n.edges]] <- list(fromNode=alt.names[i],toNode=alt.names[j],
                                       weight=tom.matrix.small[i,j],
                                       direction="undirected",
                                       fromAltName=rownames(tom.matrix.small)[i],
                                       toAltName=colnames(tom.matrix.small)[j])
          n.edges <- n.edges + 1
          if(!nodes.stored[i]){
            gene = rownames(tom.matrix.small)[i]
            if(gene %in% disease.genes)
              type = "disease"
            else if(gene %in% seed.genes)
              type = "seed"
            else
              type = "context"

            node.data[[n.nodes]] <- list(nodeName=alt.names[i],
                                         altName=rownames(tom.matrix.small)[i],
                                         biotype=biotypes[names(biotypes) == rownames(tom.matrix.small)[i]],
                                         module=net$moduleColors[names(net$moduleColors) %in% rownames(tom.matrix.small)[i]],
                                         adjacency=adjacencies[names(adjacencies) %in%
                                                                 rownames(tom.matrix.small)[i]],
                                         mm=mms[names(mms) %in% rownames(tom.matrix.small)[i]],
                                         type = type)
            nodes.stored[i] <- TRUE
            n.nodes <- n.nodes + 1
          }

          if(!nodes.stored[j]){
            gene = rownames(tom.matrix.small)[j]
            if(gene %in% disease.genes)
              type = "disease"
            else if(gene %in% seed.genes)
              type = "seed"
            else
              type = "context"
            gene.biotype = biotypes[names(biotypes) %in% colnames(tom.matrix.small)[j]]
            node.data[[n.nodes]] <- list(nodeName=alt.names[j],
                                         altName=rownames(tom.matrix.small)[j],
                                         biotype=gene.biotype,
                                         module=net$moduleColors[names(net$moduleColors) %in% colnames(tom.matrix.small)[j]],
                                         adjacency=adjacencies[names(adjacencies) %in%
                                                                 colnames(tom.matrix.small)[j]],
                                         mm=mms[names(mms) %in% rownames(tom.matrix.small)[j]],
                                         type = type)

            nodes.stored[j] <- TRUE
            n.nodes <- n.nodes + 1
          }
        }
      }

    }
  }

  #Now we have all nodes/edges, lets save them
  edge.matrix <- matrix(nrow=length(edge.data),ncol=length(edge.cols))
  colnames(edge.matrix) <- names(edge.data[[1]])
  for(i in 1:nrow(edge.matrix)){
    edge.matrix[i,] <- unlist(edge.data[[i]])
  }

  node.matrix <- matrix(nrow=length(node.data),ncol=length(node.cols))
  colnames(node.matrix) <- names(node.data[[1]])
  for(i in 1:nrow(node.matrix)){
    tryCatch(node.matrix[i,] <- unlist(node.data[[i]]),
             error = function(e)
             {
               print(paste0("Processing we get error ",
                            e$message))
             })
  }
  if(permutations > 0){
    estimates = cbind(names(estimates),estimates)
    colnames(estimates) = c("nodeName","p.value")
    write.table(estimates,paste0(folder,"estimates.txt"),
                quote=FALSE,row.names=F,col.names=T)
  }

  edges.file = paste0(folder,"edges.txt")
  nodes.file = paste0(folder,"nodes.txt")
  write.table(edge.matrix,
              edges.file,
              quote=FALSE,row.names=FALSE)
  write.table(node.matrix,
              nodes.file,
              quote=FALSE,row.names=FALSE)

  tab.to.return = NULL
  names(knn.data) = fromAny2GeneName(names(knn.data))
  for(gene in names(knn.data)){
    tab.to.return = cbind(tab.to.return,as.vector(fromAny2GeneName(knn.data[[gene]])))
  }
  colnames(tab.to.return) = names(knn.data)
  write.csv(tab.to.return,paste0(folder,"knn.csv"),
            quote=FALSE,row.names=FALSE)

  return(list(edges=edges.file,nodes=nodes.file))

}

createTOM = function(expr.data.file,
                     beta=11,
                     save.as=NULL,
                     net.type="signed",
                     debug=F){

  stopifnot(beta > 0 & beta < 40)
  stopifnot(net.type == "signed" | net.type == "unsigned")

  if(typeof(expr.data.file) == "character"){
    print(paste0("Creating matrix ",save.as," from expression data ",expr.data.file))
    expr.data <- readRDS(expr.data.file)
  }else{
    expr.data <- expr.data.file
  }
  cat("Creating TOM for",ncol(expr.data),"genes and",nrow(expr.data),"samples, beta",
      beta,"and type",net.type,"\n")
  if(debug)
    expr.data = expr.data[,1:1000]
  adjacency = adjacency(expr.data, power = beta, type = net.type )
  print("Adjacency matrix created")
  # Topological Overlap Matrix (TOM)
  # Turn adjacency into topological overlap
  print("Creating TOM")
  TOM = TOMsimilarity(adjacency)
  colnames(TOM) = colnames(expr.data)
  rownames(TOM) = colnames(TOM)

  if(!is.null(save.as)){
    cat("Saving TOM at",save.as,"\n")
    saveRDS(TOM,save.as)
  }
  else return(TOM)
}

getModuleTOM = function(tissue,
                        out.path,
                        module,
                        only.file=F,
                        which.one="CoExpROSMAP",
                        package=which.one,
                        instfolder="extdata"){
  if(is.null(out.path))
    out.path = system.file("", instfolder,
                           package = package)

  netf = getNetworkFromTissue(tissue=tissue,
                              which.one=which.one,
                              only.file = T)

  tfile = paste0(out.path,"/",netf,".",module,".tom.rds")

  if(only.file){
    return(tfile)
  }

  if(!file.exists(tfile))
    return(NULL)
  return(readRDS(tfile))
}

getModuleTOMs = function(tissue,beta,out.path,
                         which.one="CoExpROSMAP",
                         package=which.one,
                         instfolder="extdata"){

  if(is.null(out.path))
    out.path = system.file("", instfolder,
                           package = package)

  netf = getNetworkFromTissue(tissue=tissue,
                              which.one=which.one,
                              only.file = T)

  expr.data = getExprDataFromTissue(tissue=tissue,
                                    which.one=which.one)
  if(is.null(expr.data))
    return(expr.data)

  net = getNetworkFromTissue(tissue=tissue,
                             which.one=which.one)

  #Debug
  #rmask = sample(1:length(net$moduleColors),2000)
  #expr.data = expr.data[,rmask]
  #net$moduleColors = net$moduleColors[rmask]

  netf = basename(netf)

  tom = createTOM(expr.data.file = expr.data,
                  beta=beta)

  if(is.null(out.path))
    out.path = system.file("", instfolder,
                           package = package)

  modules = unique(net$moduleColors)
  lapply(modules,function(x){
    mask = net$moduleColors == x
    saveRDS(tom[mask,mask],paste0(out.path,"/",netf,".",x,".tom.rds"))
  })
}
