preservationOneWay <- function(network,
                               expr.data.files=NULL,
                               tissues=c("snig","putm"),
                               permutations=200,
                               maxModuleSize=5000,
                               maxGoldModuleSize=400,
                               randomSeed=1){

  cat("Entering preservation\n")
  n1.shortname = tissues[1]
  n2.shortname = tissues[2]

  if(typeof(network) == "character"){
    print(paste0("Reading network ",network))
    network1 <- readRDS(network)
  }else{
    network1 = network
  }


  print( tissue1 <- tissues[1] )
  print( tissue2 <- tissues[2] )
  print(paste0(tissue1," vs. ", tissue2))

  if( tissues[1] == tissues[2] ) stop("Can't do a preservation against self")
  options(stringsAsFactors = FALSE)

  print(paste0("Reading expression data for tissue ",tissues[1]))
  if(typeof(expr.data.files[[1]]) == "character"){
    expression.data1 <- readRDS(expr.data.files[[1]])
    cat(expr.data.files[[1]],"\n")

  }else{
    expression.data1 = expr.data.files[[1]]
  }
  print(expression.data1[1:5,1:5])

  print(paste0("Reading expression data for tissue ",tissues[2]))
  if(typeof(expr.data.files[[2]]) == "character"){
    expression.data2 <- readRDS(expr.data.files[[2]])
    cat(expr.data.files[[2]],"\n")
  }else
    expression.data2 = expr.data.files[[2]]
  print(expression.data2[1:5,1:5])


  ## Prepare the data

  #First we check for good genes and samples
  intersect.g = intersect(colnames(expression.data1),colnames(expression.data2))
  expression.data1 = expression.data1[,match(intersect.g,colnames(expression.data1))]
  expression.data2 = expression.data2[,match(intersect.g,colnames(expression.data2))]

  network1 = network1$moduleColors[match(intersect.g,names(network1$moduleColors))]
  network2 = network1

  cat("We'll use",ncol(expression.data1),"genes for the preservation analysis\n")

  multiExpr <- list()
  multiExpr [[1]] <- list(data = expression.data1)
  multiExpr [[2]] <- list(data = expression.data2)

  names(multiExpr) <- c(tissues[1], tissues[2])

  checkSets(multiExpr, checkStructure = FALSE, useSets = NULL)
  multiColor <- list( network1) #, network2 )
  names(multiColor) <- c(tissues[1]) #, tissues[2])
  ## Run the preservation statistics and save
  enableWGCNAThreads()
  print( WGCNAnThreads() )
  system.time( {
    mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                             nPermutations = permutations,
                             networkType="signed",
                             maxGoldModuleSize=maxGoldModuleSize,
                             ## no. of permutations
                             randomSeed=randomSeed,
                             verbose=3,
                             maxModuleSize=maxModuleSize)
  })

  getPreservationStatisticsOneWay(network=network,
                                  tissues=tissues,
                                  presRes=mp)

}


getPreservationStatisticsOneWay <- function(tissues,
                                            presRes){

  Zsummary <- NULL
  MedianRank <- NULL

  tissue1 = paste0("ref.",tissues[1])
  tissue2 = paste0("inColumnsAlsoPresentIn.",tissues[2])
  Z.tmp.list <- list(NULL)
  MR.tmp.list <- list(NULL)
  cat(tissue1, "\t", tissue2, "\n")
  fn	<- paste(tissue1, "vs", tissue2,sep="")
  mp = presRes
  return(list(Z=as.data.frame(mp$preservation$Z[[tissue1]][[tissue2]],stringsAsFactors=F),
              logp=as.data.frame(mp$preservation$log.pBonf[[tissue1]][[tissue2]],stringsAsFactors=F)))

}

getZsummaryPress = function(tissues,presRes,module,statistic="Zsummary.pres"){
  pData = getPreservationStatisticsOneWay(tissues=tissues,presRes=presRes)$Z
  mask = rownames(pData) == module
  return(pData[mask,statistic])
}

getMeanZsummaryPress = function(tissue,tissues,
                                package="CoExp10UKBEC",
                                folder="micro19K",
                                module,
                                statistic="Zsummary.pres"){
  path = system.file(package = package)
  if(path == ""){
    cat("Can´t find package ",package,"\n")
    return(NULL)

  }
  path = paste0(path,"/",folder)
  means = 0
  for(tother in tissues){
    #THAL.mic.net.19K.rds_vs_WHMT.mic.net.19K.rds_preserv.rds
    if(tissue < tother)
      fpres = paste0(path,"/",tissue,".mic.net.19K.rds_vs_",tother,".mic.net.19K.rds_preserv.rds")
    else
      fpres = paste0(path,"/",tother,".mic.net.19K.rds_vs_",tissue,".mic.net.19K.rds_preserv.rds")

    tlist = c(tissue,tother)
    cat("Reading from",fpres,"\n")
    partial = getZsummaryPress(tissues=tlist,
                               presRes=readRDS(fpres),
                               module=module,
                               statistic=statistic)
    print(partial)
    means = means + partial
  }
  return(means/length(tissues))
}

