

addNetworkToDDBB = function(netf,folder,
                            which.one,
                            tissue,
                            rewrite=F,
                            filter=c("GO","KEGG","REAC","HP"),
                            ensembl=TRUE,
                            do.go=T,
                            do.ct=T,
                            exclude.iea=T,
                            correction.method="gSCS",
                            organism = "hsapiens",
                            out.file=paste0(netf,"_gprofiler.csv")){

  if(do.go)
    go = getGProfilerOnNet(net.file=netf,
                                filter=filter,
                                ensembl=ensembl,
                                exclude.iea=exclude.iea,
                                correction.method=correction.method,
                           organism=organism,
                                out.file=out.file)
  else
    out.file = ""

  if(do.ct){
    ct = annotateByCellType(tissue=tissue,
                            which.one="new",
                            net.in=netf,
                            legend=tissue,
                            threshold=20,
                            plot.file=NULL,
                            return.processed=F)
    ctfile = paste0(netf,"_celltype.csv")
    write.csv(ct,ctfile)

  }else
    ctfile = ""

  addNet(which.one,
         tissue,
         netf,
         ctfile,
         out.file,"")



}



getNetworkCategories = function(){
  if(!is.null(coexp.nets))
    return(coexp.nets$which.one)
  return(NULL)
}

findNet = function(which.one,tissue){
  return(which(coexp.nets$tissue == tissue & coexp.nets$which.one == which.one))
}

saveDDBB = function(fileout){
  if(exists("coexp.nets")){
    if(nrow(coexp.nets) > 0){
      cat("Saving",nrow(coexp.nets),"network references in",fileout)
      write.table(coexp.nets,fileout,quote=F,sep="\t",col.names=T,row.names=F)
    }
  }
}

loadDDBB = function(filein,outtmp="/tmp/tempddbb.txt"){
  if(exists("coexp.nets"))
    saveDDBB(outtemp)
  coexp.nets <<- read.delim(filein,header=T,sep="\t")
  cat("Loading",nrow(coexp.nets),"from",filein,"\n")
}

addNet = function(which.one,tissue,netfile,ctfile,gofile,exprdatafile,overwrite){

  if(!exists("coexp.nets"))
    initDb(mode="empty")

  if(sum(coexp.nets$tissue == tissue & coexp.nets$which.one == which.one) == 0){
    cat("Adding new network",tissue,"to the category",which.one,"to the database\n")
    coexp.nets[nrow(coexp.nets) + 1,] <<- c(which.one,
                                            tissue,
                                            netfile,ctfile,gofile,exprdatafile)
  }else{
    if(overwrite){
      mask = coexp.nets$tissue == tissue & coexp.nets$which.one == which.one
      coexp.nets[mask,] <<- c(which.one,
                                              tissue,
                                              netfile,ctfile,gofile,exprdatafile)
    }

    cat("Network",tissue,"in category",which.one,"already exists, we can´t overwrite\n")
  }
}



getAvailableNetworks = function(category="rosmap"){
  if(!is.null(coexp.nets))
    return(coexp.nets$tissue[coexp.nets$which.one == category])

  return(NULL)
}

initDbGTEX = function(mandatory=F){
  if(!dir.exists("supplementary/rdsnets/gtexv6/") & mandatory)
    stop("GTEx networks are not available\n")

  if(!dir.exists("supplementary/rdsnets/gtexv6/")){
    cat("GTEx networks won´t be available, you have to install them\n")
    return
  }

  nets = getGTExTissues()
  for(net in nets){
    addNet("gtexv6",net,getGTExNet(net),
                                            paste0(getGTExNet(net),"_celltype.csv"),
                                            paste0(getGTExNet(net),"_gprof.csv"),
                                            paste0("supplementary/rdsnets/gtexv6/",
                                                   getGTExNet(net),".resids.rds"))
  }
}

initExonic = function(mandatory=F){

  if(!dir.exists("supplementary/rdsnets/exonic/") & mandatory)
    stop("RytenLab putamen and Snigra networks are not available\n")

  if(!dir.exists("supplementary/rdsnets/exonic/")){
    cat("RytenLab putamen and Snigra won´t be available, you have to install them\n")
    return
  }
  addNet("exonic","SNIG","supplementary/rdsnets/exonic/netSNIG.7.12.it.30.rds",
         "supplementary/rdsnets/exonic/netSNIG.7.12.it.30.rds.cell.type.csv",
         "supplementary/rdsnets/exonic/netSNIG.7.12.it.30.rds_gprofiler.csv",
         "supplementary/rdsnets/exonic/resids.SNIG.7.rds")
  addNet("exonic","PUTM","supplementary/rdsnets/exonic/netPUTM.8.13.it.30.rds",
         "supplementary/rdsnets/exonic/netPUTM.8.13.it.30.rds.cell.type.csv",
         "supplementary/rdsnets/exonic/netPUTM.8.13.it.30.rds_gprofiler.csv",
         "supplementary/rdsnets/exonic/resids.PUTM.8.rds")

}

initUKBECMicroarray = function(mandatory=F){

  if(!dir.exists("supplementary/rdsnets/micro19K/") & mandatory)
    stop("RytenLab 10 regions microarray based networks are not available\n")

  if(!dir.exists("supplementary/rdsnets/micro19K/")){
    cat("RytenLab 10 regions microarray based networks won´t be available, you have to install them\n")
    return
  }

  #Init microarray nets
  nets = getMicTissues()
  coexp.nets[["micro19K"]] <<- NULL
  coexp.data[["micro19K"]] <<- NULL

  for(net in nets){
    addNet("micro19K",net,paste0("supplementary/rdsnets/micro19K/net",net,
                                          ".12.signed.it.20.rds_cor_pca.rds"),
           paste0("supplementary/rdsnets/micro19K/net",net,
                  ".12.signed.it.20.rds",".USER_terms.csv"),
           paste0("supplementary/rdsnets/micro19K/net",net,
                  ".12.signed.it.20.rds_cor_pca.rds_gprofiler.csv"),
           paste0("supplementary/rdsnets/micro19K/",net,
                                          ".mic.expr.data.19K.rds"))
  }
}


initROSMAP = function(mandatory=F){

  if(!dir.exists("rdsnets/rosmap/") & mandatory)
    stop("ROSMAP based networks are not available\n")

  if(!dir.exists("rdsnets/rosmap/")){
    cat("ROSMAP based networks won´t be available, you have to install them\n")
    return
  }
  addNet("rosmap","cogdxad","rdsnets/rosmap/netad.8.it.50.rds",
         "rdsnets/rosmap/netad.8.it.50.rds.celltype.csv",
         "rdsnets/rosmap/netad.8.it.50.rds_gprofiler.csv",
         "rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.ad.rds")
  addNet("rosmap","cogdxprobad","rdsnets/rosmap/netprobad.11.it.50.rds",
         "rdsnets/rosmap/netprobad.11.it.50.rds.celltype.csv",
         "rdsnets/rosmap/netprobad.11.it.50.rds_gprofiler.csv",
         "rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.probad.rds")
  addNet("rosmap","cogdxnotad","rdsnets/rosmap/netnotad.8.it.50.rds",
         "rdsnets/rosmap/netnotad.8.it.50.rds.celltype.csv",
         "rdsnets/rosmap/netnotad.11.it.50.rds_gprofiler.csv",
         "rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.not.rds")
  addNet("rosmap","all","supplementary/rdsnets/rosmap/netROSMAPSingle.6.it.50.rds",
         "supplementary/rdsnets/rosmap/netROSMAPSingle.6.it.50.rds.celltype.csv",
         "supplementary/rdsnets/rosmap/netROSMAPSingle.6.it.50.rds_gprofiler.csv",
         "supplementary/rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.rds")
}

initDb = function(mode = "FULL",...){

  coexp.nets  <<- as.data.frame(list(which.one=character(0),
                                     tissue=character(0),
                                     net=character(0),
                                     ctype=character(0),
                                     go=character(0),
                                     data=character(0)),
                                stringsAsFactors=F)

  if(mode == "GTEX"){
    initDbGTEX(...)
    return
  }

  if(mode == "EXONIC"){
    initExonic(...)
    return
  }

  if(mode == "MICROARRAY"){
    initUKBECMicroarray(...)
    return
  }

  #Init gtexv6 nets
  if(mode == "FULL"){

    initDbGTEX(...)
    #Init exonic RNA-seq nets
    initExonic(...)
    initUKBECMicroarray(...)


  }
  cat("These are the networks available now\n")
  print(coexp.nets)
}

getCovariates = function(tissue,which.one,cov){
  if(which.one == "rosmap"){
    expr.data = getExprDataFromTissue(tissue=tissue,which.one=which.one)
    key = read.csv(paste0("rdsnets/rosmap/ROSMAP_IDkey.csv"))
    covs = read.csv(paste0("rdsnets/rosmap/ROSMAP_clinical.csv"))
    ids = rownames(expr.data)
    mask = key$projid[match(ids,key$mrna_id)]
    nonmatchingids = ids[is.na(mask)]
    goodids = NULL
    for(id in nonmatchingids){
      subids = str_split(id,"_")
      recid = paste0(subids[[1]][1],"_",subids[[1]][2])
      goodids = c(goodids,recid)
    }
    mask[is.na(mask)] = key$projid[match(goodids,key$mrna_id)]
    samples = mask
    #samples = rosmap.fromRNAseqID2ProjectID(rownames(expr.data))
    gender = as.factor(covs$msex[match(samples,covs$projid)])
    pmi = as.numeric(covs$pmi[match(samples,covs$projid)])
    braaksc = as.factor(covs$braaksc[match(samples,covs$projid)])
    cogdx = as.factor(covs$cogdx[match(samples,covs$projid)])
    educ = as.numeric(covs$educ[match(samples,covs$projid)])
    ceradsc = as.factor(covs$ceradsc[match(samples,covs$projid)])
    age = as.character(covs$age_death[match(samples,covs$projid)])

    #Impute
    pmi[is.na(pmi)] = mean(pmi[!is.na(pmi)])
    age[grep("90\\+",age)] = "90"
    age = as.numeric(age)
    race = as.factor(covs$race[match(samples,covs$projid)])

    batch = str_split(rownames(expr.data),"_")
    batch = as.factor(unlist(lapply(batch,function(x){return(x[[3]])})))

    toreturn = data.frame(batch,gender,pmi,age,race,braaksc,cogdx,educ,ceradsc)

    rownames(toreturn) = rownames(expr.data)
    return(toreturn)
  }else if(which.one == "gtexv6"){
    expr.data = getExprDataFromTissue(tissue=tissue,which.one=which.one)
    covs = read.table(paste0("supplementary/rdsnets/gtexv6/gtexv6covariates.txt"),
                      header=T,
                      sep="\t")
    ids = rownames(expr.data)
    covs$ID = gsub("-","\\.",covs$ID)
    return(covs[match(ids,covs$ID),])
  }else if(which.one == "exonic"){
    expr.data = getExprDataFromTissue(tissue=tissue,which.one=which.one)
    ids = rownames(expr.data)
    covs = read.csv(paste0("supplementary/rdsnets/exonic/all_samples_data.csv"),stringsAsFactors=F)
    covs = covs[match(ids,covs$A.CEL_file),c(2:6,8,9)]
    colnames(covs) = c("Age","PMI","RIN","Gender","tissue","causeofdeath","batch")
    return(covs)
  }else if(which.one == "micro19K"){
    expr.data = getExprDataFromTissue(tissue=tissue,which.one=which.one)
    ids = rownames(expr.data)
    covs = read.csv(paste0("supplementary/rdsnets/exonic/all_samples_data.csv"),stringsAsFactors=F)
    covs$U.SD_No = gsub("/","_",covs$U.SD_No)
    covs = covs[match(ids,covs$U.SD_No),c(2:6,8,9)]
    colnames(covs) = c("Age","PMI","RIN","Gender","tissue","causeofdeath","batch")
    return(covs)
  }
  return(NULL)
}


getMicTissues = function(){
  #	load("/Users/juanbot/Dropbox/kcl/WGCNA/data_no_outliers.RData")
  #	tissues = names(dat.clean)
  #	rm(dat.clean)
  return(c("CRBL", "FCTX", "HIPP", "MEDU", "OCTX", "PUTM", "SNIG", "TCTX",
           "THAL", "WHMT"))
}

getGTExNet = function(tissue){

  the.dir = paste0("supplementary/rdsnets/gtexv6/")
  files = list.files(the.dir)

  net.file = files[grep(paste0("net",tissue,"\\.\\d+\\.it\\.50\\.rds$"),files)]
  if(length(net.file) == 0)
    return(NULL)

  return(paste0(the.dir,net.file))
}

getGTExTissues = function(){

  the.dir = paste0("supplementary/rdsnets/gtexv6/")
  files = list.files(the.dir)

  net.files = files[grep("net\\w+\\.\\d+\\.it\\.50\\.rds$",files)]
  net.files = gsub("net","",net.files)
  net.files = gsub(".\\d+\\.it\\.50\\.rds$","",net.files)

  return(net.files)
}


getInfoFromTissue = function(which.one,tissue,what,...){
  if(what == "net")
    return(getNetworkFromTissue(which.one=which.one,tissue=tissue))
  if(what == "expr")
    return(getExprDataFromTissue(which.one=which.one,tissue=tissue))
  if(what == "go")
    return(getGOFromTissue(which.one=which.one,tissue=tissue))
  if(what == "ctype")
    return(getCellTypeFromTissue(which.one=which.one,tissue=tissue))

}

#' getNetworkFromTissue - Get predefined network or create yours
#'
#' @param tissue
#' @param which.one
#' @param only.file
#' @param genes
#' @param modules
#'
#' @return
#' @export
#'
#' @examples
getNetworkFromTissue = function(tissue="SNIG",which.one="exonic",only.file=F,genes=NULL,modules){

  if(which.one == "new"){
    cat("New network, reading from",tissue,"\n")
    return(readRDS(tissue))

  }

  if(which.one == "onthefly"){
    stopifnot(!is.null(genes))
    net = NULL
    net$moduleColors = modules
    names(net$moduleColors) = genes
    net$MEs = NULL
    return(net)
  }

  file = coexp.nets$net[findNet(which.one,tissue)]
  if(only.file)
    return(file)
  else{
    if(file.exists(file))
      return(readRDS(file))
  }
  return(NULL)

  cat(paste0("When getting the network for tissue ",tissue," we don't know ",which.one,"\nReturning NULL\n"))
  return(NULL)
}

getNetworkEigengenes = function(tissue,which.one){
  mes = getNetworkFromTissue(tissue=tissue,which.one=which.one)$MEs
  class(mes) = c("data.frame","meigengenes")
  return(mes)
}

getModulesMostRelevantGenes = function(tissues="SNIG",
                                       modules="red",which.ones,
                                       n=10,cutoff=-1){
  toreturn = NULL
  for(tissue in tissues){
    i = which(tissues == tissue)
    expr.data.file = getExprDataFromTissue(tissue=tissue,
                                           which.one=which.ones[i],only.file=T)
    mm = getMM(which.one=which.ones[i],tissue=tissue,
               table.format = T,genes=NULL,
               expr.data.file=expr.data.file)

    mm = mm[mm$module == modules[i],]
    mm = mm[order(mm$mm,decreasing=TRUE),]
    ncount = 0
    if(n > 0){
      mask = 1:n
      ncount = n
    }
    else{
      mask = mm$mm > cutoff
      ncount = sum(mask)
    }
    toreturn = rbind(toreturn,c(rep(tissue,ncount),rep(modules[i],ncount),rep(which.ones[i],ncount),
                                names(mm$mm[mask]),mm$mm[mask]))

  }
  toreturn
}

getModuleMostRelevantGenes = function(tissue="SNIG",
                                      module="red",which.one,n=10,
                                      cutoff=-1,
                                      expr.data.file=NULL){

  mm = getMM(which.one=which.one,tissue=tissue,table.format = T,
             genes=NULL,
             expr.data.file=expr.data.file)

  mm = mm[mm$module == module,]
  mm = mm[order(mm$mm,decreasing=TRUE),]
  if(n > 0)
    return(mm[1:n,])
  return(mm[mm$mm > cutoff,])
}

getExprDataFromTissue = function(tissue="SNIG",which.one="rnaseq",only.file=F){

  if(which.one == "new"){
    cat("Reading new expression data from",tissue,"\n")
    if(only.file)
      return(tissue)
    return(readRDS(tissue))
  }

  file = coexp.nets$data[findNet(which.one,tissue)]
  if(length(file)){
    if(only.file)
      return(file)
    else{
      if(file.exists(file)){
        if(which.one == "rosmap"){
          expr.data = data.frame(readRDS(file))
          colnames(expr.data) = gsub("\\.[0-9]+","",colnames(expr.data))
          return(expr.data)
        }
        return(readRDS(file))
      }

    }
    return(NULL)

  }

  cat(paste0("When getting the expression data for tissue ",tissue," we don't know ",which.one,"\nReturning NULL\n"))
  return(NULL)
}

getModulesFromTissue = function(tissue="SNIG",which.one="rnaseq",in.cluster=F){
  return(unique(getNetworkFromTissue(tissue,which.one)$moduleColors))
}

getGenesFromModule = function(tissue="SNIG",which.one="rnaseq",module="black"){
  net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  if(is.null(module)) return(names(net$moduleColors))
  return(names(net$moduleColors)[net$moduleColors == module])
}

getModuleFromGenes = function(tissue="SNIG",which.one="exonic",genes){
  net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  return(net$moduleColors[match(genes,names(net$moduleColors))])
}

getGOFromTissues = function(tissues,which.ones,modules){
  n = length(modules)
  toreturn = NULL
  for(i in 1:n){
    cat("Working with",which.ones[i],", ",tissues[i]," and ",modules[i],"\n")
    localtr = getGOFromTissue(tissue=tissues[i],
                              which.one=which.ones[i],
                              module=modules[i])
    key = paste0(which.ones[i],"_",tissues[i],"_",modules[i])
    toreturn = rbind(toreturn,cbind(rep(key,nrow(localtr)),localtr))
  }
  colnames(toreturn)[1] = "key"
  return(toreturn)
}


getGOFromTissue = function(tissue="SNIG",which.one="rnaseq",module=NULL){
  gofile = coexp.nets$go[findNet(which.one,tissue)]
  if(length(gofile)){
      go = read.csv(gofile,stringsAsFactors=F)
      if(!is.null(module))
        return(go[go$query.number == module,])
      return(go)
  }

  cat(paste0("When getting the network GO for tissue ",tissue," we don't know ",
             which.one,"\nReturning NULL\n"))
  return(NULL)
}

getCellTypeFromTissue = function(tissue="SNIG",which.one="rnaseq",module=NULL){

  ctfile = coexp.nets$ctype[findNet(which.one,tissue)]
  if(length(ctfile)){
      ct = data.frame(read.csv(ctfile,stringsAsFactors=F))
      if(!is.null(module)){
        if(which.one == "micro19K"){
          ct = ct[ct$InputCategories == module,]
          ctvec = ct[,4]
          names(ctvec) = ct[,2]

        }else{
          ctvec = ct[,module]
          names(ctvec) = ct[,1]

        }
        return(ctvec)
      }

      else return(ct)
  }


  cat(paste0("When getting the network cell type signals for tissue ",tissue,
             " we don't know ",which.one,"\nReturning NULL\n"))
  return(NULL)
}

