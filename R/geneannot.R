

reportExists = function(experiment,tissue,module){
  if(exists("coexp.report.go.all"))
    return(!is.null(coexp.report.go.all[[experiment]][[tissue]][[module]]))
  return(FALSE)
}

reportGet = function(experiment,tissue,module){
  if(reportExists(experiment,tissue,module))
    return(coexp.report.go.all[[experiment]][[tissue]][[module]])
  return(NULL)
}

reportAdd = function(experiment,tissue,module,report){
  if(!exists("coexp.report.go.all")){
    coexp.report.go.all <<- NULL
    coexp.report.go.all[[experiment]] <<- NULL
    coexp.report.go.all[[experiment]][[tissue]] <<- NULL
    coexp.report.go.all[[experiment]][[tissue]][[module]] <<- report
  }else{
    if(is.null(coexp.report.go.all)){
      coexp.report.go.all[[experiment]] <<- NULL
      coexp.report.go.all[[experiment]][[tissue]] <<- NULL
      coexp.report.go.all[[experiment]][[tissue]][[module]] <<- report
    }else if(is.null(coexp.report.go.all[[experiment]])){
      coexp.report.go.all[[experiment]][[tissue]] <<- NULL
      coexp.report.go.all[[experiment]][[tissue]][[module]] <<- report
    }else #if(is.null(report.go.all[[experiment]][[tissue]])){
      coexp.report.go.all[[experiment]][[tissue]][[module]] <<- report
  }
}

reportOnFETGenes = function(tissue,genes,which.one,net){

  #Everything will be turned into Ensembl
  if(is.null(net))
    net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  names(net$moduleColors) =  fromAny2Ensembl(names(net$moduleColors))
  en.genes = fromAny2Ensembl(genes)

  #Modules within genes
  mods = net$moduleColors[match(en.genes,names(net$moduleColors))]
  mods[is.na(mods)] = "void"
  final.report = NULL
  mult.mod.genes = NULL

  p.val.mods = rep(1,length(mods))
  table.mods = table(mods)
  table.mods = table.mods[names(table.mods) != "void"]
  total.specific = sum(mods != "void")
  total.net = length(net$moduleColors)
  for(module in names(table.mods)){
    n.module = sum(net$moduleColors == module)
    n.module.and.specific = table.mods[[module]]
    p.val = testGeneSet(n.module,n.module.and.specific,total.specific,total.net)
    p.val.mods[mods == module] = signif(p.val$p.value,4)
  }

  the.table = cbind(rep(tissue,length(genes)),rep(which.one,length(genes)),
                    genes,en.genes,mods,p.val.mods)
  colnames(the.table) = c("tissue","which.one","gene","ensgene","module","p.val.mods")
  rownames(the.table) = genes
  return(the.table)
}

globalReportOnGenes = function(tissues,
                               categories,
                               genes){
  stopifnot(length(tissues) > 0)
  stopifnot(length(tissues) == length(categories))

  allreport = NULL
  singcats = unique(categories)
  for(category in singcats){
    cat("Working on",category,"\n")
    mask = categories == category
    localtissues = tissues[mask]
    cat("Working on",paste0(localtissues,collapse=", "),"\n")
    report = reportOnGenesMultipleTissue(tissues=localtissues,
                                         genes=genes,
                                         which.one=category)$report

    allreport = rbind(allreport,report)
  }

  #Get the bonferroni factor
  ntests = 0
  pvals = NULL
  cats = tiss = mods = NULL
  for(category in singcats){
    mask = categories == category
    localtissues = tissues[mask]
    for(ltissue in localtissues){
      localreport = allreport[allreport$category == category & allreport$tissue == ltissue,]
      modules = unique(localreport$module)
      for(module in modules){
        cats = c(cats,category)
        tiss = c(tiss,ltissue)
        mods = c(mods,module)
        pvals = c(pvals,unique(allreport[allreport$category == category &
                                    allreport$tissue == ltissue &
                                    allreport$module == module,]$fisher))
      }
    }
  }

  fdrpvalues = p.adjust(as.numeric(pvals),method = "fdr")
  bonfpvalues = p.adjust(as.numeric(pvals),method = "bonferroni")
  allreport = cbind(allreport,rep(0,nrow(allreport)),rep(0,nrow(allreport)))

  mask = (ncol(allreport) - 1):ncol(allreport)
  colnames(allreport)[mask] = c("FDR","Bonferroni")
  for(mytest in 1:length(cats)){
    indexes = which(allreport$category == cats[mytest] & allreport$tissue == tiss[mytest] &
                                allreport$module == mods[mytest])
    newdata = matrix(nrow=length(indexes),ncol=2,data=c(fdrpvalues[mytest],bonfpvalues[mytest]),byrow=T)
    allreport[indexes,mask] = newdata
  }
  return(allreport)

}



reportOnGenes = function(tissue,
                         genes,
                         silent=F,
                         which.one="signedrnaseq",
                         alt.probes=NULL,
                         ens.correction=NULL,
                         gwases=NULL){

  genes = unique(genes)
  #Everything will be turned into Ensembl
  net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  names(net$moduleColors) =  fromAny2Ensembl(names(net$moduleColors))
  expr.dataf = getExprDataFromTissue(tissue=tissue,which.one=which.one,only.file = T)
  expr.data = readRDS(expr.dataf)
  print(head(colnames(expr.data)))
  colnames(expr.data) = fromAny2Ensembl(colnames(expr.data))
  en.genes = fromAny2Ensembl(genes)
  gene.names = fromAny2GeneName(genes)
  print(en.genes)
  if(!is.null(ens.correction)){
    for(i in seq(from=1,to=length(ens.correction),by=2)){
      mask = which(en.genes == ens.correction[i])
      if(sum(mask) > 0){
        en.genes[mask] = ens.correction[i+1]
        cat("Substituting", ens.correction[i],"by",ens.correction[i+1],"\n")
      }
    }
  }

  final.report = NULL
  mult.mod.genes = NULL
  mms = getMM(net=net,
              expr.data.file=expr.data,
              tissue=tissue,
              genes=en.genes)
  modules = na.omit(unique(getModuleFromGenes(which.one=which.one,tissue=tissue,genes=en.genes)))

  for(i in 1:nrow(mms)){
    gene = mms[i,2]
    mod = mms[i,3]

    cat("Working on gene",gene,"\n")

    #mod = unique(mod)
    if(length(mod) >= 1){
      cat("Gene found in module",mod,"in ",tissue,"\n")
      report = reportOnModule(tissue,mod,
                                which.one=which.one)
        mm = mms[i,4]
        report = c(gene = gene,
                   ensgene = mms[i,1],
                   mm = signif(as.numeric(mm),4),
                   report)

        if(!silent){
          print(paste0("Gene ",gene, " is in ",tissue," in ",mod, " (MM ",
                       mm,")"))
          print(paste0("Report for module ",mod))
          print(report)
        }

        #Get the evidence for cell specific genes based on WGCNA userLists
        wgcna.cell.cats = getFriendlyNameForCategories(getCategoriesForGene(gene))
        if(length(wgcna.cell.cats) == 0){
          wgcna.cell.cats = "void"
        }

        if(sum(mms[,1] == gene) > 1){
          j = sum(mms[1:i,1] == gene)
          final.report[[paste0(gene,"_version_",j)]] = report

        }

        else
          final.report[[gene]] = report
    }else{
      cat("Gene not found in",tissue," gene set",which.one,"\n")
      mod = "void"
      mm = 0
      report = NULL
    }
  }
  #And now we add the fisher exact test results

  #But before that we convert to dataframe
  the.table = NULL
  genes.with.report = names(final.report[unlist(lapply(final.report,function(x){
    return(!is.null(x))}))])
  for(gene in genes.with.report){
    report = final.report[[gene]]
    report$go.report = paste0(report$go.report,collapse=", ")
    report$pd.genes = paste0(report$pd.genes,collapse=", ")
    report.t = as.data.frame(report,stringsAsFactors=F)
    the.table = rbind(the.table,report.t)
  }

  #Now we add the VEGAS results if any
  if(!is.null(gwases)){
    modules = the.table$module
    vegastoadd = NULL
    for(module in modules){
      addbymodule = NULL
      for(gwas in gwases){
        newsignal = gsa.getSignal(which.one=which.one,tissue=tissue,gwas=gwas,module=module)
        if(!is.null(newsignal)){
          pvalue = newsignal$pvalue
          addbymodule = cbind(addbymodule,pvalue)
          colnames(addbymodule)[ncol(addbymodule)] = paste0(gwas,"vegaspvalue")
        }

      }
      if(!is.null(addbymodule)){
        vegastoadd = rbind(vegastoadd,addbymodule)
      }

    }
    if(!is.null(vegastoadd))
      the.table = cbind(the.table,vegastoadd)
  }

  total.specific = length(genes.with.report)
  total.net = length(net$moduleColors)
  mods = the.table$module

  #It could be that there are more than 1 modules per gene
  mult.mod.index = grep(",",mods)
  for(single.index in mult.mod.index){
    #Get the first module
    few.modules = mods[single.index]
    few.modules = str_split(few.modules,", ")
    mods[single.index] = few.modules[[1]][[1]]
  }

  if(!is.null(mods)){
    p.val.mods = vector(mode="numeric",length=length(mods))
    names(p.val.mods) = mods
    table.mods = table(mods)
    table.mods = table.mods[names(table.mods) != "void"]
    for(module in names(table.mods)){
      n.module = sum(net$moduleColors == module)
      n.module.and.specific = table.mods[[module]]
      p.val = testGeneSet(n.module,n.module.and.specific,total.specific,total.net)
      p.val.mods[names(p.val.mods) == module] = signif(p.val$p.value,4)
    }
    the.table = cbind(the.table,p.val.mods)
    the.table = as.data.frame(the.table,stringsAsFactors=F)
    rownames(the.table) = names(final.report[genes.with.report])
    colnames(the.table)[ncol(the.table)] = "fisher"

  }
  cat("These genes are found in more than 1 module:",mult.mod.genes,"\n")
  return(list(report=the.table,mult.genes=mult.mod.genes))
}

reportOnGenesMultipleTissue = function(tissues,genes,silent=F,
                                       which.one="signedrnaseq",
                                       alt.probes=NULL,
                                       out.file=NULL,
                                       gwases=NULL){

  out.table = NULL

  problematic.genes = NULL
  for(tissue in tissues){

    tissue.report = reportOnGenes(tissue=tissue,
                                  genes=genes,
                                  silent=silent,
                                  which.one=which.one,
                                  alt.probes=alt.probes,
                                  gwases=gwases)

    if(!is.null(tissue.report$report)){
      problematic.genes[[tissue]] = tissue.report$mult.genes
      tissue.report = tissue.report$report
      tissue.report.ready = cbind(rep(tissue,nrow(tissue.report)),tissue.report)
      colnames(tissue.report.ready) = c("tissue",colnames(tissue.report))
      if(is.null(out.table)){
        out.table = tissue.report.ready

      }else{
        out.table = rbind(out.table,tissue.report.ready)
      }
    }
  }

  #We have to correct p.val.mods here!
  # if(mod.level.correct){
  #   modules = 0
  #   for(tissue in tissues)
  #     modules = modules + length(getModulesFromTissue(tissue=tissue,
  #                                                     which.one=which.one))
  #   bonf = out.table$fisher * modules
  #   bonf[bonf > 1] = 1
  #   out.table = cbind(out.table,bonf)
  #
  #   colnames(out.table)[ncol(out.table)] = "bonfpval"
  # }

  out.table = cbind(rep(which.one,nrow(out.table)),out.table)
  colnames(out.table)[1] = "category"


  if(!is.null(out.file))
    write.csv(out.table,out.file)
  return(list(report=out.table,mult.genes=problematic.genes))
}

functionalReportOnModule = function(tissue="SNIG",module,which.one="rnaseq"){

  gprof.file = findGO(which.one=which.one,tissue=tissue)

  go <- read.csv(gprof.file,stringsAsFactors=FALSE)
  go$X = NULL
  return(go[go$query.number == module,])
}

reportOnModule = function(tissue="SNIG",
                          module,
                          which.one="rnaseq",
                          how.many=5,
                          simple.cell.types=F,
                          cell.type.threshold=0.05,
                          ctcollapse=T){

  tmp.report = reportGet(experiment=which.one,tissue=tissue,module=module)
  if(!is.null(tmp.report))
    return(tmp.report)

  include.pd = T
  pd.genes = "void"
  if(include.pd){
    pd.genes = read.table(paste0(system.file("", "", package = "CoExpNets"),"/pd_genes.txt"),
                          stringsAsFactors=F,header=F)$V1
    genes.in.module = getGenesFromModule(tissue=tissue,which.one=which.one,module=module)
    stopifnot(!is.null(genes.in.module))
    if(substr(genes.in.module[1],1,4) == "ENSG")
      genes.in.module = fromEnsembl2GeneName(genes.in.module)
    mask = pd.genes %in% genes.in.module
    if(sum(mask) > 0)
      pd.genes = pd.genes[mask]
    else
      pd.genes = "void"
  }
  go = getGOFromTissue(which.one=which.one,tissue=tissue)
  cat("done GO\n")

  p.values <- go$p.value[go$query.number == module & (go$domain %in% "BP")]
  terms <- go$term.name[go$query.number == module & (go$domain %in% "BP")]
  term.ids = go$term.id[go$query.number == module & (go$domain %in% "BP")]

  if("IC" %in% colnames(go)){
    ics = go$IC[go$query.number == module & (go$domain %in% "BP")]
    p.values = p.values[!is.na(ics)]
    terms = terms[!is.na(ics)]
    ics = na.omit(ics)
  }

  else
    ics = NULL
  go.report = ""
  if(length(p.values) > 0){
    go.count = length(p.values)
    if(how.many > go.count)
      how.many = go.count
    the.order <- order(-log10(p.values),decreasing=TRUE)
    if(is.null(ics))
      go.report = paste0(terms[the.order][1:how.many]," ",term.ids[the.order][1:how.many],
                         " (p-value ",signif(p.values[the.order][1:how.many],4),")")
    else{
      p.values = p.values[the.order]
      terms = terms[the.order]
      term.ids = term.ids[the.order]
      ics = ics[the.order]

      p.values.f = p.values[ics > 2.5]
      terms.f = terms[ics > 2.5]
      term.ids.f = term.ids[ics > 2.5]

      if(length(p.values.f) > 0){
        if(length(p.values.f) < how.many)
          how.many = length(p.values.f)
        go.report = paste0(terms.f[1:how.many]," ",term.ids.f[1:how.many],
                           " (p-value ",signif(p.values.f[1:how.many],4),")")
      }else
        go.report = paste0(terms[1:how.many], " ",term.ids.f[1:how.many],
                           " (p-value ",signif(p.values[1:how.many],4),")")
    }
    go.report.size = go.count
    go.report.signif = sum(-log10(p.values))
  }else{
    go.report = "void"
    go.report.size = go.report.signif = 0
  }

  cat("Now cell type")
  #We include cell types here
  cell.types = getCellTypeFromTissue(which.one=which.one,tissue=tissue)


  cell.type.report = "void"
  if(any(colnames(cell.types) %in% module)){
    cell.types = cell.types[order(cell.types[,module],decreasing=F),]
    cell.type.row = rownames(cell.types)[which(cell.types[,module] <= cell.type.threshold)]
    if(length(cell.type.row) > 0){
      cell.type.report = NULL
      if(!simple.cell.types){

        for(i in 1:length(cell.type.row)){
          if(ctcollapse)
            cell.type.report = paste0(cell.type.report," ",
                                      getFriendlyNameForColumnCategories(cell.type.row[i]),
                                      " (p-value ",signif(cell.types[cell.type.row[i],module],4),"). ")
          else
            cell.type.report = c(cell.type.report,
                                 paste0(getFriendlyNameForColumnCategories(cell.type.row[i]),
                                        " (p-value ",signif(cell.types[cell.type.row[i],module],4),")"))

        }
      }else{
        cell.type.report = unique(getFriendlyNameForColumnCategories(names(cell.type.row)))
      }
    }
  }
  #preservation = getZSummary(tissue=tissue,module=module,which.one=which.one)
  preservation = NA
  if(is.na(preservation))
    preservation = "void"

  to.return = list(module=module,
                   size=length(getGenesFromModule(tissue=tissue,which.one=which.one,
                                                  module=module)),
                   go.report=go.report,
                   pd.genes=pd.genes,
                   preservation=preservation,
                   cell.type.pred=cell.type.report)
  reportAdd(which.one,tissue,module,to.return)
  return(to.return)


}


getZSummary = function(tissue="SNIG",
                       modules=c("red"),
                       which.one="signedrnaseq"){

  if(is.null(modules))
    modules = unique(getNetworkFromTissue(tissue=tissue,which.one=which.one)$moduleColors)

  source = 1
  target = 2
  if(which.one == "micro19K")
  {
    source = 2
    target = 1
    pr.file = paste0("supplementary/rdsnets/micro19K/PUTM.mic.net.19K.rds_vs_SNIG.mic.net.19K.rds_preserv.rds")
  }else if(which.one == "exonic"){
    pr.file = paste0("supplementary/rdsnets/exonic/SNIG_vs_PUTM_preserv.rds")
  }else if(which.one == "exonic.vs.gtex"){
    pr.file = paste0("supplementary/rdsnets/gtexv6/",
                     tissue,".UKBEC_vs_",tissue,".GTEx_preserv.rds")
    stats = readRDS(pr.file)
    loc.modules = rownames(stats$preservation$Z[[1]][[2]])
    return(stats$preservation$Z[[1]][[2]][match(modules,loc.modules),"Zsummary.pres"])
  }else if(which.one == "exonic.vs.micro"){
    pr.file = paste0("supplementary/rdsnets/exonic/",
                     tissue,"-RNAseq_vs_",tissue,"-microarray_preserv.rds")
    stats = readRDS(pr.file)
    loc.modules = rownames(stats$preservation$Z[[1]][[2]])
    return(stats$preservation$Z[[1]][[2]][match(modules,loc.modules),"Zsummary.pres"])
  }else if(which.one == "gtexv6"){
    return(NA)
  }else if(which.one == "rosmap"){
    return(NA)
  }else{
    stop(paste0("In getZSummary which.one=",which.one," not known"))
  }

  stats = readRDS(pr.file)

  if(tissue == "SNIG"){
    loc.modules = rownames(stats$preservation$Z[[source]][[target]])
    return(stats$preservation$Z[[source]][[target]][match(modules,loc.modules),"Zsummary.pres"])
  }else if (tissue == "PUTM"){
    loc.modules = rownames(stats$preservation$Z[[target]][[source]])
    return(stats$preservation$Z[[target]][[source]][match(modules,loc.modules),"Zsummary.pres"])
  }else return(NA)
}


getZsummary10MicTissues1Module = function(tissue,module){
  tissues = getMicTissues()
  i = which(tissues == tissue)
  if(i == 1)
    lefttissues = NULL
  else
    lefttissues = tissues[1:(i - 1)]
  if(i == 10)
    righttissues = NULL
  else
    righttissues = tissues[(i+1):10]
  prvals = NULL

  tnet = getNetworkFromTissue(tissue=tissue,which.one="micro19K")
  outtissues = NULL
  for(othertissue in lefttissues){
    othernet = getNetworkFromTissue(tissue=othertissue,which.one="micro19K")
    file.in = paste0("supplementary/rdsnets/micro19K/",
                     othertissue,".mic.net.19K.rds_vs_",tissue,
                     ".mic.net.19K.rds_preserv.rds","Zsummary.csv")
    if(!file.exists(file.in))
      getPreservationStatistics(networks=list(othernet,tnet),
                                tissues=c(othertissue,tissue),
                                paste0("supplementary/rdsnets/micro19K/",othertissue,".mic.net.19K.rds_vs_",tissue,
                                       ".mic.net.19K.rds_preserv.rds"))

    pr = read.csv(file=file.in,stringsAsFactors=F)
    if(sum(pr$tissue == tissue & pr$color == module) == 0)
      stop(paste0("When checking preservation of ",module," from tissue ",tissue,
                  " at tissue ",othertissue," we cant find the module"))
    prvals = c(prvals,pr[pr$tissue == tissue & pr$color == module,othertissue])
    outtissues = c(outtissues,othertissue)
  }

  for(othertissue in righttissues){
    othernet = getNetworkFromTissue(tissue=othertissue,which.one="micro19K")
    file.in = paste0("supplementary/rdsnets/micro19K/",
                     tissue,".mic.net.19K.rds_vs_",othertissue,
                     ".mic.net.19K.rds_preserv.rds","Zsummary.csv")
    if(!file.exists(file.in))
      getPreservationStatistics(networks=list(tnet,othernet),
                                tissues=c(tissue,othertissue),
                                paste0("supplementary/rdsnets/micro19K/",
                                       tissue,".mic.net.19K.rds_vs_",othertissue,
                                       ".mic.net.19K.rds_preserv.rds"))
    pr = read.csv(file=file.in,stringsAsFactors=F)
    if(sum(pr$tissue == tissue & pr$color == module) == 0)
      stop(paste0("When checking preservation of ",module," from tissue ",tissue,
                  " at tissue ",othertissue," we cant find the module"))
    prvals = c(prvals,pr[pr$tissue == tissue & pr$color == module,othertissue])
    outtissues = c(outtissues,othertissue)
  }
  names(prvals) = c(lefttissues,righttissues)
  return(list(min=min(prvals),max=max(prvals),mean=mean(prvals),all=prvals,tissues=outtissues))
}

getPreservationStatistics = function(networks,tissues,prev.data.file){

  net.objects <- list()

  if(typeof(networks[[1]]) == "character"){
    for(index in 1:length(networks)){
      print(paste0("Loading network ",networks[index]," for tissue ",tissues[index]))
      net.objects[[tissues[index]]] <- readRDS(networks[[index]])
    }
  }else{
    for(index in 1:length(networks)){
      print(paste0("Loading network for tissue ",tissues[index]))
      net.objects[[tissues[index]]] <- networks[[index]]
    }

  }



  Zsummary <- NULL
  MedianRank <- NULL

  for(tissue1 in tissues){
    Z.tmp.list <- list(NULL)
    MR.tmp.list <- list(NULL)
    for(tissue2 in tissues){
      ## loop through columns
      cat(tissue1, "\t", tissue2, "\n")
      if(tissue2 == tissue1) next() ## skip self-self comparison
      if( tissue1 < tissue2 ){
        fn	<- paste(tissue1, "vs", tissue2,sep="")
        ref <- 1
        test <- 2
      } else {
        ## swap around
        fn <- paste(tissue2, "vs", tissue1,sep="")
        ref <- 2
        test <- 1
      }
      mp = readRDS(prev.data.file)
      ### other statistics can be exstracted here
      #these statistics for instance can be exstracted:
      #Zsummary.pres","Zconnectivity.pres","Zdensity.pres","Z.cor.kIM"

      Z.tmp <- mp$preservation$Z[[ref]][[test]][,"Zsummary.pres",drop=F]
      rownames(Z.tmp) <- paste(tissue1, "_", rownames(Z.tmp),sep="")
      colnames(Z.tmp) <- tissue2
      Z.tmp.list[[ tissue2 ]] <- Z.tmp

      MR.tmp <- mp$preservation$observed [[ ref ]][[ test ]][ , "medianRank.pres", drop=F]
      rownames(MR.tmp) <- paste(tissue1, "_", rownames(MR.tmp),sep="")
      colnames(MR.tmp) <- tissue2
      MR.tmp.list[[ tissue2 ]] <- MR.tmp
      rm(mp, Z.tmp, MR.tmp, tissue2, fn, ref, test)
    }

    Z.tmp.list <- Z.tmp.list[ which(names(Z.tmp.list) != "") ]
    Z.tmp.list <- cbind( do.call("cbind", Z.tmp.list), tissue1=NA )
    colnames(Z.tmp.list)[ colnames(Z.tmp.list)=="tissue1" ] <- tissue1

    Zsummary <- rbind(Zsummary, Z.tmp.list[ , tissues])
    MR.tmp.list <- MR.tmp.list[ which(names(MR.tmp.list) != "") ]
    MR.tmp.list <- cbind( do.call("cbind", MR.tmp.list), tissue1=NA )
    colnames(MR.tmp.list)[ colnames(MR.tmp.list)=="tissue1" ] <- tissue1
    MedianRank <- rbind(MedianRank, MR.tmp.list[ , tissues])
  }

  # 230 i.e. ten extra
  # gold module consists of 1000 randomly selected genes
  print(dim(Zsummary))
  w <- grep("gold", rownames(Zsummary))
  round( Zsummary[w, ], 2 )

  MedianRank[w, ]
  Zsummary <- Zsummary[ -w, ]
  MedianRank <- MedianRank[ -w, ]

  # get the real module sizes
  moduleColors <- sapply(net.objects, function(obj) obj$moduleColors)
  mat = NULL
  mat[[tissues[1]]] = table(moduleColors[,1])
  mat[[tissues[2]]] = table(moduleColors[,2])

  #mat <- apply( moduleColors[, 1:length(tissues)], 2, function(x) table(x) )
  mat <- data.frame(as.matrix(unlist(mat)))
  rownames(mat) <- gsub("[.]","_",rownames(mat))
  kk <- match(rownames(mat),rownames(Zsummary))
  Zsummary <- Zsummary[kk,]
  kk <- match(rownames(mat),rownames(MedianRank))
  MedianRank <- MedianRank[kk,]

  Zsummary <- cbind( mat[,1],Zsummary)
  colors.col <- NULL
  tissue.col <- NULL
  for(index in 1:length(rownames(mat))){
    values <- str_split(rownames(mat)[index],"_")
    colors.col[index] <- values[[1]][2]
    tissue.col[index] <- values[[1]][1]
  }
  Zsummary <- cbind(tissue.col,colors.col,Zsummary)
  colnames( Zsummary ) <- c("tissue","color","size",tissues)

  MedianRank <- cbind( mat[,1],MedianRank )
  colnames( MedianRank ) <- c("size",tissues)
  write.csv( Zsummary, file=paste0(prev.data.file,"Zsummary.csv"), quote=F )
  write.csv( MedianRank, file=paste0(prev.data.file,"MedianRank.csv"), quote=F )

}


UserGOenrichment <- function(net.file,net.name=NULL,which.one="exonic",do.plot=F){
  #Variable initialization
  if(typeof(net.file) == "character")
    net <- readRDS(net.file)
  else
    net = net.file
  if(is.null(net.name))
    net.name <- str_extract(net.file,"network[0-9]+.+rds$")
  out.file <- paste0(net.name,".USER_terms.csv")

  cat("Passing",length(net$moduleColors)," gene symbols to user enrichment\n")

  cat("Writting everything to",out.file,"\n")

  if(which.one != "braineac")
    genes = fromAny2GeneName(names(net$moduleColors))
  else
    genes = newmic.fromtID2GeneName(names(net$moduleColors))

  enrichments = WGCNA::userListEnrichment(genes,
                                          net$moduleColors,
                                          fnIn=NULL,
                                          useBrainLists = TRUE,
                                          useBrainRegionMarkers=T,
                                          useImmunePathwayLists=T,
                                          useBloodAtlases=T,
                                          useStemCellLists=T,
                                          usePalazzoloWang=T,
                                          nameOut=out.file,
                                          outputCorrectedPvalues=TRUE)



  if(do.plot){
    pdf(paste0(net.file,".go_user_zscores.pdf"),width=8,height=6)
    plotGOZScoresUserList(enrichments,net,
                          title=paste0("User enrichment. Sum z-s (p-val < 10e-2), ",
                                       net.name))
    dev.off()
  }
  write.csv(enrichments$sigOverlaps,out.file)
  #saveGOTermDefinition(enrichments,out.file)

}



annotate = function(net,subnets=T,organism="hsapiens",ensembl=F){
  if(typeof(net) == "character")
    net = readRDS(net)
  if(is.null(names(net$moduleColors)))
    names(net$moduleColors) = names(net$adjacency)
  cat("Annotating network",basename(net$file),"\n")
  net$go = getGProfilerOnNet(net.file=net,out.file=NULL,ensembl=ensembl)
  notHuman = ifelse(organism == "hsapiens",F,T)
  net$ct = genAnnotationCellType(net.in=net,return.processed=F,which.one="new",notHuman=notHuman)
  # if(subnets){
  #   subnetsgo = NULL
  #   subnets = net$subnets
  #   for(job in names(subnets)){
  #     cat("Annotating",job,"net\n")
  #     localnet = NULL
  #     nparts = length(subnets[[job]]$partitions)
  #     localnet$moduleColors = subnets[[job]]$partitions[[nparts]]
  #     subnetsgo[[job]] = getGProfilerOnNet(net.file=net,...)
  #   }
  # }
  net
}


#This function generates a csv file with the Gene Ontology enrichment
#signals. It is based on the gProfileR R package (see documentation for the
#package). Main parameters are
#net.file			: full path name for the net file
#filter				: ontologies to test for enrichemtn (see gProfileR docs)
#exclude.iea		: do not use Inferred electronic annotations (see gProfileR docs)
#correction.method	: method for multiple testing correction (see gProfileR docs)
#out.file			: full name for the csv file where to store results
#incluster			: ignore that one
getGProfilerOnNet <- function(net.file,
                              filter=c("GO","KEGG","REAC"),
                              ensembl=TRUE,
                              exclude.iea=T,
                              organism = "hsapiens",
                              correction.method="gSCS",
                              out.file=NULL,...){

  if(typeof(net.file) == "character")
    net <- readRDS(net.file)
  else
    net = net.file

  modules = unique(net$moduleColors)

  background <- names(net$moduleColors)
  if(ensembl)
    background <- fromEnsembl2GeneName(background)

  all.genes <- NULL
  for(module in modules){

    genes <- names(net$moduleColors)[net$moduleColors == module]
    if(ensembl)
      genes <- fromEnsembl2GeneName(genes)
    all.genes[[module]] <- genes
  }

  go <- gProfileR::gprofiler(all.genes,
                             correction_method=correction.method,
                             #custom_bg=background,
                             src_filter=filter,
                             organism=organism,
                             exclude_iea=exclude.iea)
  if(nrow(go) == 0)
    return(go)
  #png_fn = paste0(out.file,".png"),
  #no_isects=T)
  go.out = cbind(go,rep(NA,nrow(go)))
  colnames(go.out) = c(colnames(go),"IC")

  #	#Now we will add our own column with the information content of each term
  #	#So it is possible to evaluate them
  cat("Generating IC-BP")
  loadICsGOSim("BP")
  mask = go$domain == "BP"
  bp.terms.go = go$term.id[mask]
  bp.ic.go = IC[bp.terms.go]
  go.out$IC[mask] = bp.ic.go

  cat("Generating IC-MF")
  loadICsGOSim("MF")
  mask = go$domain == "MF"
  mf.terms.go = go$term.id[mask]
  mf.ic.go = IC[mf.terms.go]
  go.out$IC[mask] = mf.ic.go

  cat("Generating IC-CC")
  loadICsGOSim("CC")
  mask = go$domain == "CC"
  cc.terms.go = go$term.id[mask]
  cc.ic.go = IC[cc.terms.go]
  go.out$IC[mask] = cc.ic.go
  #
  #	data("ICsMFhumanall")
  #	data("ICsCChumanall")
  if(!is.null(out.file)){
    write.csv(go.out,out.file)
    cat("Obtained",nrow(go),"terms saved at",out.file,"\n")
  }else{
    cat("Obtained",nrow(go),"terms\n")
  }
  return(go.out)
}

loadICsGOSim = function(onto){
  stopifnot(onto %in% c("BP","MF","CC"))
  if(onto == "BP")
    load(system.file("", "ICsBPhumanall.rda", package = "CoExpNets"),.GlobalEnv)
  else if(onto == "MF")
    load(system.file("", "ICsMFhumanall.rda", package = "CoExpNets"),.GlobalEnv)
  else
    load(system.file("", "ICsCChumanall.rda", package = "CoExpNets"),.GlobalEnv)
}

getICReport = function(gprof){
  ontologies = c("BP","CC","MF")
  ics = NULL
  for(onto in ontologies){
    loadICsGOSim(onto)
    terms = gprof$term.id[gprof$domain == onto]
    ics[[onto]] = IC[terms]
    names(ics[[onto]]) = terms
  }
  return(ics)
}

plotGOZScoresUserList <- function(goresult,net,title="Sum z-scores "){
  #Transforming p values in log10 scale
  logs <- -log10(goresult$sigOverlaps$CorrectedPvalues)
  modules <- unique(goresult$sigOverlaps$InputCategories)

  mod.colors <- as.character(modules)

  print(modules)
  z.scores <- NULL
  n.terms <- NULL
  mod.labels <- NULL
  c.modules <- 0


  for(module in modules){
    indexes <- which(goresult$sigOverlaps$InputCategories == module)
    print(indexes)
    good.p.values <- logs[indexes]
    c.modules <- c.modules + 1
    z.scores[c.modules] <- sum(good.p.values)
    n.terms[c.modules] <- length(good.p.values)
    mod.labels[c.modules] <- paste0(n.terms[c.modules],"/",
                                    table(net$moduleColors == module)[2])
  }

  data.to.plot <- as.data.frame(cbind(1:c.modules,z.scores,mod.colors,mod.labels,n.terms),
                                stringsAsFactors=FALSE)
  data.to.plot <- data.to.plot[order(-z.scores,-n.terms),]

  x <- barplot(height=as.numeric(data.to.plot$z.scores),names.arg=data.to.plot$mod.labels,
               width=as.numeric(data.to.plot$n.terms),col=data.to.plot$mod.colors,
               main=title, xaxt="n",
               ylab="-log10(pval)")
  text(cex=1, x=x-2, y=-1.25, pos=1, offset=2, data.to.plot$mod.labels, xpd=TRUE, srt=45)
  legend("topright", title="Module colors",
         legend=data.to.plot$mod.colors,
         fill=data.to.plot$mod.colors,cex=0.5)
}

getFriendlyNameForCategories = function(cats){
  friendly.cats = cats
  existing.cats = cats[cats %in% coexp.cell.categories]
  friendly.cats[cats %in% coexp.cell.categories] =
    coexp.long.friendly.cell.categories[match(existing.cats,coexp.cell.categories)]
  return(friendly.cats)
}

getCategoriesForGene = function(gene){
  #Load the gene lists
  if (!(exists("BrainLists")))
    BrainLists = NULL
  data("BrainLists", envir = sys.frame(sys.nframe()))
  return(BrainLists[,2][BrainLists[,1] == gene])
}

getFriendlyNameForColumnCategories = function(cats){
  friendly.cats = cats
  existing.cats = cats[cats %in% coexp.cell.types]

  friendly.cats[cats %in% coexp.cell.types] =
    coexp.long.friendly.cell.categories[match(existing.cats,coexp.cell.types)]
  return(friendly.cats)
}

#This method is useful for generating cell type enrichment per module from a network
#It will generate a csv ending with celltype.csv with a squared matrix (modules in rows
#and tested cell type markers in columns)
#All can be left to its default value except done for the net.in parameter which
#must have a full path for the network RDS file.
#The csv will be generated at the same folder than the network
#The pdf file with the signals will be generated at plot.file location (can be used as a prefix)
#for the pdf file as in the default
cellTypeByModule = function(tissue="None",
                            do.plot=T,
                            which.one="new",
                            plot.file=NULL,
                            net.in=NULL,
                            legend=NULL,
                            threshold=20,
                            return.processed=T,
                            use.grey=F,
                            display.cats=NULL){ #This last flag FALSE when you want the raw p-values

  cat("Entering cellTypeByModule with",which.one,"\n")
  #if(do.plot & is.null(plot.file))
  #  stop(paste0("You have to specify the plot file if you want a heatmap plotted\n"))


  if(is.null(net.in)){
    net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
    modules = getModulesFromTissue(tissue,which.one)
    enrichment = getCellTypeFromTissue(which.one=which.one,tissue=tissue)
    enrf = findUserCT(which.one=which.one,tissue=tissue)
    userEnrichment = data.frame()
    if(!is.null(enrf)){
      if(file.exists(enrf))
        userEnrichment = read.csv(enrf)
    }
    is.any.network = F
  }else{
    if(typeof(net.in) == "character")
      net = readRDS(net.in)
    else
      net = net.in
    modules = unique(net$moduleColors)
    is.any.network = T
    file.name = paste0(net.in,".USER_terms.csv")
    #if(!file.exists(file.name)){
    cat("Generating new user enrichment into ",file.name,"\n")
    UserGOenrichment(net.file=net,
                     net.name=net.in)
    cat("Done generating new User enrichment\n")
    #}
    userEnrichment = read.csv(file.name,stringsAsFactors=F)

  }
  if(!use.grey)
    modules = modules[modules != "grey"]


  external.ref = read.csv(paste0(system.file("", "", package = "CoExpNets"),
                                 "/cell_type_TableS1.csv"),stringsAsFactors=F)

  all.cell.types = c(coexp.cell.types,
                     coexp.ukbec.cell.types,
                     paste0("External-",colnames(external.ref)))
  cell.type.data = matrix(ncol=length(modules),nrow=length(all.cell.types))
  cell.type.data[,] = 1
  rownames(cell.type.data) = c(getFriendlyNameForColumnCategories(coexp.cell.types),
                               coexp.ukbec.cell.types,colnames(external.ref))
  colnames(cell.type.data) = modules

  #WGCNA enrichment
  if(nrow(userEnrichment) > 0){
    for(module in modules){
      i = match(module,modules)
      for(cell.type in coexp.cell.types){
        j = match(cell.type,coexp.cell.types)
        p.val = as.numeric(userEnrichment$CorrectedPvalues[userEnrichment$InputCategories %in% module &
                                                             userEnrichment$UserDefinedCategories %in% coexp.cell.categories[j]])
        if(length(p.val) > 0)
          cell.type.data[j,i] = p.val #-log10(p.val)
      }
    }
  }

  input.file.names = c(paste0(system.file("", "cell_type_TableS1.csv", package = "CoExpNets"),"/",
                              coexp.ukbec.cell.files),
                       paste0(system.file("", "cell_type_TableS1.csv", package = "CoExpNets"),"/",
                              "cell_type_TableS1.csv.",colnames(external.ref),".txt"))
  cell.types.i.f = c(coexp.ukbec.cell.types,colnames(external.ref))

  last.cell.types = c(coexp.ukbec.cell.types,colnames(external.ref))

  if(nrow(userEnrichment) > 0)
    for(i in 1:nrow(userEnrichment)){
      #print(enrichment)
      module = userEnrichment$InputCategories[i]
      #print(module)
      for(j in 1:length(input.file.names)){
        if(grepl(input.file.names[j],userEnrichment$UserDefinedCategories[i]))
          break
      }
      category = cell.types.i.f[j]
      cell.type.data[category,module] = userEnrichment$CorrectedPvalues[i]
    }


  #But before, rename them
  rownames(cell.type.data) = c(getFriendlyNameForColumnCategories(coexp.cell.types),
                               coexp.ukbec.cell.types,paste0(colnames(external.ref),"-External"))
  unprocessed.cell.type.data = cell.type.data
  cell.type.data = cell.type.data[,apply(cell.type.data,2,function(x){ any(x < 1)}),drop=FALSE]
  cell.type.data = -log10(cell.type.data)
  cell.type.data[is.infinite(cell.type.data)] = max(cell.type.data[!is.infinite(cell.type.data)])

  if(threshold > 0)
    cell.type.data[cell.type.data > threshold] = threshold
  #Order by cell type

  cell.type.data = cell.type.data[order(rownames(cell.type.data)),,drop=FALSE]

  #Now filter
  if(!is.null(display.cats)){
    cell.type.data = cell.type.data[rownames(cell.type.data) %in% display.cats,]
  }


  if(do.plot & (ncol(cell.type.data) >= 2)){
    if(is.null(legend))
      legend = tissue
    if(!is.null(plot.file))
      pdf(plot.file,width=15,height=8)
    gplots::heatmap.2(cell.type.data[apply(cell.type.data,1,function(x){ any(x > 2)}),],
                      trace="none",
                      col=heat.colors(100)[100:1],
                      cexCol=0.7,
                      cexRow=0.7,
                      Rowv=F,
                      Colv=T,
                      main=paste0("Cell type enrichment for ",legend),
                      key=F,srtCol=45,dendrogram="none",
                      margins=c(6,20))

    if(!is.null(plot.file))
      dev.off()
  }
  if(return.processed)
    return(cell.type.data)
  return(unprocessed.cell.type.data)
}

genAnnotationCellType = function(tissue="None",
                                 which.one="new",
                                 markerspath=system.file("ctall", "",
                                                         package = "CoExpNets"),
                                 net.in=NULL,
                                 legend=NULL,
                                 doheatmap=F,
                                 notHuman=F,
                                 plot.file=NULL,
                                 threshold=20,
                                 return.processed=T,
                                 getMarkerSize=F,
                                 getOnlyOverlap=F,
                                 getMyOwnTest=F,
                                 applyFDR=F,
                                 display.cats=NULL){ #This last flag FALSE when you want the raw p-values

  cat("Entering genAnnotationCellType with",which.one,"\n")

  if(which.one == "new"){
    if(typeof(net.in) == "character")
      net = readRDS(net.in)
    else
      net = net.in
  }else{
    net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  }
  modules = unique(net$moduleColors)
  if(notHuman)
    names(net$moduleColors) = toupper(names(net$moduleColors))
  #So the 1st heatmap
  #will have a column for each cell.types element and a row for
  #each module and we will show -log10(p-values) in a scale

  files = list.files(path=markerspath,full.names = T)
  files = files[grep(pattern=".txt$",files)]
  if(getMarkerSize){
    sizes = lapply(files,function(x){
      nrow(read.delim(x,header=T))
    })
    return(cbind(markerset=files,genes=sizes))
  }


  markernames = gsub(".txt","",basename(files))


  ctypes = markernames
  ctypedata = matrix(ncol=length(modules),nrow=length(files))
  ctypedata[,] = 1
  rownames(ctypedata) = ctypes
  colnames(ctypedata) = modules

  if(getOnlyOverlap){
    myfunction = Vectorize(function(m,f){
      submarkers = read.delim(f,stringsAsFactors=F,header=T)[,1]
      #cat("Module is",m,"file is",f,"\n")
      #print(names(net$moduleColors)[net$moduleColors == m])
      #print(submarkers)
      return(sum(submarkers %in% names(net$moduleColors)[net$moduleColors == m]))
    })

    overlap = outer(modules,files,myfunction)
    colnames(overlap) = gsub(".txt","",basename(files))
    rownames(overlap) = modules
    if(applyFDF){
      overlap = p.adjust(overlap,method="fdr")
    }
    return(overlap)

  }

  if(getMyOwnTest){
    myfunction = Vectorize(function(m,f){
      submarkers = read.delim(f,stringsAsFactors=F,header=T)[,1]
      ov = sum(submarkers %in% names(net$moduleColors)[net$moduleColors == m])
      nm = sum(net$moduleColors == m)

      #cat("Module is",m,"file is",f,"\n")
      #print(names(net$moduleColors)[net$moduleColors == m])
      #print(submarkers)
      return(testGeneSet(n.module=nm,total.specific=length(submarkers),
                         total.net=length(net$moduleColors),n.module.and.specific=ov)$p.value)
    })

    overlap = outer(modules,files,myfunction)
    colnames(overlap) = gsub(".txt","",basename(files))
    rownames(overlap) = modules
    return(overlap)

  }



  all.gene.names = names(net$moduleColors)
  all.gene.names = fromAny2GeneName(names(net$moduleColors))

  tmp.enr.f = paste0(markerspath,"enrichment_tmp.csv")
  if(file.exists(tmp.enr.f))
    file.remove(tmp.enr.f)
  ukbec.en = WGCNA::userListEnrichment(all.gene.names,net$moduleColors,files,
                                       nameOut=tmp.enr.f)
  if(file.exists(tmp.enr.f)){
    enrichment = read.csv(tmp.enr.f,
                          stringsAsFactors=F)

    for(i in 1:nrow(enrichment)){
      module = enrichment$InputCategories[i]
      for(j in 1:length(files)){
        if(grepl(ctypes[j],enrichment$UserDefinedCategories[i]))
          break
      }
      category = ctypes[j]
      ctypedata[category,module] = enrichment$CorrectedPvalues[i]
    }
  }
  rownames(ctypedata) = gsub("cell_type_TableS1.csv.","",rownames(ctypedata))
  uctypedata = ctypedata
  ctypedata = ctypedata[,apply(ctypedata,2,function(x){ any(x < 1)}),drop=FALSE]
  ctypedata = -log10(ctypedata)
  ctypedata[is.infinite(ctypedata)] = max(ctypedata[!is.infinite(ctypedata)])

  if(threshold > 0)
    ctypedata[ctypedata > threshold] = threshold
  #Order by cell type

  ctypedata = ctypedata[order(rownames(ctypedata)),,drop=FALSE]

  if(doheatmap){
    if((ncol(ctypedata) >= 2)){
      if(is.null(legend))
        legend = tissue
      if(!is.null(plot.file))
        pdf(plot.file,width=15,height=8)
      heatmap.2(ctypedata[apply(ctypedata,1,function(x){ any(x > 2)}),],
                trace="none",
                col=heat.colors(100)[100:1],
                cexCol=0.7,
                cexRow=0.7,
                Rowv=F,
                Colv=T,
                main=paste0("Cell type enrichment for ",legend),
                key=F,srtCol=45,dendrogram="none",
                margins=c(6,20))

      if(!is.null(plot.file))
        dev.off()
    }
  }

  if(return.processed)
    return(ctypedata)
  return(uctypedata)
}


