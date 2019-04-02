

#Use system.file("extdata", "2010.csv", package = "testdat", mustWork = TRUE)
#' Title
#'
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
fromAny2Ensembl = function(genes){
  if(substr(genes[1],1,1) == "X"){
    return(fromGeneName2Ensembl(fromXtIDToGeneSymbols19K(genes)))
  }else if(substr(genes[1],1,4) != "ENSG"){
    return(fromGeneName2Ensembl(genes))
  }

  #Still one possibility, does it have versions?
  genes = gsub("\\.[0-9]+","",genes)

  return(genes)
}

#' Title
#'
#' @param genes
#' @param use38
#'
#' @return
#' @export
#'
#' @examples
fromGeneName2EnsemblBM = function(genes,use38=T){

  if(use38){
    ensembl <- useMart(host="www.ensembl.org",
                       biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_name"
  }else{
    ensembl <- useMart(host="jun2013.archive.ensembl.org",
                       biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_id"
  }

  attributes <- c(external.gene.att,"ensembl_gene_id")
  genes.with.name = getBM(attributes=attributes, filters="hgnc_symbol", values=genes,mart=ensembl)
  cat("From",length(genes),"gene IDs we got",nrow(genes.with.name),"genes with Ensemble name\n")
  #if(nrow(genes.with.name) >= length(genes))
  thematch = match(genes,genes.with.name$external_gene_name)
  outgenes = genes.with.name[,2][thematch]
  outgenes[is.na(thematch)] = genes[is.na(thematch)]
  cat("fromGeneName2EnsemblBM, couldn't convert",sum(is.na(thematch)),"genes\n")
  return(outgenes)



  #Lets give it another try
  left.genes = genes[!(genes %in% genes.with.name$hgnc_symbol)]
  cat("Trying now with the left",length(left.genes),"genes in Vega\n")
  attributes <- c("clone_based_vega_gene_name","ensembl_gene_id")
  genes.with.vega = getBM(attributes=attributes, filters="clone_based_vega_gene_name",
                          values=left.genes,mart=ensembl)
  cat("From",length(left.genes),"gene IDs we got",nrow(genes.with.vega),"genes with VEGA clones\n")
  genes[match(genes.with.vega[,1],genes)] = genes.with.vega[,2]
  genes[match(genes.with.name[,1],genes)] = genes.with.name[,2]
  return(genes)
}

#' Title
#'
#' @param genes
#' @param use38
#'
#' @return
#' @export
#'
#' @examples
fromEnsembl2GeneNameBM = function(genes,use38=T){
  if(use38){
    ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_name"
  }else{
    ensembl <- useMart(host="jun2013.archive.ensembl.org",
                       biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_id"
  }

  attributes <- c("ensembl_gene_id",external.gene.att)
  genes.with.name = getBM(attributes=attributes, filters="ensembl_gene_id", values=genes,mart=ensembl)
  cat("From",length(genes),"Ensembl IDs we got",nrow(genes.with.name),"genes with external gene name\n")
  thematch = match(genes,genes.with.name$ensembl_gene_id)
  outgenes = genes.with.name[,2][thematch]
  if(sum(is.na(thematch)) == 0)
    return(outgenes)
  cat("Couldn't conver",sum(is.na(thematch)),"genes\n")
  outgenes[is.na(thematch)] = genes[is.na(thematch)]
  return(outgenes)


  thematch = match(genes,genes.with.name$external_gene_name)
  outgenes = genes.with.name[,2][thematch]
  outgenes[is.na(thematch)] = genes[is.na(thematch)]
  cat("fromGeneName2EnsemblBM, couldn't convert",sum(is.na(thematch)),"genes\n")
  return(outgenes)

}


#' Title
#'
#' @param genes
#' @param table.file
#' @param ignore.unknown
#' @param which.are.unknown
#'
#' @return
#' @export
#'
#' @examples
fromEntrez2Ensembl <- function(genes,
                               table.file=biotype.all.file,
                               ignore.unknown=FALSE,
                               which.are.unknown=FALSE){

  function.table <- read.table(table.file,header=TRUE,stringsAsFactors=FALSE,sep=",")
  gene.names <- function.table$ensembl_gene_id[match(genes,function.table$external_gene_id)]

  if(which.are.unknown)
    return(is.na(gene.names))

  if(ignore.unknown){
    gene.names <- na.omit(gene.names)
  }else{
    unnamed.genes.count <- sum(is.na(gene.names))
    if(unnamed.genes.count > 0){
      print(paste0("Genes without name ",unnamed.genes.count))
      if(unnamed.genes.count > 5)
        print(paste0("Genes ",paste0(genes[is.na(gene.names)][1:5],
                                     collapse=", "), " and more... don't have a name"))
      else
        print(paste0("Genes ",paste0(genes[is.na(gene.names)],collapse=", "), " don't have a name"))
      gene.names[is.na(gene.names)] <- genes[is.na(gene.names)]
    }

  }

  return(gene.names)
}

#' Title
#'
#' @param genes
#' @param ignore.unknown
#' @param which.are.unknown
#'
#' @return
#' @export
#'
#' @examples
fromGeneName2Ensembl <- function(genes,ignore.unknown=FALSE,which.are.unknown=FALSE){
  if(!exists("coexp.utils.gene.names.conv.table"))
    loadConvTable()

  gene.names <- coexp.utils.gene.names.conv.table$Ensembl[match(genes,coexp.utils.gene.names.conv.table$Gene)]

  if(which.are.unknown)
    return(is.na(gene.names))

  if(ignore.unknown){
    gene.names <- na.omit(gene.names)
  }else{
    unnamed.genes.count <- sum(is.na(gene.names))
    if(unnamed.genes.count > 0){
      print(paste0("Genes without name ",unnamed.genes.count))
      print(paste0("Genes ",paste0(genes[is.na(gene.names)][1:unnamed.genes.count],collapse=", "),
                   " don't have a name"))
      gene.names[is.na(gene.names)] <- genes[is.na(gene.names)]
    }
  }
  return(gene.names)
}

#' Title
#'
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
fromAny2GeneName = function(genes){
  if(substr(genes[1],1,1) == "X"){
    return(fromXtIDToGeneSymbols19K(genes))
  }else if(substr(genes[1],1,4) == "ENSG"){
    return(fromEnsembl2GeneName(genes))
  }
  return(genes)
}

#' Title
#'
#' @param xids
#'
#' @return
#' @export
#'
#' @examples
fromXtIDToGeneSymbols19K = function(xids){

  trans.table <- read.csv(paste0(system.file("", "", package = "CoExpNets"),
                                 "annot_19K.csv"),
                          stringsAsFactors=F)
  gene.symbols = trans.table$Gene_Symbol[match(xids,trans.table$XtID)]
  gene.symbols[is.na(gene.symbols)] = xids[is.na(gene.symbols)]
  return(gene.symbols)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
loadConvTable = function(){

  the.table = read.table(system.file("", "EnsemblNamesGRCh38.txt", package = "CoExpNets"),
                         header=F,stringsAsFactors=F)
  colnames(the.table) = c("Ensembl","Gene")
  coexp.utils.gene.names.conv.table <<- the.table
}

fromEnsembl2GeneName <- function(genes,ignore.unknown=FALSE,which.are.unknown=FALSE){

    loadConvTable()

  gene.names <- coexp.utils.gene.names.conv.table$Gene[match(genes,coexp.utils.gene.names.conv.table$Ensembl)]

  if(which.are.unknown)
    return(is.na(gene.names))

  if(ignore.unknown){
    gene.names <- na.omit(gene.names)
  }else{
    unnamed.genes.count <- sum(is.na(gene.names))
    if(unnamed.genes.count > 0){
      cat("Genes without name ",unnamed.genes.count,"\n")
      if(unnamed.genes.count > 3)
        unnamed.genes.count = 10
      cat("Genes ",paste0(genes[is.na(gene.names)][1:unnamed.genes.count],collapse=", "), " and more... don't have a name\n")
      gene.names[is.na(gene.names)] <- genes[is.na(gene.names)]
    }
  }
  return(gene.names)
}

