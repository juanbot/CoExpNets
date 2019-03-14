
# Introduction to using co-expression networks with CoExp suite

This document is an introduction to co-expression networks. We will learn how to use co-expression networks (GCNs). A gene co-expression network, as we use them in CoExp, is a network of genes created with the WGCNA (Weighted Co-expression Networks) R package, see <https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/> in combination with the km2gcn R package, see <https://www.ncbi.nlm.nih.gov/pubmed/28403906>. This package has been materialised into the CoExpNets package, hosted in GitHub, <http://github.com/juanbot/CoExpNets>. 

This document has two purposes. The first one is to serve as a tutorial on how to create and annotate GCNs from gene expression data adequately cleaned for such purpose. For this, we will need the CoExpNets package, that can be installed as follows

```r
devtools::install_github("juanbot/CoExpNets")
```

The second one is to illustrate how to use the database of GCNs that come in separate resources, for research purposes. Currently, there are three packages available, `CoExpGTEx`, `CoExp10UKBEC` and `CoExpROSMAP` with 47, 10 and 4 co-expression networks respectively. They can be found at <http://github.com/juanbot/CoExpGTEx>, <http://github.com/juanbot/CoExp10UKBEC> and <http://github.com/juanbot/CoExpROSMAP> respectively.


# CoExpNets: resource for co-expression networks and gene expression profiling

For the rest of things we´ll use here, we need to install the ROSMAP package, together with CoExpNets, so besides installing CoExpNets as indicated above, we have to do 

```r
devtools::install_github("juanbot/CoExpROSMAP")
```

No we can start using the software by doing

```r
library(CoExpNets)
library(CoExpROSMAP)
CoExpROSMAP::initDb()
```

And now we are ready to start. A number of networks have been added to the database.


If you want to compile this documentation itself in your computer in order to check whether your installation was correct, then the rest of chunks of code for this documentation will be run in order for the documentation to be generated so please change the first line of the chunk above accordingly before compiling again.

```r
for(data.family in getNetworkCategories()){
	cat(paste0("Family of gene expression dataset ",data.family," available\n"))
  tissues = getAvailableNetworks(data.family)
	for(tissue in tissues)
		cat(paste0("Gene expression profiling dataset for tissue ",tissue,"\n is ",
		           getExprDataFromTissue(which.one=data.family,tissue=tissue,only.file = T)," available\n"))
}
```

In order to understand this we need to check the page <https://www.synapse.org/#!Synapse:syn3219045> for information on the ROS/MAP study and how to access its data. This is a study on multi-omics data from frontal cortex samples of donors with various levels of severity of Alzheimer's disease.

Let's get the sample covariates so we can have a look at the values for the samples.
```r
#In ROS/MAP we have four networks, one for all the samples and three different 
#networks segreating by different levels of severity for cognitive decline 
#(not AD, probable AD and surely AD, AD stands for Alzheimer's disease)
getAvailableNetworks("CoExpROSMAP")
#One of the covariates of this dataset is cognitive decline (cogdx) and it is value
#between 0 and 5 being 5 highest severity
covs = CoExpROSMAP::getCovariates(tissue="allsamples",which.one="CoExpROSMAP")
oldpar = par()
par(mfrow=c(2,2))
barplot(table(covs$gender),main="Gender in ROS/MAP")
hist(covs$age,main="Age in ROS/MAP")
barplot(table(covs$cogdx),main="Cognitive decline in ROS/MAP")
barplot(table(covs$braaksc),main="Braak stage in ROS/MAP")
par(oldpar)
```


In the same way, we can check that there is a co-expression network available for each of these datasets. 
```r
for(data.family in getNetworkCategories()){
	cat(paste0("Family of gene expression dataset ",data.family," available\n"))
  tissues = getAvailableNetworks(data.family)
	for(tissue in tissues)
		cat(paste0("Co-expression network for tissue ",tissue,"\n is ",
		           getNetworkFromTissue(which.one=data.family,
		                                tissue=tissue,only.file = T),
		           " available\n"))
}
```

## Accessing gene expression data

So, for example, if we wanted to access the gene expression profiling of `ad` tissue, we can do
```r
expr.data = getExprDataFromTissue(which.one="CoExpROSMAP",tissue="ad")
cat("The dataset has ",dim(expr.data)," samples and genes respectively\n")
cat("First gene names are",paste0(colnames(expr.data)[1:3],collapse=","),"\n")
cat("First sample IDs are",paste0(rownames(expr.data)[1:3],collapse=", "),"\n")
```

Note that gene IDs depend on the gene expression profiling family. For ROS/MAP all gene names are [Ensembl](http://www.ensembl.org/index.html) IDs. We keep the familiy dependent sample ID too when possible.

We can convert between Ensembl IDs and gene symbols by using 
```r
fromEnsembl2GeneName(gsub("\\.[0-9]+$","",colnames(expr.data)[1:5]))
```
And back to Ensembl IDs
```r
fromGeneName2Ensembl(fromEnsembl2GeneName(gsub("\\.[0-9]+$","",
                                               colnames(expr.data)[1:5])))
```
Note that this is done by using internal files. The most reliable manner in which we can manipulate Ensembl IDs is going directly to biomart. But at the moment this falls outside of the scope of this tutorial.

## Accessing network data and annotations

We can access a network and get its genes

```r
net = getNetworkFromTissue(which.one="CoExpROSMAP",tissue="ad")
cat("The network has ",length(unique(net$moduleColors))," modules\n")
cat("The module size is\n")
print(table(net$moduleColors))
cat("The 1st gene's ID is ",names(net$moduleColors)[1]," and it belongs to module",
    net$moduleColors[1],"\n")
```

Or we can visually check their corresponding size and color as follows

```{r}
CoExpNets::plotModSizes(tissue="ad",which.one="CoExpROSMAP")
```



We can also access network properties by gene, e.g. the module membership of each gene within its module by using

```r
expr.data = getExprDataFromTissue(which.one="CoExpROSMAP",tissue="ad")
getMM(which.one="CoExpROSMAP",tissue="ad",
      genes=gsub("\\.[0-9]+$","",colnames(expr.data)[1:5]))
```

Or we can represent an histogram with all MM values for module yellow

```r
mms = getMM(which.one="CoExpROSMAP",tissue="ad",
                  genes=getGenesFromModule(which.one="CoExpROSMAP",
                                                 tissue="ad",
                                           module="yellow"))
hist(mms,main="MM values for black module in ROS/MAP AD subjects GCN")
```

We can also get the top 5 Hub genes for that very same module

```r
fromEnsembl2GeneName(names(mms[order(mms,decreasing=T)][1:5]))
```

We can also access gene information by modules. What are the modules available in the AD subjects GCN?

```r
getModulesFromTissue(which.one="CoExpROSMAP",tissue="ad")
```

Which are the genes within the tan module?

```r
genes = getGenesFromModule(which.one="CoExpROSMAP",
                           tissue="ad",
                                 module="tan")
cat("Tan module has",length(genes),"genes\n")
```

WGCNA has a nice feature and this is that one can visually inspect module similarities in terms of gene expression by using the eigengene. An eigengene is the 1st PCA component of expression and we can construct a dendrogram and see how the eigengenes are organised within such graph. Modules hanging from the same root will be more similar than others far away in the dendrogram.

```r
#We can get the eigengenes in the form of a matrix with samples in rows and columns for modules 
#There is one eigengene for each module
tissue = "ad"
which.one = "CoExpROSMAP"
egs = getNetworkEigengenes(which.one=which.one,tissue=tissue)
cat(tissue,"network was obtained from",nrow(egs),
    "samples and has",ncol(egs),"modules\n")
#And now we plot the EGs
plotEGClustering(which.one=which.one,tissue=tissue)
```

O we can simply plot its correlation

```r
library(corrplot)
corrplot(cor(egs),tl.srt=45, 
         tl.col="black",
         type="upper",
         tl.cex=0.45,
         order="hclust",
			  title=paste0("Eigengene correlations for ",tissue," within ",which.one),
			mar=c(0,3,3.5,2))
```

Note that the correlation plot will apply a clustering process too so the modules with highest correlation will appear together. There will be a correspondence between closeness at the dendrogramm and closeness at the triangular matrix (check for example thistle1, darkmagenta and greenyellow modules in both plots).

Finally, we can check the relation of modules with their covariates.
```r
library(WGCNA)
#We whole network has all samples and some modules will capture differences between
#genes in cases and controls, lets see the whole picture
covlist = colnames(getCovariates(tissue="allsamples",
                                 which.one="CoExpROSMAP"))
corWithCatTraits(which.one="CoExpROSMAP",tissue="allsamples",
                 covlist=covlist,covs=getCovariates(tissue="allsamples"))
```

So interesting module for understanding AD would be those that show any significant relation with assessed states of cognitive decline, e.g. darkmagenta, darkturquoise and royalblue. Note how PMI, age, batch, gender and race show no correlation whatsoever as they were applied for correction of the gene expression before creating the network.

##Making sense of the networks

Thanks to annotation of the modules, we can obtain information on what the networks reflect in terms of the biology of the genes. There are two types of annotations. But both of them work at the module level. Given a network compound of modules $m_1, m_2, ..., m_n$, the unit of annotation is the module. In this regard, each module $m_i$ will be annotated from enrichment analyses against manually curated ontologies like the Gene Ontology, KEGG and REACTOME, but also with cell type enrichment by using both in-house built cell markers and those coming from WGCNA.

Let us use, in this section, the ROS/MAP "all" network.

```r
net = getNetworkFromTissue(which.one="CoExpROSMAP",tissue="allsamples")
cat("Network built with",length(net$moduleColors),"genes and",nrow(net$MEs),"samples\n")
cat("The list of modules is",paste0(unique(net$moduleColors),collapse=", "),"\n")
```

Let us focus on the `royalblue` module. What is the enrichment of that module in terms of the biology of its genes?

```r
report=CoExpNets::reportOnModule(tissue="allsamples",
                      which.one="CoExpROSMAP",
                      module="turquoise",
                      ctcollapse=F,
                      how.many=30)
```

Now we can access the report's fields

```r
names(report)
```

This is a comprehensive report on a module in which we can access a lot of information that can tell us about gene set function 

```r
#How many genes do we have in the module
report$size
#What are the top five annotations for that module coming from the Gene Ontology
report$go.report
```

The module hs 359 genes, and the annotation DB based enrichment has 7 significant terms. The field `go.report` gives us the top five most significant annotation DB based enrichment terms that show, at a glimpse, what the module genes do. In this case, this module is clearly doing something related to DNA repair and chromatin (terms are ordered by significance). 

All in all, we can access all the DB based annotation separately as in

```r
report=functionalReportOnModule(tissue="allsamples",
                                which.one="CoExpROSMAP",
                                module="turquoise")
dim(report)
names(report)
table(report$domain)
```

As we see, there are 26 terms, 18 coming from the Gene Ontology Biological Process database, 4 from KEGG and 4 from REACTOME pathway databases. So for example, in order to get information about REACTOME and KEGG terms

```r
report[report$domain %in% c("rea","keg"),c("term.id","term.name","p.value")]
```


This data shows the corresponding REACTOME pathways that come up from the analysis, check [gProfileR](https://cran.r-project.org/web/packages/gProfileR/index.html) package documentation.

Finally, we can also generate a plot for the cell type

```r
results = cellTypeByModule(tissue="allsamples",which.one="CoExpROSMAP",do.plot=T)
```

