

# Preparing gene expression data for co-expression networks

# Selecting samples and genes

Let us suppose we want to work on the [ROS/MAP](https://www.synapse.org/#!Synapse:syn3219045) data. This is a very nice dataset with RNA-seq gene expression data (FPKM values) on post-morten Cortex samples obtained from two different cohorts of elderly people with very similar habits (i.e. they were monks and nuns so within each group of people they had the same life style).

We will see how we can correct for batch effect, gender, sex and PMI (post morten interval) but also for unknow covariates by using surrogate variable analysis (SVA).

We use the file `examples/ROSMAP_RNAseq_FPKM_gene.tsv` as a starting point. This is the gene expression profiling of 640 RNA-seq cortex samples with 55889 genes profiled there. We will only focus on the most expressed genes across those samples. So the first thing to do so you can have the expression data ready is to filter and keep only those genes expressed with more than 2 RPKM values across at least 80% of the samples. Note this is something already done for you. Therefore, you will find `examples/rosmap_rpkms.rds` ready for you.


```r
rpkms = read.table(paste0(gdp.rosmap(),"/ROSMAP_RNAseq_FPKM_gene.tsv"),
                       header=T,row.names=1)
colnames(rpkms) = gsub("X","",colnames(rpkms))
rpkms = rpkms[,-1]
expressed.genes = rowSums(rpkms > 2)  > (0.8 * ncol(rpkms))
rpkms = rpkms[expressed.genes,]
```

As a result of this, we will get 640 samples and 11304 genes in total.
Normally, we would apply a lower cut-off for expression. However, we want to keep it simple in this tutorial so we are more stringent to get an easier to work with dataset.

```r
covs = getCovariates(tissue="all",which.one="rosmap")
common = intersect(rownames(covs),colnames(rpkms))
covs = covs[common,]
rpkms = rpkms[,common]
stopifnot(identical(colnames(rpkms),rownames(covs)))
saveRDS(rpkms,"~/tmp/rosmap_rpkms.rds")
```


#Visualizing the samples

Now the data is ready to work with. The rows of covariates (samples) are aligned with the columns of the expression data (the same samples). Before doing any other thing, we would like to explore ROS/MAP samples at the level of the covariates in order to check for evident batch effects. We won't be using all samples for this representation (it takes a while to compute), only 100 randomly chosen. We will be using the nice method `limma::plotMDS` which takes gene expression data and plots each sample in a 2D plot, tring to mimic in the distance between samples at the plot, their log2 fold changes in expression. Which means that samples appearing together are similar in expression, and the other way around.

```r
CoExpROSMAP::initDb()
#We need limma for plotMDS
library(limma)
#We also need swamp for princomp
library(swamp)

#Load the selected genes and samples from the file produced with the example code from above
rpkms = readRDS("~/tmp/rosmap_rpkms.rds")
#Get all covariates of interest for ROS/MAP project
covs = CoExpROSMAP::getCovariates(tissue="allsamples")
#Keep only samples with available covariates
common = intersect(rownames(covs),colnames(rpkms))
covs = covs[common,]
rpkms = rpkms[,common]
stopifnot(identical(colnames(rpkms),rownames(covs)))
#We will use only batch for getting the MDS plot coloured per sample
covariate = "batch"
#MDS takes a while to compute so we will use only 100 samples to represent instead of the total
#(more than 600) ROS/MAP has
mask = sample(1:ncol(rpkms),100)
#Get one color for each batch batch category
colors = rainbow(length(unique(as.numeric(covs[,covariate]))))
finalcolors = colors[as.numeric(covs[,covariate])][mask]
#Plot samples with IDs and colour with their corresponding batch
#samples will be disposed in 2D such that the distance between samples reflects differences
#in log2 expression
plotMDS(rpkms[,mask],col=finalcolors,
				main=paste0("ROS/MAP MDS using ",covariate))
legend("topright",fill=colors,
				legend=levels(covs[,covariate]))
```

And we can see that batch effect is evident across samples because samples with the same color (i.e. they come from the same batch) are usually appearing together.  We can also check possible Gender biases, with a similar code.

```r
covariate = "gender"
mask = sample(1:ncol(rpkms),100)
colors = rainbow(length(unique(as.numeric(covs[,covariate]))))
finalcolors = colors[as.numeric(covs[,covariate])][mask]
plotMDS(rpkms[,mask],col=finalcolors,
				main=paste0("ROS/MAP MDS using ",covariate))
legend("topright",fill=colors,
				legend=levels(covs[,covariate]))
```

And the situation is more subtle (it usually is) but can be appreciated in some of the samples.


#Detecting possible outliers

We will see later how to correct this in the data so we clean this effect. But before doing this we will check for the existence of outliers. We do it by using a simple test based on the comparison of each sample z-score (outlier minus the mean, divided by sd) with the sample distribution. Check `help(grubbs.test)`.


```r
library(outliers)
library(flashClust)
outlist = NULL
#Traspose data frame and get a distance matrix between pairs of samples
sampdist = as.matrix(dist(CoExpNets::trasposeDataFrame(rpkms,F)))
result = grubbs.test(apply(sampdist,1,sum))
#Is the test significant at 95% confidence level?
if(result$p.value < 0.05){
  cat("We should drop out",names(result$p.value),"\n")
  plot(flashClust(dist(CoExpNets::trasposeDataFrame(rpkms,F)), method = "average"),cex=0.5,
       main=paste0("Detected outlier ",names(result$p.value)," in ROS/MAP cases + ctrls"))
}else
  cat("No plausible outlier detected\n")


```


As a visual test of what the statistical test is saying, we can represent each sample in a dendrogram through hierarchical clustering (we know how to do that) and we will see that effectively this sample appears as an extreme. However, one may doubt about it being an outlier. It will depend on how happy we are with the number of samples we have. We will actually keep it.


#Data normalisation

And now below we prepare the data again for normalization of gene expression across samples, by using quantile normalization and save results. We will need samples at columns and genes at rows. Quantile normalization is based on the notion of the q-q plot in which two samples are identical if their q-q plot is a straight diagonal line. The paper is available at the material for these seminars.

```r
dim(rpkms)
rpkms = readRDS("~/tmp/rosmap_rpkms.rds")
dnames = dimnames(rpkms)

#Select five samples randomly
mask = sample(1:length(dnames[[2]]),5)
#Get five different colors to distinguish the plots
colors = rainbow(5)
#Generate the plot with the gene expression distribution for the
#1st sample
plot(density(as.numeric(rpkms[,mask[1]])),main="Before normalising",
     xlim=c(0,800),colors[1])
#Plot the rest of samples. Now you can see the values are not normalised
#across samples because the distribution plots are different
ret = lapply(mask[2:length(mask)],function(x){
  i = which(mask == x)
  lines(density(as.numeric(rpkms[,x])),col=colors[i])
  })

#Normalize samples
rpkms_qn = preprocessCore::normalize.quantiles(as.matrix(rpkms))
dimnames(rpkms_qn) = dnames
#Plot exactly the same. Now you will see plots are identical
plot(density(as.numeric(rpkms_qn[,mask[1]])),main="After normalising",xlim=c(0,800),colors[1])
ret = lapply(mask[2:length(mask)],function(x){
  i = which(mask == x)
  lines(density(as.numeric(rpkms_qn[,x])),col=colors[i])
  })
```

Below, we prepare the code for ComBat <https://academic.oup.com/biostatistics/article/8/1/118/252073>. ComBAT is based on Bayes and it is robust when the batch size is small.
We will use `covs$batch` as the batch covariate for samples. By using `mod=model.matrix(~1,data=covs)` we correct for the batch effect taking into account the rest of covariates so ComBAT can keep their potential effect on the data.

```r
rpkms_combat = sva::ComBat(dat=rpkms_qn,batch=covs$batch,
                       mod=model.matrix(~1,data=covs[,c("gender","age","pmi")]),
                      par.prior=T)

covariate="batch"
mask = sample(1:ncol(rpkms_combat),100)
colors = rainbow(length(unique(as.numeric(covs[,covariate]))))
finalcolors = colors[as.numeric(covs[,covariate])][mask]
plotMDS(rpkms_combat[,mask],col=finalcolors,
				main=paste0("ROS/MAP (batch corrected) MDS using ",covariate))
legend("topleft",fill=colors,
				legend=levels(covs[,covariate]))

```

And clearly, now the batch effect is not so apparent. Another way of checking whether ComBat has worked fine is to generate a list of PCAs and correlate them with covariates.

What follows is an adaptation of plots from the `swamp` R package to visually inspect the effect of ComBat correction. If we compare the two plots, we will see that the one created with `rpkms_combat` does not show any principal component correlated with batch. Uncorrected data, 'rpkms_qn' shows abundant correlation signals across most of the axes.

```r
pcres = prince(as.matrix(rpkms_qn),covs[,colnames(covs) != "diseasestatus"],top=20)
CoExpNets::princePlot(prince=pcres,main="ROS/MAP All samples, batch uncorrected FPKMQN")

pcres = prince(as.matrix(rpkms_combat),covs[,colnames(covs) != "diseasestatus"],top=20)
CoExpNets::princePlot(prince=pcres,main="ROS/MAP All samples, batch corrected FPKMQN")


```


As we see in the first heat-map, batch has a strong effect on the data that we have to hopefully correct with ComBat. In the plot that follows the 1st one, we see those correlations are not present anymoare. Still we have correlations with gender, age and PMI it would be desiderable to remove. Note that variables like cognitive decline (cogdx) shouldn't be corrected because these represent effects we want to study and they are not genuine features on the samples but assessed by an expert during the lifetime of the corresponding subject. As we still need to apply corrections, we will residual correct the data.

But before doing that, we will apply SVA analysis to try and find hidden effects in the data so we get a set of residuals which are clean enough to get a good co-expression network.

```r
#We define the model we want to test by
mm = model.matrix(~ gender + pmi + age + race,data=covs)
nullmm = model.matrix(~ 1,data=covs)
cat("Launching svaseq\n")
print(Sys.time())
svas = sva::svaseq(dat=as.matrix(rpkms_combat),mod=mm,mod0=nullmm)
print(Sys.time())

```

As we see, the `svaseq` method needs a null model in `mod0` and the model we want to test, `mod`. The null model is a trivial model which assumes there is no significant effect of any biological covariate in the data. The model we want to test includes `gender`, `pmi`, `age` and `race` in the data. SVA uses both model to look for additional covariates, the unknown ones, which affect all genes in a sample.
After the `standard` 5 iterations of the method, we can see `svaseq` generates 8 additional axes to be added to the regression model to apply to the batch corrected data. See more information on SVA plust ComBat at [the SVA paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3307112/).


#Correcting the data by a linear regression approach

Now, in `svas` we have the axes ready to be used for data correction. Being these axes surrogate variables or hidden effects we can't account for explicitly because there is no covariate defining them.

But before doing that, it would be important to check whether the discovered axes are not modeling some phenomena of interest in the data. If this was the case and we use the corresponding axis to correct the data, such phenomena of interest would be removed out of the data. We don't want that to happen. In our case, the variable of interest would be the one referring to cognitive decline, `cogdx`. What we do is to create a heat-map with correlations between axes and covariates and we're interested on the significance of correlation, not the magnitude.

```r
#Now we control whether the surrogate variables
#have some relation with the covariates used in the mm model
#through a heat map of correlation p-values
numeric.covs = covs[,colnames(covs) != "diseasestatus"]
#We will create a matrix for the heatmap putting at rows the
#covariates we know and at columns the new SVAs we just discovered
#In that way we can show with a heat map the correlation of all SVAs
#against all covariates
linp = matrix(ncol=svas$n.sv,nrow=ncol(numeric.covs))
rownames(linp) = colnames(numeric.covs)
colnames(linp) = paste0("SV",1:svas$n.sv)
linp[] = 0
for(cov in 1:ncol(numeric.covs)){
  for(sva in 1:svas$n.sv){
    if(svas$n.sv == 1)
      axis = svas$sv
    else
      axis = svas$sv[,sva]

    linp[cov,sva] = cor.test(as.numeric(numeric.covs[,cov]),axis)$p.value
  }
}
#We transform the p-values into a log10(pvalue) scale
#And we put a limit of -10 for the highest significance value
smallest = -10
linp10 <- log10(linp)
linp10 <- replace(linp10, linp10 <= smallest, smallest)
tonote <- signif(linp, 1)
gplots::heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none",
          trace = "none", symbreaks = F, symkey = F, breaks = seq(-20, 0, length.out = 100),
          key = T, colsep = NULL, rowsep = NULL, sepcolor = "black",
          sepwidth = c(0.05, 0.05), main = "Corr. of SVs and covs., Combat-FPKM QN",
          labCol = colnames(linp10),
          labRow = colnames(covs),
          xlab = "Surrogate variables")
```


We could remove SV4 from the data correction process as it clearly shows some significant correlation with the phenomena we want to study.

On the other hand, remember that in `rpkms_combat` we have the successfully batch corrected gene expression data with genes at the rows. In order to correct the data, we will have to act at gene level. This means that we will apply a regression model to each gene, modeling gene expression as an additive model of biological covariates (age, gender, pmi, race) plus the SVA axes we want to apply, i.e. all but 4. So, for the ith gene we will have the following model
$$y_i = w_1 age + w_2 gender + w_3 pmi + w_4 race + w_5 sva_1 + w_6 sva_2 + w_7 sva_3 +w_8 sva_5+ \ldots + w_{11} sva_8 + w_{12} + b_i,$$
where $y_i$ is gene expression for the ith gene and we represent cleaned data with $w_{12}$. Thus, applying a linear regression for each gene, we will keep only the residuals, i.e. differences between left and right part of the expression, as in
$$w_{12} = y_i - (w_1 age + w_2 gender + w_3 pmi + w_4 race + w_5 sva_1 + w_6 sva_2 + w_7 sva_3 +w_8 sva_5+ \ldots + w_{11} sva_8 + b_i)$$

So we use `apply` at each row of the gene expression data variable with a linear model where $y$ correspond to a vector or gene expression values for the corresponding gene, as the dependent variable. Independent variables would be `age`, `gender`, `age`, `race` and the 7 SVA axes.

Here is the code, simple, isn't it?

```r
covs.rs = covs[,match(c("gender","pmi","age","race"),colnames(covs))]
print("Getting residuals")
resids <- apply(rpkms_combat, 1, function(y){
  lm( y ~ . , data=cbind(covs.rs,svas$sv[,-4]))$residuals
})
```

Now we have to finally check that the residuals are clean of any effect.

```r
pcres = prince(as.matrix(CoExpNets::trasposeDataFrame(resids,F)),covs.rs,top=20)
CoExpNets::princePlot(prince=pcres,
                 main="ROS/MAP: residual cor. with gender, pmi, age, race")

```

Clearly clean... Last thing is storing the data. Now we are ready to generate a co-expression network.

```r
saveRDS(resids,"~/tmp/rosmap_residuals.rds")
```


