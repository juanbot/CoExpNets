
# CoExpNets brief instructions 

CoExpNets is based on the excellent and widely used WGNCA package, but incorporates a refinement we detail in the following paper <https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6>. It provides ways to create, access and exploit networks based on basic WGCNA features and ideas but we extend them a little bit and make it accessible through a Web interface.

An interesting part of this is that we make available a series of packages like CoExpROSMAP <https://github.com/juanbot/CoExpROSMAP> and CoExpGTEx <https://github.com/juanbot/CoExpGTEx> with already created networks you can use and exploit through the main package, CoExpNets.

## Install the development version from GitHub:

Simply do this
```r
devtools::install_github('juanbot/CoExpNets')
```

And that will be all. More help to come soon.
In the meantime, you can access the tutorials in the package.

# Available tutorials

* ![Tutorial for creation of CoExp networks](inst/tutorials/Tutorial_1.md)

* ![Tutorial for using CoExp networks suite](inst/tutorials/Tutorial_2.md)

* ![Tutorial for preparing expression data before creating the networks](inst/tutorials/Tutorial_3.md)


# Credits

The development of this suite of packages is leaded by Juan A. Botía. Co-expression network creation is leaded by Juan A. Botía and Mina Ryten from the Ryten Lab but many people contributed in some way including Jana Vandrovcova, Mar Matarin, Paola Forabosco, Conceisao Bettencourt, Seb Guelfi, Sonia García-Ruiz.

Sonia García-Ruíz is deploying all these co-expression networks at this Web page from the Ryten Lab

<https://snca.atica.um.es/coexp/Run/Catalog/>

which is currently under construction.

If you want to check examples of collaborations in which we used CoExpNets for annotation of genes, their function, cell type and main co-expressed genes, check the following PUBMED references

* <https://www.ncbi.nlm.nih.gov/pubmed/30328509>
* <https://www.ncbi.nlm.nih.gov/pubmed/29365066>
* <https://www.ncbi.nlm.nih.gov/pubmed/29127725>
* <https://www.ncbi.nlm.nih.gov/pubmed/28899015>
* <https://www.ncbi.nlm.nih.gov/pubmed/28575651>
* <https://www.ncbi.nlm.nih.gov/pubmed/26912063>
* <https://www.ncbi.nlm.nih.gov/pubmed/26707700>

And if we use the resource please cite us, with the GitHub URL and also this paper

<https://www.ncbi.nlm.nih.gov/pubmed/28403906>
