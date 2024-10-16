# pgxRpi

Welcome to our R wrapper package for Progenetix REST API that leverages the capabilities of [Beacon v2](https://docs.genomebeacons.org/) specification. Please note that a stable internet connection is required for the query functionality. This package is aimed to simplify the process of accessing oncogenomic data from [Progenetix](https://progenetix.org/) database via the Beacon v2 API with some extensions (BeaconPlus). 

You can install this package using either of the following methods:

### From Bioconductor

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("pgxRpi")
```

### From Github for the latest development version 

```r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("progenetix/pgxRpi")
```

For accessing metadata of biosamples/individuals, or learning more about filters, get started from the vignette [Introduction_1_loadmetadata](https://bioconductor.org/packages/devel/bioc/vignettes/pgxRpi/inst/doc/Introduction_1_loadmetadata.html).

For accessing CNV variant data, get started from this vignette [Introduction_2_loadvariants](https://bioconductor.org/packages/devel/bioc/vignettes/pgxRpi/inst/doc/Introduction_2_loadvariants.html).

For accessing CNV frequency data, get started from this vignette [Introduction_3_loadfrequency](https://bioconductor.org/packages/devel/bioc/vignettes/pgxRpi/inst/doc/Introduction_3_loadfrequency.html).

For processing local pgxseg files, get started from this vignette [Introduction_4_process_pgxseg](https://bioconductor.org/packages/devel/bioc/vignettes/pgxRpi/inst/doc/Introduction_4_process_pgxseg.html).

If you encounter problems, try to reinstall the latest version. If reinstallation doesn't help, please contact us.
