# [preciseTAD](https://dozmorovlab.github.io/preciseTAD/)

<!-- [![Travis build
status](https://travis-ci.com/stilianoudakis/preciseTAD.svg?branch=master)](https://travis-ci.com/stilianoudakis/preciseTAD) -->

<!-- badges: start -->

[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

## Overview

preciseTAD provides functions to predict the location of boundaries of topologically associated domains (TADs) and chromatin loops at base-level resolution. As an input, it takes BED-formatted genomic coordinates of domain boundaries detected from low-resolution Hi-C data, and coordinates of high-resolution genomic annotations from ENCODE or other consortia. preciseTAD employs several feature engineering strategies and resampling techniques to address class imbalance, and trains an optimized random forest model for predicting low-resolution domain boundaries. Translated on a base-level, preciseTAD predicts the probability for each base to be a boundary. Density-based clustering and scalable partitioning techniques are used to detect precise boundary regions and summit points. Compared with low-resolution boundaries, preciseTAD boundaries are highly enriched for CTCF, RAD21, SMC3, and ZNF143 signal and more conserved across cell lines. The pre-trained model can accurately predict boundaries in another cell line using CTCF, RAD21, SMC3, and ZNF143 annotation data for this cell line. 

The main functions (in order of implementation) are:

- `extractBoundaries()` accepts a 3-column data.frame or matrix with the chromosomal coordinates of user-defined TADs and outputs the unique boundaries
- `bedToGRangesList()` accepts a filepath containing BED files representing the coordinates of ChIP-seq defined functional genomic annotations
- `createTADdata()` accepts a set of unique boundaries and genomic annotations derived from `preciseTAD::extractBoundaries()` and `preciseTAD::bedToGRangesList()` respectively to create the data matrix used to build a model to predict TAD boundary regions
- `TADrandomForest()` a wrapper of the `randomForest` package which implements a random forest binary classification algorithm on TAD boundary data
- `preciseTAD()` which leverages a TAD boundary prediction model (i.e., random forest) and density-based clustering to predict TAD boundary coordinates at a base-level resolution

## Installation

The latest version of `preciseTAD` can be directly installed from Github:

``` r
devtools::install_github('dozmorovlab/preciseTAD', build_vignettes = TRUE)
# Or, install the developmental version
# install_github('stilianoudakis/preciseTAD', build_vignettes = TRUE)
library(preciseTAD)
```

Alternatively, the package can be installed from Bioconductor (to be submitted):

``` r
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("preciseTAD")
# library(preciseTAD)
```

## Usage

Below is a brief workflow of how to implement `preciseTAD` on binned data from CHR1 to get precise base pair coordinates of TAD boundaries for a 10mb section of CHR 22. For more details, including the example how to use the pre-trained model, see `vignette("preciseTAD")`

First, you need to obtain called TAD boundaries using an established TAD-caller. As an example, consider the [Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead) TAD-caller, a part of the juicer suite of tools developed by the Aiden Lab. Arrowhead outputs a .txt file with the chromosomal start and end coordinates of their called TADs. As an example, we have provided Arrowhead TADs for GM12878 at 5kb resolution.

``` r
data("arrowhead_gm12878_5kb")
head(arrowhead_gm12878_5kb)
```

The unique boundaries for CHR1 and CHR22 can be extracted as:

``` r
bounds <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb, preprocess=FALSE, CHR=c("CHR1","CHR22"), resolution=5000)
```

Next, you will need to download cell line-specific ChIP-seq data in the form of BED files from [ENCODE](https://www.encodeproject.org/chip-seq-matrix/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&assay_title=TF%20ChIP-seq&status=released). Once, you have downloaded your preferred list of functional genomic annotations, store them in a specific file location. These files can then be converted into a GRangesList object and used for downstream modeling using the following command:

``` r
path = "pathToBEDfiles"
tfbsList <- bedToGRangesList(filepath=path, bedList=NULL, bedNames=NULL, pattern = "*.bed", signal=4)
```

As an example, we have already provided a GRangesList object with a variety of transcription factor binding sites specific to the GM12878 cell line. Once you load it in, you can see the list of transcription factors using the following:

``` r
data("tfbsList")
names(tfbsList)
```

For the purposes of this example, let's focus only on CTCF, RAD21, SMC3, and ZNF143 transcription factors.

``` r
tfbsList_filt <- tfbsList[names(tfbsList) %in% c("Gm12878-Ctcf-Broad","Gm12878-Rad21-Haib","Gm12878-Smc3-Sydh","Gm12878-Znf143-Sydh")]
```

Now, using the “ground-truth” boundaries and the following TFBS, we can build the data matrix that will be used for predictive modeling. The following command creates the training data from CHR1 and reserves the testing data from CHR22. We specify 5kb sized genomic bins (to match the resolution used to call the original TADs), a distance-type feature space, and apply random under-sampling (RUS) on the training data only.

``` r
set.seed(123)
tadData <- createTADdata(bounds.GR          = bounds,
                         resolution         = 5000,
                         genomicElements.GR = tfbsList_filt,
                         featureType        = "distance",
                         resampling         = "rus",
                         trainCHR           = "CHR1",
                         predictCHR         = "CHR22")
```

We can now implement our machine learning algorithm of choice to predict TAD-boundary regions. Here, we opt for the random forest algorithm.

``` r
set.seed(123)
tadModel <- TADrandomForest(trainData  = tadData[[1]],
                            testData   = tadData[[2]],
                            tuneParams = list(mtry=2,
                                              ntree=500,
                                              nodesize=1),
                            cvFolds    = 3,
                            cvMetric   = "Accuracy",
                            verbose    = TRUE,
                            model      = TRUE,
                            importances= TRUE,
                            impMeasure = "MDA",
                            performances = TRUE)
# The model itself
tadModel[[1]]
                            
# Variable importances (mean decrease in accuracy)
tadModel[[2]]

# Model performance metrics
tadModel[[3]]
```

Lastly, we take our TAD-boundary region predictive model and use it to make predictions on a 10mb section of CHR22:35,000,000-45,000,000.

``` r
# Run preciseTAD
set.seed(123)
pt <- preciseTAD(genomicElements.GR = tfbsList_filt,
                 featureType        = "distance",
                 CHR                = "CHR22",
                 chromCoords        = list(35000000,45000000),
                 tadModel           = tadModel[[1]],
                 threshold          = 1.0,
                 flank              = NULL,
                 verbose            = TRUE,
                 parallel           = 2,
                 DBSCAN_params      = list(5000,3))
                 
# View preciseTAD predicted boundary coordinates between CHR22:35mb-45mb
pt[[2]]
```
