
<!-- README.md is generated from README.Rmd. Please edit that file -->

# preciseTAD: A transfer learning framework for 3D domain boundary prediction at base-pair resolution

<!-- [![Travis build
status](https://travis-ci.com/stilianoudakis/preciseTAD.svg?branch=master)](https://travis-ci.com/stilianoudakis/preciseTAD) -->
<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

Stilianoudakis, Spiro, and Mikhail G. Dozmorov. “preciseTAD: A machine
learning framework for precise 3D domain boundary prediction at
base-level resolution.” bioRxiv (2020).
<https://doi.org/10.1101/2020.09.03.282186>

Predicted preciseTAD boundary points (PTBPs) and regions (PTBRs) for 60
cell lines are available
[here](https://drive.google.com/drive/folders/15Rc6PhrrBjThwE-5dSyNX-ILELaUu6uG?usp=sharing).

## Overview

preciseTAD provides functions to predict the location of boundaries of
topologically associated domains (TADs) and chromatin loops at
base-level resolution. As an input, it takes BED-formatted genomic
coordinates of domain boundaries detected from low-resolution Hi-C data,
and coordinates of high-resolution genomic annotations from ENCODE or
other consortia. preciseTAD employs several feature engineering
strategies and resampling techniques to address class imbalance, and
trains an optimized random forest model for predicting low-resolution
domain boundaries. Translated on a base-level, preciseTAD predicts the
probability for each base to be a boundary. Density-based clustering and
scalable partitioning techniques are used to detect precise boundary
regions and summit points. Compared with low-resolution boundaries,
preciseTAD boundaries are highly enriched for CTCF, RAD21, SMC3, and
ZNF143 signal and more conserved across cell lines. The pre-trained
model can accurately predict boundaries in another cell line using CTCF,
RAD21, SMC3, and ZNF143 annotation data for this cell line.

The main functions (in order of implementation) are:

-   `extractBoundaries()` accepts a 3-column data.frame or matrix with
    the chromosomal coordinates of user-defined domains and outputs the
    unique boundaries. The second and third columns are the domain
    anchor centers.
-   `bedToGRangesList()` accepts a filepath containing BED files
    representing the coordinates of ChIP-seq defined functional genomic
    annotations
-   `createTADdata()` accepts a set of unique boundaries and genomic
    annotations derived from `extractBoundaries()` and
    `bedToGRangesList()`, respectively, to create the data matrix used
    to build a model to predict domain boundary regions
-   `TADrandomForest()` a wrapper of the `randomForest` package which
    implements a random forest binary classification algorithm on domain
    boundary data
-   `preciseTAD()` which leverages a domain boundary prediction model
    (i.e., random forest) and density-based clustering to predict TAD
    boundary coordinates at a base-level resolution

## Installation

`preciseTAD` can be installed from Bioconductor:

``` r
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("preciseTAD")
library(preciseTAD)
#> 
```

The latest version of `preciseTAD` can be directly installed from
Github:

``` r
devtools::install_github("dozmorovlab/preciseTAD", build_vignettes = TRUE)
library(preciseTAD)
```

## Usage

Below is a brief workflow of how to implement `preciseTAD` on binned
data from CHR1 to get precise base pair coordinates of TAD boundaries
for a 10mb section of CHR 22. For more details, including the example
how to use the pre-trained model, see `vignette("preciseTAD")`

First, you need to obtain called TAD boundaries using an established
TAD-caller. As an example, consider the
[Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead)
TAD-caller, a part of the juicer suite of tools developed by the Aiden
Lab. Arrowhead outputs a .txt file with the chromosomal start and end
coordinates of their called TADs. As an example, we have provided
Arrowhead TADs for GM12878 at 5kb resolution.

``` r
data("arrowhead_gm12878_5kb")
head(arrowhead_gm12878_5kb)
#>   V1        V2        V3
#> 1  1  49375000  50805000
#> 2  1  16830000  17230000
#> 3  1 163355000 164860000
#> 4  1 231935000 233400000
#> 5  1 149035000 149430000
#> 6  1   3995000   5505000
```

The unique boundaries for CHR1 and CHR22 can be extracted as:

``` r
bounds <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb, filter = FALSE, CHR = c("CHR1", "CHR22"), resolution = 5000)
bounds
#> GRanges object with 1901 ranges and 0 metadata columns:
#>         seqnames            ranges strand
#>            <Rle>         <IRanges>  <Rle>
#>     787     chr1     815000-815001      *
#>    9196     chr1     890000-890001      *
#>      37     chr1     915000-915001      *
#>    8923     chr1   1005000-1005001      *
#>     275     chr1   1015000-1015001      *
#>     ...      ...               ...    ...
#>    8379    chr22 50815000-50815001      *
#>   16788    chr22 50935000-50935001      *
#>    8382    chr22 50945000-50945001      *
#>   16791    chr22 51060000-51060001      *
#>   16670    chr22 51235000-51235001      *
#>   -------
#>   seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Next, you will need to download cell line-specific ChIP-seq data in the
form of BED files from
[ENCODE](https://www.encodeproject.org/chip-seq-matrix/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&assay_title=TF%20ChIP-seq&status=released).
Once, you have downloaded your preferred list of functional genomic
annotations, store them in a specific file location. These files can
then be converted into a GRangesList object and used for downstream
modeling using the following command:

``` r
path <- "pathToBEDfiles"
tfbsList <- bedToGRangesList(filepath = path, bedList = NULL, bedNames = NULL, pattern = "*.bed", signal = 4)
```

As an example, we have already provided a GRangesList object with a
variety of transcription factor binding sites specific to the GM12878
cell line. Once you load it in, you can see the list of transcription
factors using the following:

``` r
data("tfbsList")
names(tfbsList)
#>  [1] "Gm12878-Atf3-Haib"       "Gm12878-Cfos-Sydh"      
#>  [3] "Gm12878-Cmyc-Uta"        "Gm12878-Ctcf-Broad"     
#>  [5] "Gm12878-Egr1-Haib"       "Gm12878-Ets1-Haib"      
#>  [7] "Gm12878-Gabp-Haib"       "Gm12878-Jund-Sydh"      
#>  [9] "Gm12878-Max-Sydh"        "Gm12878-Mazab85725-Sydh"
#> [11] "Gm12878-Mef2a-Haib"      "Gm12878-Mxi1-Sydh"      
#> [13] "Gm12878-P300-Sydh"       "Gm12878-Pol2-Haib"      
#> [15] "Gm12878-Pu1-Haib"        "Gm12878-Rad21-Haib"     
#> [17] "Gm12878-Rfx5-Sydh"       "Gm12878-Sin3a-Sydh"     
#> [19] "Gm12878-Six5-Haib"       "Gm12878-Smc3-Sydh"      
#> [21] "Gm12878-Sp1-Haib"        "Gm12878-Srf-Haib"       
#> [23] "Gm12878-Taf1-Haib"       "Gm12878-Tr4-Sydh"       
#> [25] "Gm12878-Yy1-Sydh"        "Gm12878-Znf143-Sydh"
```

For the purposes of this example, let’s focus only on CTCF, RAD21, SMC3,
and ZNF143 transcription factors.

``` r
tfbsList_filt <- tfbsList[names(tfbsList) %in% c("Gm12878-Ctcf-Broad", "Gm12878-Rad21-Haib", "Gm12878-Smc3-Sydh", "Gm12878-Znf143-Sydh")]
```

Now, using the “ground-truth” boundaries and the following TFBS, we can
build the data matrix that will be used for predictive modeling. The
following command creates the training data from CHR1 and reserves the
testing data from CHR22. We specify 5kb sized genomic bins (to match the
resolution used to call the original TADs), a distance-type feature
space, and apply random under-sampling (RUS) on the training data only.

``` r
set.seed(123)
tadData <- createTADdata(bounds.GR          = bounds,
                         resolution         = 5000,
                         genomicElements.GR = tfbsList_filt,
                         featureType        = "distance",
                         resampling         = "rus",
                         trainCHR           = "CHR1",
                         predictCHR         = "CHR22"
)
```

We can now implement our machine learning algorithm of choice to predict
TAD-boundary regions. Here, we opt for the random forest algorithm.

``` r
set.seed(123)
tadModel <- TADrandomForest(trainData    = tadData[[1]],
                            testData     = tadData[[2]],
                            tuneParams   = list(mtry = 2, ntree = 500, nodesize = 1),
                            cvFolds      = 3,
                            cvMetric     = "Accuracy",
                            verbose      = TRUE,
                            model        = TRUE,
                            importances  = TRUE,
                            impMeasure   = "MDA",
                            performances = TRUE)
#> Loading required package: lattice
#> Loading required package: ggplot2
#> + Fold1: mtry=2, ntree=500, nodesize=1 
#> - Fold1: mtry=2, ntree=500, nodesize=1 
#> + Fold2: mtry=2, ntree=500, nodesize=1 
#> - Fold2: mtry=2, ntree=500, nodesize=1 
#> + Fold3: mtry=2, ntree=500, nodesize=1 
#> - Fold3: mtry=2, ntree=500, nodesize=1 
#> Aggregating results
#> Fitting final model on full training set
# The model itself
tadModel[[1]]
#> Random Forest 
#> 
#> 3190 samples
#>    4 predictor
#>    2 classes: 'No', 'Yes' 
#> 
#> No pre-processing
#> Resampling: Cross-Validated (3 fold) 
#> Summary of sample sizes: 2126, 2127, 2127 
#> Resampling results:
#> 
#>   MCC        ROC        Sens    Spec       Pos Pred Value  Neg Pred Value
#>   0.4679235  0.7990197  0.7072  0.7598457  0.747661        0.7224703     
#>   Accuracy   Kappa    
#>   0.7335413  0.4670627
#> 
#> Tuning parameter 'mtry' was held constant at a value of 2
#> Tuning
#>  parameter 'ntree' was held constant at a value of 500
#> Tuning
#>  parameter 'nodesize' was held constant at a value of 1

# Variable importances (mean decrease in accuracy)
tadModel[[2]]
#>                 Feature Importance
#> 4 `Gm12878-Znf143-Sydh`   75.98200
#> 3   `Gm12878-Smc3-Sydh`   62.50681
#> 2  `Gm12878-Rad21-Haib`   62.21715
#> 1  `Gm12878-Ctcf-Broad`   32.61855

# Model performance metrics
tadModel[[3]]
#>              Metric  Performance
#> 1                TN 6.400000e+03
#> 2                FN 6.700000e+01
#> 3                FP 2.954000e+03
#> 4                TP 2.390000e+02
#> 5             Total 9.660000e+03
#> 6       Sensitivity 7.810458e-01
#> 7       Specificity 6.841993e-01
#> 8             Kappa 8.363202e-02
#> 9          Accuracy 6.872671e-01
#> 10 BalancedAccuracy 7.326225e-01
#> 11        Precision 7.485124e-02
#> 12              FPR 3.158007e-01
#> 13              FNR 1.036029e-02
#> 14              NPV 9.896397e-01
#> 15              MCC 1.732169e-01
#> 16               F1 1.366105e-01
#> 17              AUC 7.983352e-01
#> 18           Youden 4.652450e-01
#> 19            AUPRC 9.383824e-02
```

Lastly, we take our TAD-boundary region predictive model and use it to
make predictions on a 10mb section of CHR22:35,000,000-45,000,000.

``` r
# Run preciseTAD
set.seed(123)
pt <- preciseTAD( genomicElements.GR = tfbsList_filt,
                  featureType        = "distance",
                  CHR                = "CHR22",
                  chromCoords        = list(35000000, 45000000),
                  tadModel           = tadModel[[1]],
                  threshold          = 1.0,
                  verbose            = FALSE,
                  parallel           = 2,
                  DBSCAN_params      = list(30000, 3),
                  slope              = 5000,
                  genome             = "hg19")

# View preciseTAD predicted boundary coordinates between CHR22:35mb-45mb
pt[[2]]
#> GRanges object with 64 ranges and 0 metadata columns:
#>        seqnames            ranges strand
#>           <Rle>         <IRanges>  <Rle>
#>    [1]    chr22 35337875-35342203      *
#>    [2]    chr22 35413232-35426376      *
#>    [3]    chr22 35543040-35549073      *
#>    [4]    chr22 35620677-35634681      *
#>    [5]    chr22 35740783-35747642      *
#>    ...      ...               ...    ...
#>   [60]    chr22 43261091-43270096      *
#>   [61]    chr22 43425520-43490629      *
#>   [62]    chr22 43774902-43791526      *
#>   [63]    chr22 44272171-44283944      *
#>   [64]    chr22 44382230-44396240      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
