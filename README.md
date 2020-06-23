
<!-- README.md is generated from README.Rmd. Please edit that file -->

# preciseTAD

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/stilianoudakis/preciseTAD.svg?branch=master)](https://travis-ci.com/stilianoudakis/preciseTAD)
<!-- badges: end -->

`preciseTAD` is a collection of functions that transforms topologically
associating domain(TAD)-calling into a supervised learning framework for
more precise boundary prediction. To do so, a binary classifier is built
on user-defined “ground truth” TAD boundaries and ChIP-seq peak
functional genomic annotations. `preciseTAD` provides functionality for
creating the response vector and feature space necessary for building
such a model. `preciseTAD` is random forest based and leverages the
ensemble classifier using density-based clustering (DBSCAN) and scalable
partitioning around medoids (CLARA) on predicted probabilities of base
pair coordinates to more precisely predict TAD boundaries at base pair
resolution.

The main functions (in order of implementation) are:

  - `extractBoundaries` accepts a 3-column data.frame or matrix with the
    chromosomal coordinates of user-defined TADs and outputs the unique
    boundaries
  - `bedToGRangesList` accepts a filepath containing BED files
    representing the coordinates of ChIP-seq defined functional genomic
    annotations
  - `createTADdata` accepts a set of unique boundaries and genomic
    annotations derived from `preciseTAD::extractBoundaries` and
    `preciseTAD::bedToGRangesList` respectively to create the data
    matrix used to built a model to predict TAD boundary regions
  - `TADrandomForest` a wrapper of the `randomForest` package which
    implements a random forest binary classification algorithm on TAD
    boundary data
  - `preciseTAD` which leverages a TAD boundary prediction model
    (i.e. random forest) and density-based clustering to predict TAD
    boundary coordinates at base pair resolution

## Installation

First make sure you have all dependencies installed in R.

    install.packages(c('pbapply', 
                       'doSNOW', 
                       'cluster',
                       'foreach',
                       'parallel',
                       'bigmemory',
                       'IRanges',
                       'GenomicRanges',
                       'dbscan'))

The latest version of `preciseTAD` can be directly installed from
Github:

    devtools::install_github('stilianoudakis/preciseTAD', build_vignettes = TRUE)
    library(preciseTAD)

Alternatively, the package can be installed from Bioconductor (to be
submitted):

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("preciseTAD")
    library(preciseTAD)

## Usage

Below is a brief workflow of how to implement `preciseTAD` on binned
data from CHR1 to get precise base pair coordinates of TAD boundaries
for a 10mb section of CHR 22. For more details visit the `preciseTAD`
vignette.

First you need to obtain called TAD boundaries using an established
TAD-caller. As an example, consider the
[ARROWHEAD](https://github.com/aidenlab/juicer/wiki/Arrowhead)
TAD-caller, a part of the juicer suite of tools developed by the Aiden
Lab. ARROWHEAD outputs a .txt file with the chromosomal start and end
coordinates of their called TADs. As an example, we have provided
ARROWHEAD TADs for GM12878 at 5 resolution. Simply load the data in to
view:

    data("arrowhead_gm12878_5kb")
    head(arrowhead_gm12878_5kb)

The unique boundaries for CHR1 and CHR22 can be extracted by running the
following line:

    bounds <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb, preprocess=FALSE, CHR=c("CHR1","CHR22"),resolution=5000)

Next you will need to download cell line specific ChIP-seq data in the
form of BED files from
[ENCODE](https://www.encodeproject.org/chip-seq-matrix/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&assay_title=TF%20ChIP-seq&status=released).
Once, you have downloaded your preferred list of functional genomic
annotations, store them in a specific file location. These files can
then be converted into a GRangesList object and used for downstream
modelling using the following command:

    path = "pathToBEDfiles"
    tfbsList <- bedToGRangesList(filepath=path, pattern = "*.bed", signal=4)

As an example, we have already provided a GRangesList object with a
variety of transcription factor binding sites specific to the GM12878
cell line. Once you load it in, you can see the list of transcription
factors using the following:

    data("tfbsList")
    names(tfbsList)

For the purposes of this example, lets focus only on CTCF, RAD21, SMC3,
and ZNF143 transcription factors.

    tfbsList <- tfbsList[names(tfbsList) %in% c("Gm12878-Ctcf-Broad","Gm12878-Rad21-Haib","Gm12878-Smc3-Sydh","Gm12878-Znf143-Sydh")]

Now, using the “ground-truth” boundaries and the following TFBS, we can
build the data matrix that will be used for predictive modelling. The
following command creates the training data from CHR1 and reserves the
testing data from CHR22. We specify 5kb sized genomic bins (to match the
resolution used to call the original TADs), a distance-type feature
space, and to apply random under-sampling (RUS) on the training data
only.

    tadData <- createTADdata(bounds.GR=bounds.GR,
                             resolution=5000,
                             genomicElements.GR=tfbsList,
                             featureType="distance",
                             resampling="rus",
                             trainCHR="CHR1",
                             predictCHR="CHR22",
                             seed=123)

We can now implement our machine learning algorithm of choice to predict
TAD-boundary regions. Here, we opt for the random forest algorithm.

    tadModel <- TADrandomForest(trainData=tadData[[1]],
                                testData=tadData[[2]],
                                tuneParams=list(mtry=2,
                                                ntree=500,
                                                nodesize=1),
                                cvFolds=3,
                                cvMetric="Accuracy",
                                verbose=TRUE,
                                seed=123,
                                model=TRUE,
                                importances=TRUE,
                                impMeasure="MDA",
                                performances=TRUE)
                                
    #variable importances
    tadModel[[2]]
    
    #model performance metrics
    tadModel[[3]]

Lastly, we take our TAD-boundary region predictive model and use it to
make predictions on a 10mb section of CHR22:35,000,000-45,000,000.

    #first make sure to restrict the called boundaries to the specific chromosome used to make predictions on
    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   preprocess=FALSE,
                                   CHR="CHR22",
                                   resolution=5000)
                                   
    #run preciseTAD
    pt <- preciseTAD(bounds.GR=bounds.GR,
                     genomicElements.GR=tfbsList,
                     featureType="distance",
                     CHR="CHR22",
                     chromCoords=list(35000000,45000000),
                     tadModel=tadModel[[1]],
                     threshold=1.0,
                     flank=NULL,
                     verbose=TRUE,
                     seed=123,
                     parallel=TRUE,
                     cores=4,
                     splits=4,
                     DBSCAN=TRUE,
                     DBSCAN_params=list(5000,3),
                     method.Clust=NULL,
                     PTBR=FALSE,
                     CLARA=TRUE,
                     method.Dist="euclidean",
                     samples=100,
                     juicer=FALSE)
                     
    #view preciseTAD predicted boundary coordinates between CHR22:35mb-45mb
    pt[[2]]
    
    #view original called TAD boundary coordinates from ARROWHEAD between CHR22:35mb-45mb
    pt[[3]]
