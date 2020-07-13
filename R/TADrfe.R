#' A wrapper function passed to \code{caret::rfe} to apply recursive feature
#' elimination (RFE) on binned domain data as a feature reduction technique for
#' random forests. Backward elimination is performed from p down to 2, by
#' powers of 2, where p is the number of features in the data.
#'
#' @param trainData Data frame, the binned data matrix to built a random forest
#' classifiers (can be obtained using \code{\link{createTADdata}})
#' @param tuneParams List, providing \code{ntree} and \code{nodesize}
#' parameters to feed into \code{\link{randomForest}}
#' @param cvFolds Numeric, number of k-fold cross-validation to perform
#' @param cvMetric Character, performance metric to use to choose optimal
#' tuning parameters (one of either "Kappa", "Accuracy", "MCC","ROC","Sens",
#' "Spec", "Pos Pred Value", "Neg Pred Value"). Default is "Accuracy"
#' @param verbose Logical, controls whether or not details regarding modeling
#' should be printed out (default is TRUE)
#' @param seed Numeric, controls randomization incurred during data splitting
#' from cross-validation (default is 123)
#'
#' @return A list containing: 1) the performances extracted at each of the k
#' folds and, 2) Variable importances among the top features at each step of
#' RFE. For 1) `Variables` - the best subset of features to consider at each
#' iteration, `MCC` (Matthews Correlation Coefficient), `ROC` (Area under the
#' receiver operating characteristic curve), `Sens` (Sensitivity), `Spec`
#' (Specificity), `Pos Pred Value` (Positive predictive value), `Neg Pred Value`
#' (Negative predictive value), `Accuracy`, and the corresponding standard
#' deviations across the cross-folds. For 2) `Overall` - the variable
#' importance, `var` - the feature name, `Variables` - the number of features
#' that were considered at each cross-fold, and `Resample` - the cross-fold
#' @export
#'
#' @importFrom ModelMetrics mcc
#' @import randomForest caret e1071
#'
#' @examples
#' \dontrun{
#' # Read in ARROWHEAD-called TADs at 5kb
#' data(arrowhead_gm12878_5kb)
#'
#' #Extract unique boundaries
#' bounds.GR <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb,
#'                                preprocess = FALSE,
#'                                CHR = "CHR22",
#'                                resolution = 5000)
#'
#' # Read in GRangesList of 26 TFBS
#' data(tfbsList)
#'
#' # Create the binned data matrix for CHR22 using:
#' # 5 kb binning,
#' # oc-type predictors from 26 different TFBS from the GM12878 cell line, and
#' # random under-sampling
#' tadData <- createTADdata(bounds.GR = bounds.GR,
#'                          resolution = 5000,
#'                          genomicElements.GR = tfbsList,
#'                          featureType = "oc",
#'                          resampling = "rus",
#'                          trainCHR = "CHR22",
#'                          predictCHR = NULL,
#'                          seed=123)
#'
#' # Perform RFE for fully grown random forests with 100 trees using 5-fold CV
#' # Evaluate performances using accuracy
#' rfe_res <- TADrfe(trainData = tadData[[1]],
#'                   tuneParams = list(ntree = 100, nodesize = 1),
#'                   cvFolds = 5,
#'                   cvMetric = "Accuracy",
#'                   verbose = TRUE,
#'                   seed = 123)
#' }
TADrfe <- function(trainData, tuneParams = list(ntree = 500, nodesize = 1),
                   cvFolds = 5, cvMetric = "Accuracy", verbose = FALSE, seed = 123) {

    #CHECKING FUNCTION INPUTS#

    if (class(trainData) != "data.frame") {
        print("trainData is not a data.frame object!")
        return(0)
    }
    for (i in 1:2) {
        if (lapply(tuneParams, class)[[i]] != "numeric") {
            print("at least 1 component of tuneParams is not a numeric object!")
            return(0)
        }
    }
    if (class(cvFolds) != "numeric") {
        print("cvFolds is not a numeric object!")
        return(0)
    }
    if (class(cvMetric) != "character") {
        print("cvMetric is not a character object!")
        return(0)
    }
    if (!(cvMetric %in% c("Kappa", "Accuracy", "MCC", "ROC", "Sens", "Spec",
                          "Neg Pred Value"))) {
        print("cvMetric must be one of either 'Kappa', 'Accuracy', 'MCC','ROC','Sens','Spec','Pos Pred Value','Neg Pred Value'!")
        return(0)
    }
    if (class(verbose) != "logical") {
        print("verbose is not a logical object")
        return(0)
    }
    if (class(seed) != "numeric") {
        print("seed is not a numeric object!")
        return(0)
    }

    #Establishing summary function#

    predictiveValues <- function(data, lev = NULL, model = NULL, ...) {
        PPVobj <- posPredValue(data[, "pred"], data[, "obs"])
        NPVobj <- negPredValue(data[, "pred"], data[, "obs"])
        out <- c(PPVobj, NPVobj)
        # out <- c(NPVobj)
        names(out) <- c("Pos Pred Value", "Neg Pred Value")
        # names(out) <- c('Neg Pred Value')
        out
    }

    allSummary <- function(data, lev = NULL, model = NULL) {
        lvls <- levels(data$obs)

        # mcc
        mcc <- mcc(ifelse(data$obs == lev[2], 0, 1), data[, lvls[1]], cutoff = 0.5)

        # roc
        b1 <- twoClassSummary(data, lev, model)

        # auprc & f1 c1 <- prSummary(data, lev, model)

        # ppv & npv
        d1 <- predictiveValues(data, lev, model)

        # accuracy & kappa
        e1 <- defaultSummary(data, lev, model)

        out <- c(mcc, b1, d1, e1)
        names(out)[1] <- c("MCC")
        out
    }

    #PERFORMING RFE#

    rfectrl <- rfFuncs
    rfectrl$summary <- allSummary
    rfFuncs$fit <- function(x, y, param, first, last, ...) {
        randomForest::randomForest(x, y, ntree = tuneParams[[1]], nodesize = tuneParams[[2]],
                                   importance = TRUE, ...)
    }
    control <- rfeControl(functions = rfectrl, method = "cv", number = cvFolds,
                          verbose = ifelse(verbose == TRUE, TRUE, FALSE), allowParallel = FALSE)
    control$returnResamp <- "final"

    n = dim(trainData)[2] - 1
    z <- numeric()
    x = 0
    i = 1
    while (x < n) {
        x = 2^(i)
        i = i + 1
        z <- c(z, x)
    }
    z[length(z)] <- n

    set.seed(seed)
    tadModel <- rfe(trainData[, -1], trainData[, 1], metric = cvMetric,
                    sizes = z, rfeControl = control)

    rfeModelResultsList <- list(tadModel$results, tadModel$variables[,
                                                                     -c(1, 2)])
    names(rfeModelResultsList) <- c("CVPerformances", "Importances")

    return(rfeModelResultsList)
}
