#' A wrapper function passed to \code{caret::train} to apply a random forest
#' classification algorithm built and tested on user-defined binned domain
#' data from \code{\link{createTADdata}}.
#'
#' @param trainData Data frame, the binned data matrix to built random forest
#' classifiers (can be obtained using \code{\link{createTADdata}})
#' @param testData Data frame, the binned data matrix to test random forest
#' classifiers (can be obtained using \code{\link{createTADdata}}). The first
#' column must be a factor with postitive class "Yes". Default in NULL in which
#' case no performances are evaluated.
#' @param tuneParams List, providing \code{mtry}, \code{ntree}, and
#' \code{nodesize} parameters to feed into \code{\link{randomForest}}. Default
#' is list(mtry = ceiling(sqrt(ncol(trainData) - 1)), ntree = 500, nodesize = 1).
#' If multiple values are provided, then a grid search is performed to tune the
#' model.
#' @param cvFolds Numeric, number of k-fold cross-validation to perform in
#' order to tune the hyperparameters
#' @param cvMetric Character, performance metric to use to choose optimal
#' tuning parameters (one of either "Kappa", "Accuracy", "MCC", "ROC", "Sens",
#' "Spec", "Pos Pred Value", "Neg Pred Value"). Default is "Accuracy"
#' @param verbose Logical, controls whether or not details regarding modelling
#' should be printed out (default in TRUE)
#' @param seed Numeric, controls randomization incurred during data splitting
#' from cross-validation (default is 123)
#' @param model Logical, whether to keep the model object. Default is TRUE
#' @param importances Logical, whether to extract variable importances. Default
#' is TRUE
#' @param impMeasure Character, indicates the variable importance measure to
#' use (one of either "MDA" (mean decrease in accuracy) or "MDG" (mean decrease
#' in gini)). Ignored if importances = FALSE
#' @param performances Logical, indicates whether various performance metrics
#' should be extracted when validating the model on the test data. Ignored if
#' testData = NULL
#'
#' @return A list containg: 1) a train object from \code{caret} with model
#' information, 2) a data.frame of variable importances for each feature
#' included in the model, and 3) a data.frame of various performance metrics
#' @export
#'
#' @importFrom ModelMetrics mcc
#' @importFrom PRROC pr.curve
#' @importFrom pROC roc
#' @importFrom pROC auc
#' @importFrom stats predict
#' @import randomForest caret e1071
#'
#' @examples
#' \dontrun{
#' # Read in ARROWHEAD-called TADs at 5kb
#' data(arrowhead_gm12878_5kb)
#'
#' # Extract unique boundaries
#' bounds.GR <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb,
#'                                preprocess = FALSE,
#'                                CHR = c("CHR21", "CHR22"),
#'                                resolution = 5000)
#'
#' # Read in GRangesList of 26 TFBS
#' data(tfbsList)
#'
#' # Create the binned data matrix for CHR1 (training) and CHR22 (testing)
#' #using 5 kb binning, distance-type predictors from 26 different TFBS from
#' the GM12878 cell line, and random under-sampling
#' tadData <- createTADdata(bounds.GR = bounds.GR,
#'                          resolution = 5000,
#'                          genomicElements.GR = tfbsList,
#'                          featureType = "distance",
#'                          resampling = "rus",
#'                          trainCHR = "CHR21",
#'                          predictCHR = "CHR22",
#'                          seed = 123)
#'
#' # Perform random forest using TADrandomForest by tuning mtry over 10 values
#' using 3-fold CV
#' tadModel <- TADrandomForest(trainData = tadData[[1]],
#'                             testData = tadData[[2]],
#'                             tuneParams = list(mtry = c(2,5,8,10,13,16,18,21,24,26),
#'                                             ntree = 500,
#'                                             nodesize = 1),
#'                             cvFolds = 3,
#'                             cvMetric = "Accuracy",
#'                             verbose = TRUE,
#'                             seed = 123,
#'                             model = TRUE,
#'                             importances = TRUE,
#'                             impMeasure = "MDA",
#'                             performances = TRUE)
#' }
TADrandomForest <- function(trainData,
                            testData=NULL,
                            tuneParams=list(mtry=ceiling(sqrt(ncol(trainData)-1)),
                                            ntree=500,
                                            nodesize=1),
                            cvFolds=5,
                            cvMetric="Accuracy",
                            verbose=FALSE,
                            seed=123,
                            model=TRUE,
                            importances=TRUE,
                            impMeasure="MDA",
                            performances=FALSE){

    ##########################
    #CHECKING FUNCTION INPUTS#
    ##########################

    if(class(trainData)!="data.frame"){print("trainData is not a data.frame object!"); return(0)}
    if(!is.null(testData)){
        if(class(testData)!="data.frame"){print("testData is not a data.frame object!"); return(0)}
    }
    for(i in 1:3){
        if(lapply(tuneParams,class)[[i]]!="numeric"){print("at least 1 component of tuneParams is not a numeric object!"); return(0)}
    }
    if(class(cvFolds)!="numeric"){print("cvFolds is not a numeric object!"); return(0)}
    if(class(cvMetric)!="character"){print("cvMetric is not a character object!"); return(0)}
    if(!(cvMetric %in% c("Kappa","Accuracy","MCC","ROC","Sens","Spec","Neg Pred Value"))){print("cvMetric must be one of either 'Kappa', 'Accuracy', 'MCC','ROC','Sens','Spec','Pos Pred Value','Neg Pred Value'!"); return(0)}
    if(class(verbose)!="logical"){print("verbose is not a logical object"); return(0)}
    if(class(model)!="logical"){print("model is not a logical object!"); return(0)}
    if(class(importances)!="logical"){print("importances is not a logical object!"); return(0)}
    if(class(performances)!="logical"){print("performances is not a logical object!"); return(0)}
    if(class(seed)!="numeric"){print("seed is not a numeric object!"); return(0)}
    if(class(impMeasure)!="character"){print("impMeasure is not a character object!"); return(0)}
    if(!(impMeasure %in% c("MDA","MDG"))){print("impMeasure must be one of either 'MDA' or 'MDG'!"); return(0)}


    #################################
    #PREDICTING TAD BOUNDARY REGIONS#
    #################################

    ########################################################################
    #Establishing summary function to evaluate cross validation performance#
    ########################################################################

    predictiveValues <- function (data, lev = NULL, model = NULL,...){
        PPVobj <- posPredValue(data[, "pred"], data[, "obs"])
        NPVobj <- negPredValue(data[, "pred"], data[, "obs"])
        out <- c(PPVobj, NPVobj)
        #out <- c(NPVobj)
        names(out) <- c("Pos Pred Value", "Neg Pred Value")
        #names(out) <- c("Neg Pred Value")
        out}

    allSummary <- function(data, lev = NULL, model = NULL){
        lvls <- levels(data$obs)

        #mcc
        mcc <- mcc(ifelse(data$obs == lev[2], 0, 1), data[, lvls[1]], cutoff = .5)

        #roc
        b1 <- twoClassSummary(data, lev, model)

        #auprc & f1
        #c1 <- prSummary(data, lev, model)

        #ppv & npv
        d1 <- predictiveValues(data, lev, model)

        #accuracy & kappa
        e1 <- defaultSummary(data, lev, model)

        out <- c(mcc, b1, d1, e1)
        names(out)[1] <- c("MCC")
        out
    }

    ########################################
    #Establishing random forest using caret#
    ########################################

    customRF <- getModelInfo(model = "rf", regex = FALSE)
    customRF$rf$parameters <- data.frame(parameter = c("mtry", "ntree", "nodesize"),
                                         class = rep("numeric", 3),
                                         label = c("mtry", "ntree", "nodesize"))
    customRF$rf$grid <- function(x, y, len = NULL, search = "grid") {}
    customRF$rf$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
        randomForest(x,
                     y,
                     mtry = param$mtry,
                     ntree=param$ntree,
                     nodesize = param$nodesize,
                     importance=TRUE,
                     ...)
    }

    ##########################
    #Establishing tuning grid#
    ##########################

    tunegrid <- expand.grid(mtry=tuneParams[[1]],
                            ntree=tuneParams[[2]],
                            nodesize=tuneParams[[3]])

    ###############################
    #Establishing control function#
    ###############################

    set.seed(seed)
    seeds <- vector(mode = "list", length = (cvFolds+1))
    for(i in 1:cvFolds){
        #set.seed(seed)
        seeds[[i]]<- sample.int(n=100000, nrow(tunegrid))
    }
    #for the last model
    set.seed(seed)
    seeds[[cvFolds+1]]<-sample.int(100000, 1)

    control <- trainControl(seeds = seeds,
                            method = "cv",
                            number = cvFolds,
                            verboseIter = ifelse(verbose==TRUE, TRUE, FALSE),
                            ## Estimate class probabilities
                            classProbs = TRUE,
                            ## Evaluate performance using
                            ## the following function
                            summaryFunction = allSummary,
                            allowParallel = FALSE)


    set.seed(seed)
    tadModel <- train(y~.,
                      data=trainData,
                      method=customRF$rf,
                      metric=cvMetric,
                      tuneGrid=tunegrid,
                      trControl=control)

    #####################################
    #EXTRACTING IMPORTANCES/COEFFICIENTS#
    #####################################

    if(impMeasure=="MDA"){
        rfimpvars <- data.frame(Feature=rownames(randomForest::importance(tadModel$finalModel, scale=TRUE)),
                                Importance=as.vector(randomForest::importance(tadModel$finalModel, scale=TRUE)[,3]))
        rfimpvars <- rfimpvars[order(rfimpvars$Importance, decreasing=TRUE),]
    }else if(impMeasure=="MDG"){
        rfimpvars <- data.frame(Feature=rownames(randomForest::importance(tadModel$finalModel, scale=TRUE)),
                                Importance=as.vector(randomForest::importance(tadModel$finalModel, scale=TRUE)[,4]))
        rfimpvars <- rfimpvars[order(rfimpvars$Importance, decreasing=TRUE),]
    }else{
        rfimpvars <- NA
    }

    #################################
    #EVALUATING PERFORMANCES        #
    #################################

    if(!(is.null(testData)) & performances==TRUE){
        rfperf <- data.frame(Metric = c("TN",
                                        "FN",
                                        "FP",
                                        "TP",
                                        "Total",
                                        "Sensitivity",
                                        "Specificity",
                                        "Kappa",
                                        "Accuracy",
                                        "BalancedAccuracy",
                                        "Precision",
                                        "FPR",
                                        "FNR",
                                        "NPV",
                                        "MCC",
                                        "F1",
                                        "AUC",
                                        "Youden",
                                        "AUPRC"),
                             Performance=NA)

        pred.tadModel <- as.vector(predict(tadModel,
                                           newdata=testData[,-1],
                                           type="prob")[,"Yes"])

        fg <- pred.tadModel[testData$y == "Yes"]
        bg <- pred.tadModel[testData$y == "No"]
        pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

        pred.tadModel2 <- predict(tadModel,
                                  newdata=testData[,-1],
                                  type="raw")


        confMat <- confusionMatrix(data=pred.tadModel2, testData[,1], positive="Yes")
        TN = as.numeric(confMat$table[1,1])
        FN = as.numeric(confMat$table[1,2])
        FP = as.numeric(confMat$table[2,1])
        TP = as.numeric(confMat$table[2,2])
        rfperf[1,2] <- TN
        rfperf[2,2] <- FN
        rfperf[3,2] <- FP
        rfperf[4,2] <- TP
        rfperf[5,2] <- sum(confMat$table)
        rfperf[6,2] <- as.vector(confMat$byClass["Sensitivity"])
        rfperf[7,2] <- as.vector(confMat$byClass["Specificity"])
        rfperf[8,2] <- as.vector(confMat$overall["Kappa"])
        rfperf[9,2] <- as.vector(confMat$overall["Accuracy"])
        rfperf[10,2] <- ((TP/(TP+FN)) + (TN/(FP+TN)))/2
        rfperf[11,2] <- TP/(TP+FP)
        rfperf[12,2] <- FP/(FP+TN)
        rfperf[13,2] <- FN/(FN+TN)
        rfperf[15,2] <- TN/(TN+FN)
        rfperf[16,2] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
        rfperf[17,2] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
        rfperf[18,2] <- pROC::auc(pROC::roc(testData$y, pred.tadModel, quiet = TRUE))
        rfperf[19,2] <- (TP/(TP + FN)) + (TN/(TN + FP)) - 1
        rfperf[20,2] <- pr$auc.integral
    }

    results_list <- list()

    if(model==TRUE){results_list[[1]] = tadModel}else{results_list[[1]] = NA}
    if(importances==TRUE){results_list[[2]] = rfimpvars}else{results_list[[2]] = NA}
    if(!(is.null(testData)) & performances==TRUE){results_list[[3]] = rfperf}else{results_list[[3]] = NA}

    return(results_list)
}
