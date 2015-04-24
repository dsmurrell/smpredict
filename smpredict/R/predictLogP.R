ReplaceInfinitesWithNA <- function(d) {
  do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x), NA)))
}

PredictPropertytoCSV <- function(property = c('LogP', 'LogS'), csv.file, structures.file, error.variance=FALSE, threads = -1) {
  property <- match.arg(property)
  write.csv(PredictProperty(property, structures.file, error.variance, threads), file = csv.file)
}

PredictProperty <- function(property = c('LogP', 'LogS'), structures.file, error.variance=FALSE, threads = -1) {
  property <- match.arg(property)
  standardised.file <- tempfile("standardised", fileext=".sdf")
  name.file <- tempfile("name", fileext=".txt") # used to recover the original ordering of molecules (multicore descriptor generation messes this up)
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  StandardiseMolecules(structures.file, standardised.file, name.file = name.file, limit = -1)
  GenerateDescriptors.internal(standardised.file, descriptors.file, name.file, threads)
  PredictPropertyFromDescriptorsFile(property, descriptors.file, error.variance)
}

PredictPropertyFromDescriptorsFile <- function(property = c('LogP', 'LogS'), descriptors.file, error.variance=FALSE) {
  property <- match.arg(property)
  descriptors <- read.csv(file = descriptors.file)
  PredictPropertyFromDescriptors(property, descriptors, error.variance)
}

PredictPropertyFromDescriptors <- function(property = c('LogP', 'LogS'), descriptors, error.variance=FALSE, error.cores=1) {
  property <- match.arg(property)
  do.call(what=paste("Predict", property, "FromDescriptors", sep=''), args=list(descriptors,error.variance,error.cores))  
}

PredictLogPFromDescriptors <- function(descriptors, error.variance=FALSE, error.cores=1) {
  print("Predicting LogP")
  load(file = system.file("extdata", "logp_data.rda", package="smpredict"))
  load(file = system.file("extdata", "logp_svm.rda", package="smpredict"))
  load(file = system.file("extdata", "logp_gbm.rda", package="smpredict"))
  names <- descriptors[, 1]
  used.descriptors <- names(x.train)
  descriptors <- descriptors[, used.descriptors]
  suppressWarnings(descriptors <- apply(descriptors, 2, as.numeric))
  nrows <- nrow(descriptors)
  if(is.null(nrows)) { nrows <- 1 }
  to.impute <- rbind(descriptors, x.train)
  to.impute <- ReplaceInfinitesWithNA(to.impute)
  set.seed(777)
  imputed <- as.data.frame(impute.knn(as.matrix(to.impute), k = 10)$data)
  to.predict <- imputed[1:nrows, ]
  to.predict <- rbind(to.predict[1,], to.predict) # HACK
  x <- predict(transformation, to.predict)
  svm_pred <- predict(svm$finalModel, newdata = x)
  gbm_pred <- predict(gbm$finalModel, newdata = x, n.trees = 500)
  greedy_pred <-  gbm_pred*0.288 + svm_pred*0.712
  r <- data.frame(ID = names, smLogP = greedy_pred[2:length(greedy_pred)]) # HACK
  r$smLogP <- format(round(r$smLogP, 2), nsmall = 2)
  if(error.variance) {
    library(eve)
    estimator <- readRDS(system.file("extdata", "logp_eve.rds", package="smpredict"))
    sigmas <- PredictSigmasMC(x=x, estimator=estimator, cores=error.cores)
    r$error_variance <- format(round(sigmas, 3), nsmall = 3)
  }
  r
}

PredictLogSFromDescriptors <- function(descriptors, error.variance=FALSE, error.cores=1) {
  print("Predicting LogS")
  dataset <- readRDS(system.file("extdata", "logs_dataset.rds", package="smpredict"))
  x.train <- dataset$x.train
  transformation <- dataset$transformation
  svmRadial <- readRDS(system.file("extdata", "logs_svmRadial.rds", package="smpredict")) 
  gbm <- readRDS(system.file("extdata", "logs_gbm.rds", package="smpredict"))
  cubist <- readRDS(system.file("extdata", "logs_cubist.rds", package="smpredict"))
  names <- descriptors[, 1]
  used.descriptors <- names(x.train)
  descriptors <- descriptors[, used.descriptors]
  suppressWarnings(descriptors <- apply(descriptors, 2, as.numeric))
  nrows <- nrow(descriptors)
  if(is.null(nrows)) { nrows <- 1 }
  to.impute <- rbind(descriptors, x.train)
  to.impute <- ReplaceInfinitesWithNA(to.impute)
  set.seed(777)
  imputed <- as.data.frame(impute.knn(as.matrix(to.impute), k = 10)$data)
  to.predict <- imputed[1:nrows, ]
  to.predict <- rbind(to.predict[1,], to.predict) # HACK
  x <- predict(transformation, to.predict)
  svmRadial_pred <- predict(svmRadial$finalModel, newdata = x)
  gbm_pred <- predict(gbm$finalModel, newdata = x, n.trees = 500)
  library(Cubist)
  cubist_pred <- predict(cubist$finalModel, newdata = x)
  greedy_pred <-  gbm_pred*0.285 + svmRadial_pred*0.318 + cubist_pred*0.354
  r <- data.frame(ID = names, smLogS = greedy_pred[2:length(greedy_pred)])
  r$smLogS <- format(round(r$smLogS, 2), nsmall = 2)
  #if(error.variance) {
  #  library(eve)
  #  estimator <- readRDS(system.file("extdata", "logp_eve.rds", package="smpredict"))
  #  sigmas <- PredictSigmasMC(x=x, estimator=estimator, cores=error.cores)
  #  r$error_variance <- format(round(sigmas, 3), nsmall = 3)
  #}
  r
}
