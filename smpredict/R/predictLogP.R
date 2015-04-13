findMinDistances <- function(set, referenceSet, numberOfNeighbours) {
  n <- numberOfNeighbours
  x <- as.matrix(referenceSet)
  minDistances <- c()
  for(i in 1:nrow(set)) {
    y <- as.vector(set[i,], mode="numeric")
    y[is.na(y)] <- 0
    distances <- colSums((t(x)-y)^2) 
    minDistances <- c(minDistances, sum(sort(sqrt(distances),partial=n)[1:n])/n)
  }
  minDistances[is.na(minDistances)] <- 0
  minDistances
}

ReplaceInfinitesWithNA <- function(d) {
  do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x), NA)))
}

PredictLogPFromDescriptors <- function(descriptors, error.variance=FALSE, error.cores=1) {
  print("Predicting LogP")
  load(file = system.file("extdata", "data.rda", package="smpredict"))
  load(file = system.file("extdata", "svm.rda", package="smpredict"))
  load(file = system.file("extdata", "gbm.rda", package="smpredict"))
  names <- descriptors[, 1]
  used.descriptors <- names(x.train)
  descriptors <- descriptors[, used.descriptors]
  suppressWarnings(descriptors <- apply(descriptors, 2, as.numeric))
  nrows <- nrow(descriptors)
  to.impute <- rbind(descriptors, x.train)
  to.impute <- ReplaceInfinitesWithNA(to.impute)
  set.seed(777)
  imputed <- as.data.frame(impute.knn(as.matrix(to.impute), k = 10)$data)
  to.predict <- imputed[1:nrows, ]
  x <- predict(transformation, to.predict)
  svm_pred <- predict(svm$finalModel, newdata = x)
  gbm_pred <- predict(gbm$finalModel, newdata = x, n.trees = 500)
  greedy_pred <- svm_pred*0.712 + gbm_pred*0.288
  r <- data.frame(ID = names, smLogP = greedy_pred)
  if(error.variance) {
    library(eve)
    estimator <- readRDS(system.file("extdata", "logPeve.rds", package="smpredict"))
    sigmas <- PredictSigmasMC(x=x, estimator=estimator, cores=error.cores)
    r$sigma <- sigmas
  }
  #if(error.variance) {
  #  n <- 5
  #  p <- 1.5
  #  load(file = system.file("extdata", "spline.rda", package="smpredict"))
  #  load(file = system.file("extdata", "as.rda", package="smpredict"))
  #  outset <- to.predict
  #  for(i in 1:length(as)) {
  #    outset[,i] <- outset[,i] * (as[i]^p)
  #  }
  #  inset <- x.train
  #  for(i in 1:length(as)) {
  #    inset[,i] <- inset[,i] * (as[i]^p)
  #  }
  #  distances <- findMinDistances(outset, inset, n)
  #  sf <- 15/mean(distances)
  #  distances <- distances*sf
  #  variances <- predict(spline, distances)$y
  #  r$variance <- variances
  #}
  r
}

PredictLogPFromDescriptorsFile <- function(descriptors.file, error.variance=FALSE) {
  descriptors <- read.csv(file = descriptors.file)
  PredictLogPFromDescriptors(descriptors, error.variance)
}

PredictLogP <- function(structures.file, error.variance=FALSE, threads = -1) {
  standardised.file <- tempfile("standardised", fileext=".sdf")
  name.file <- tempfile("name", fileext=".txt") # used to recover the original ordering of molecules (multicore descriptor generation messes this up)
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  StandardiseMolecules(structures.file, standardised.file, name.file = name.file, limit = -1)
  GenerateDescriptors.internal(standardised.file, descriptors.file, name.file, threads)
  PredictLogPFromDescriptorsFile(descriptors.file, error.variance)
}

PredictLogPtoCSV <- function(csv.file, structures.file, error.variance=FALSE, threads = -1) {
  write.csv(PredictLogP(structures.file, error.variance, threads), file = csv.file)
}
