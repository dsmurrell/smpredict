#!/usr/bin/Rscript
sink("log", append=FALSE, split=TRUE)
library(smPredict)

df <- read.csv("./external_tool_input.smi", header=FALSE)

input.file <- tempfile("input", fileext=".smi")

write.table(df, input.file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

prediction <- PredictLogP(input.file, threads = 1)

write.table(prediction,"./external_tool_output.csv",row.names=FALSE,col.names=TRUE,sep=",", quote=FALSE)
