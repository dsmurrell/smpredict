#### TESTING PACKAGE ON 50 TEST STRUCTURES TAKEN FROM THE PKKB DATASET
detach(package:smpredict, unload=TRUE)
library(smpredict)

system.time(predictions <- PredictLogP(structures.file="structures.sdf", error.variance=TRUE, threads = 1))

predictions
plot(density(predictions$variance))

predictions[1,]

library(ggplot2)
p <- predictions[1,]
smLogP <- p$smLogP
sigma <- sqrt(p$variance)
x=seq(smLogP-6*sigma,smLogP+6*sigma,length=400)
A <- (1 - 0)/(1/(sigma*2*pi))*exp(-0.5*(0/sigma)^2)
y <- A*0.98*(1/(sigma*2*pi))*exp(-0.5*(x-smLogP/sigma)^2)

p <- ggplot() + xlim(min(x), max(x)) + ylim(0,1)

gaussian <- data.frame(x=x, y=y)
p <- p + geom_line(aes(x = x, y = y), size=1.5, colour = "red", data=gaussian)
plot(p)




# 21.7 seconds

logP <- read.csv("logP.csv")
m <- merge(logP, predictions)

plot(m$logP, m$smLogP)
abline(0,1)

.libPaths(c("/home/dsm38/R/x86_64-pc-linux-gnu-library/3.0",
            "/shared/shared/ubuntu-12.04/R/site-library",
            "/usr/local/lib/R/site-library",
            "/usr/lib/R/site-library",                      
            "/usr/lib/R/library"))

?.libPaths

?.onAttach