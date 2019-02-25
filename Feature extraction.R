#######set directory
setwd('D:\\Peptide prediction\\Anticancer peptides\\Dataset')
#######Load package
library(randomForest)
library(protr)
library(seqinr)
library(C50)
library(RWeka)
library(Interpol)
library(RWeka)
library(caret)
library(randomForest)
library(kernlab)
library(corrplot)
library(C50)
library(nnet)
library(e1071)
library(GA)
library(cvTools) 
library(Metrics)
library(MASS)
library(pls)
library(Interpol)
library(protr)
library(seqinr)
library(Peptides)
library(AUC)
library(ROCR)
library(mltools)
library(DMwR)
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

#######Extract AAC DPC TPC
x <- read.fasta('ZOH.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Label.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)
PCP  <- matrix(nrow = m, ncol = 531)
pse = 3
weight = 0.4
PAAC <- matrix(nrow = m, ncol = 20 + pse)

AAC <- t(sapply(A, extractAAC))
DPC <- t(sapply(A, extractDC))
TPC <- t(sapply(A, extractTC))

for(i in 1:m){ 
x = A[[i]][1]
b = AAdescriptor(x, 531,2)
PCP[i,] = Interpol(b, 531,"linear")
}

for(i in 1:m){ 
PAAC[i, ] = extractPAAC(A[[i]][1],lambda = pse, w = weight, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

col = 20+ 2*4
APAAC  <- matrix(nrow = length(A), ncol = col)
for (i in 1:length(A)){
APAAC[i,] = extractAPAAC(A[[i]][1],lambda = 4, w = 0.01, customprops = NULL)
}
