setwd('D:\\Peptide prediction\\Anticancer peptides\\Dataset')

library(caret)
library(randomForest)
library(rpart)
library(RRF)
library(inTrees)
library(party)
library(partykit)
library(xtable)

#######Extract AAC DPC TPC
x <- read.fasta('ZOH.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Label.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)
AAC <- t(sapply(A, extractAAC))

########## Generate rule from RF
D = data.frame(round(AAC,3),Class = label)
D = na.omit(D)

X <- D[,-ncol(D)]
target <- D[,ncol(D)]
rf <- RRF(X,as.factor(target),ntree=500, mtry = 7) # build an ordinary RF
treeList <- RF2List(rf)
ruleExec <- extractRules(treeList,X)
ruleExec <- unique(ruleExec)
ruleMetric <- getRuleMetric(ruleExec,X,target) # measure rules
ruleMetric <- pruneRule(ruleMetric,X,target) # prune each rule
#ruleMetric <- selectRuleRRF(ruleMetric,X,target) # rule selection
learner <- buildLearner(ruleMetric,X,target)
pred <- applyLearner(learner,X)
read <- presentRules(learner,colnames(X)) # more readable format
# format the rule and metrics as a table in latex code
data.frame(xtable(read))
