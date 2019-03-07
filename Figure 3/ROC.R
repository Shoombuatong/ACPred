setwd('D:\\Peptide prediction\\Anticancer peptides\\Dataset')
#######Load package
library(AUC)
library(ROCR)
library(prospectr)


par(mfrow = c(2,2 ),oma=c(0, 0, 0, 2))

data = read.csv("RF_5CV.csv", header = TRUE)

pred <- prediction(data[,1],data[,6])
p11 <- performance(pred, 'tpr', 'fpr')
pred <- prediction(data[,2],data[,6])
p12 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,3],data[,6])
p13 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,4],data[,6])
p14 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,5],data[,6])
p15 <- performance(pred,'tpr', 'fpr')

plot(p11, avg="vertical", lwd=2, col="orange", main = "RF (5-fold CV)",
spread.estimate="stderror",plotCI.lwd=1.7, ylab= "True positive rate", font.lab=1,cex.lab=1)
lines(x = c(0,100), y = c(0,100))
plot(p12, avg="vertical", lwd=2, col="forestgreen",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p13, avg="vertical", lwd=2, col="red",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p14, avg="vertical", lwd=2, col="blue",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p15, avg="vertical", lwd=2, col="darkorchid",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)

legend(0.5, 0.65, bty='n', xpd=NA,
       c("AmPseAAC","AAC+PseAA","AAC+AmPseAAC","PseAA+AmPseAAC","AAC+PseAAC+AmPseAAC"), lty = 1, lwd=1.7, 
       col=c("orange","forestgreen","red","blue","darkorchid"),cex=0.8)

data = read.csv("SVM_5CV2.csv", header = TRUE)

pred <- prediction(data[,1],data[,2])
p11 <- performance(pred, 'tpr', 'fpr')
pred <- prediction(data[,3],data[,4])
p12 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,5],data[,6])
p13 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,7],data[,8])
p14 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,9],data[,10])
p15 <- performance(pred,'tpr', 'fpr')

plot(p11, avg="vertical", lwd=2, col="orange",  main = "SVM (5-fold CV)",
spread.estimate="stderror",plotCI.lwd=1.7, ylab= "True positive rate", font.lab=1,cex.lab=1)
lines(x = c(0,100), y = c(0,100))
plot(p12, avg="vertical", lwd=2, col="forestgreen",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p13, avg="vertical", lwd=2, col="red",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p14, avg="vertical", lwd=2, col="blue",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p15, avg="vertical", lwd=2, col="darkorchid",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)

legend(0.5, 0.65, bty='n', xpd=NA,
       c("AmPseAAC","AAC+PseAA","AAC+AmPseAAC","PseAA+AmPseAAC","AAC+PseAAC+AmPseAAC"), lty = 1, lwd=1.7, 
       col=c("orange","forestgreen","red","blue","darkorchid"),cex=0.8)

data = read.csv("RF_LOO.csv", header = TRUE)

pred <- prediction(data[,1],data[,6])
p11 <- performance(pred, 'fpr', 'tpr')
pred <- prediction(data[,2],data[,6])
p12 <- performance(pred,'fpr', 'tpr')
pred <- prediction(data[,3],data[,6])
p13 <- performance(pred,'fpr', 'tpr')
pred <- prediction(data[,4],data[,6])
p14 <- performance(pred,'fpr', 'tpr')
pred <- prediction(data[,5],data[,6])
p15 <- performance(pred,'fpr', 'tpr')

plot(p11, avg="vertical", lwd=2, col="orange", main = "RF (LOOCV)",
spread.estimate="stderror",plotCI.lwd=1.7, ylab= "True positive rate", font.lab=1,cex.lab=1)
lines(x = c(0,100), y = c(0,100))
plot(p12, avg="vertical", lwd=2, col="forestgreen",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p13, avg="vertical", lwd=2, col="red",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p14, avg="vertical", lwd=2, col="blue",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p15, avg="vertical", lwd=2, col="darkorchid",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)

legend(0.5, 0.65, bty='n', xpd=NA,
       c("AmPseAAC","AAC+PseAA","AAC+AmPseAAC","PseAA+AmPseAAC","AAC+PseAAC+AmPseAAC"), lty = 1, lwd=1.7, 
       col=c("orange","forestgreen","red","blue","darkorchid"),cex=0.8)

data = read.csv("SVM_LOO2.csv", header = TRUE)

pred <- prediction(data[,1],data[,2])
p11 <- performance(pred, 'tpr', 'fpr')
pred <- prediction(data[,3],data[,4])
p12 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,5],data[,6])
p13 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,7],data[,8])
p14 <- performance(pred,'tpr', 'fpr')
pred <- prediction(data[,9],data[,10])
p15 <- performance(pred,'tpr', 'fpr')

plot(p11, avg="vertical", lwd=2, col="orange",  main = "SVM (LOOCV)",
spread.estimate="stderror",plotCI.lwd=1.7, ylab= "True positive rate", font.lab=1,cex.lab=1)
lines(x = c(0,100), y = c(0,100))
plot(p12, avg="vertical", lwd=2, col="forestgreen",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p13, avg="vertical", lwd=2, col="red",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p14, avg="vertical", lwd=2, col="blue",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)
plot(p15, avg="vertical", lwd=2, col="darkorchid",
spread.estimate="stderror",plotCI.lwd=1.7,add=TRUE)

legend(0.5, 0.65, bty='n', xpd=NA,
       c("AmPseAAC","AAC+PseAA","AAC+AmPseAAC","PseAA+AmPseAAC","AAC+PseAAC+AmPseAAC"), lty = 1, lwd=1.7, 
       col=c("orange","forestgreen","red","blue","darkorchid"),cex=0.8)
