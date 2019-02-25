#######set directory
setwd('D:\\Peptide prediction\\Anticancer peptides\\Dataset')
#######Load package
library(Interpol)
library(protr)
library(seqinr)

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

### Please note that PAAC and APAAC are PseAAC and Am-PseAAC, respectively.
for(i in 1:m){ 
PAAC[i, ] = extractPAAC(A[[i]][1],lambda = pse, w = weight, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

col = 20+ 2*4
APAAC  <- matrix(nrow = length(A), ncol = col)
for (i in 1:length(A)){
APAAC[i,] = extractAPAAC(A[[i]][1],lambda = 4, w = 0.01, customprops = NULL)
}
