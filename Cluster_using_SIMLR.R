library("Matrix")
library("parallel")
library("igraph")
library("grDevices")
# load SIMLR
source("/Users/koukansei/Desktop/SIMLR-SIMLR/R/SIMLR.R")
source("/Users/koukansei/Desktop/SIMLR-SIMLR/R/compute.multiple.kernel.R")
source("/Users/koukansei/Desktop/SIMLR-SIMLR/R/network.diffusion.R")
source("/Users/koukansei/Desktop/SIMLR-SIMLR/R/utils.simlr.R")
source("/Users/koukansei/Desktop/SIMLR-SIMLR/R/tsne.R")

dyn.load("/Users/koukansei/Desktop/SIMLR-SIMLR/R/projsplx_R.so")
setwd("/Users/koukansei/Desktop/SIMLR-SIMLR/R")
data <- read.csv("/Users/koukansei/Desktop/对比的方法/data/10x-5cl.csv", header = FALSE)
label <- read.csv("/Users/koukansei/Desktop/对比的方法/label/10x-5cl.csv", header = FALSE)
label <-as.numeric(unlist(label))  #想直接算NMI的必要步骤
data <- t(data)
result = SIMLR(X = data, c = max(label) - min(label) + 1)
write.table(result$S, file = "/Users/koukansei/Desktop/SIMLR-SIMLR/SIMLR-10x-5cl.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
y_p = result$y$cluster

write.table(y_p, file = "/Users/koukansei/Desktop/SIMLR-SIMLR/10x-5cl.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

