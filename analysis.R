library(igraph)
library(RedeR)

createNetworkCor <- function(expr.data, cutoff) {
    # Takes an expr.data matrix and returns a igraph where
    # nodes are the genes in expr.data. Nodes are connected
    # if their correlation coeff is larger than cor.cutoff
    expr.mat <- t(as.matrix(expr.data))
    cor.mat <- cor(expr.mat)

    adj.mat  <- (abs(cor.mat > cutoff)) * 1

    return(adj.mat)
}

getOlaps  <- function(badj.mat, k=10) {
  tf.olaps <- badj.mat %*% t(badj.mat)
  diag(tf.olaps)  <- 0
  olaps <- unlist(mapply(rep, names(table(tf.olaps)), table(tf.olaps)))
  l.olaps <- unique(sort(as.numeric(olaps), decreasing=T))

  return(tf.olaps)
}

getBiadjacencyMat <- function(adj.mat, tfs) {
    tfs.ind <- tfs[tfs=="TRUE"]
    genes.ind <- tfs[tfs=="FALSE"]

    return(adj.mat[tfs.ind, tfs])
}

getEdges <- function(adj.mat) {
    edges <- which(adj.mat==1, arr.ind = T)
    edges <- unlist(mapply(c, edges[, 1], edges[,2], SIMPLIFY = FALSE))

    return(edges)
}

load("data/annotation.RData")
expr.data <- read.table(file="data/disc_set/discovery_ExpressionMatrix_red.txt",
                        header=T, comment.char="", row.names=1)
tfs <- read.table("data/disc_set/tfs.txt")
adj.mat  <- createAdjMat(expr.data, 0.5)


rdp <- RedPort()
calld(rdp)
addGraph(rdp, g)


