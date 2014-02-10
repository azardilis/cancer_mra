library(igraph)
library(RedeR)

createNetwork <- function(dat, fun, cutoff) {
    # Takes an expr.data matrix and returns a igraph where
    # nodes are the genes in expr.data. Nodes are connected
    # if their correlation coeff is larger than cor.cutoff
    cor.mat <- fun(dat$tfs, dat$genes)
    adj.mat <- (abs(cor.mat) > cutoff) * 1
    edges <- which(adj.mat==1, arr.ind = T)
    edges <- unlist(mapply(c, edges[, 1], edges[,2], SIMPLIFY = FALSE))
    g <- graph(edges=edges,directed = FALSE)

    return(g)
}


splitSets <- function(expr.data, annotation) {
    expr.data <- cbind(expr.data, annotation$is.TF)
    tfs <- expr.data[expr.data$annotation == TRUE, ]
    genes <- expr.data[expr.data$annotation== FALSE, ]
    colnames(genes)[998] <- "annotation"
    genes$annotation <- NULL
    colnames(tfs)[998] <- "annotation"
    tfs$annotation <- NULL
    tfs.mat <- t(as.matrix(tfs))
    genes.mat <- t(as.matrix(genes))

    return(list(tfs=tfs.mat, genes=genes.mat))
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
tfs <- read.table("/home/argyris/compbio/nb/cancer_mra/data/disc_set/tfs.txt")
expr.mat <- t(as.matrix(expr.data))
cor.mat <- cor(expr.mat)



rdp <- RedPort()
calld(rdp)
addGraph(rdp, g)

    
