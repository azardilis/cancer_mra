library(igraph)
library(RedeR)

createNetworkCor <- function(expr.data, cor.cutoff) {
    # Takes an expr.data matrix and returns a igraph where
    # nodes are the genes in expr.data. Nodes are connected
    # if their correlation coeff is larger than cor.cutoff
    cor.mat <- cor(expr.data, expr.data)
    adj.mat <- (abs(cor.mat) > cor.cutoff) * 1
    cor.graph <- graph.adjacency(adj.mat)

    return(cor.graph)
}



#first load the annotation and expression data
#reduce expression data to matrix to only include annotated probes


#fake dataset
expr.data <- matrix(rnorm(25), nrow=5)
cor.graph <- createNetworkCor(expr.data, 0.75)


#visualise network using RederR
rdp <- RedPort()
calld(rdp)
addGraph(rdp,cor.graph)

    
