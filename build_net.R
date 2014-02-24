library(psych)
library(igraph)
library(corpcor)

pcorTest <- function(expr.mat) {
    r <- cor(expr.mat)
    n <- nrow(expr.mat)

    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- 2 * (1 - pt(abs(t), (n - 2)))
    p <- p.adjust(p, "holm")

    return(list(r=r, p=p))
}

corPermTest <- function(expr.mat) {
    #test using permutation tests
    expr.mat1 <- expr.mat[sample(nrow(expr.mat)), ]
    return(expr.mat1)
}

buildNet <- function(cor.res, tfs, p.cutoff) {
    sign <- cor.res$p < p.cutoff

    signf <- which(sign == FALSE)
    signt <- which(sign == TRUE)

    cor.res$r[signt] = 1
    cor.res$r[signf] = 0
    adj.mat <- cor.res$r

    badj.mat <- adj.mat[which(tfs == TRUE), which(tfs==FALSE)]
    g <- graph.incidence(badj.mat)

    return(list(g, badj.mat))
}

buildNetCor <- function(expr.data, tfs, p.cutoff) {
    expr.mat <- t(as.matrix(expr.data))
    cor.res <- corr.test(expr.mat)

    return(buildNet(cor.res, tfs, p.cutoff))
}

buildNetPCor <- function(expr.data, tfs, p.cutoff) {
    expr.mat <- t(as.matrix(expr.data))
    pcor.res <- pcorTest(expr.mat)

    return(buildNet(pcor.res, tfs, p.cutoff))
}

runNetworkBuilds <- function() {
    load("data/annotation.RData")
    expr.data <- read.table(file="data/disc_set/discovery_ExpressionMatrix_red.txt",
                            header=T, comment.char="", row.names=1)
    tfs <- read.table("data/disc_set/tfs.txt")

    cor.g <- buildNetCor(expr.data, tfs, 0.01)
    pcor.g <- buildNetPCor(expr.data, tfs, 0.01)
}

olapAnalysis <- function(badj.mat) {
    olaps <- badj.mat %*% t(badj.mat)


}


