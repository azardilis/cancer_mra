library(igraph)
library(corpcor)
library(GeneNet)

corTest <- function(expr.mat) {
    r <- cor(expr.mat)
    n <- nrow(expr.mat)

    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- 2 * (1 - pt(abs(t), (n - 2)))
    p <- p.adjust(p, "holm")

    return(list(r=r, p=p))
}

buildNet <- function(cor.res, tfs, p.cutoff) {
    sign <- cor.res$p < p.cutoff

    signf <- which(sign == FALSE)
    cor.res$r[signf] = 0
    adj.mat <- cor.res$r

    genes <- which(!(colnames(adj.mat) %in% tfs))
    badj.mat <- adj.mat[tfs, genes]
    g <- graph.incidence(badj.mat, weighted="weight")

    return(g)
}

buildNetCor <- function(expr.data, tfs, p.cutoff) {
    expr.mat <- t(as.matrix(expr.data))
    cor.res <- corTest(expr.mat)

    return(buildNet(cor.res, tfs, p.cutoff))
}

buildNetPCor <- function(expr.data, tfs, p.cutoff) {
    expr.mat <- t(as.matrix(expr.data))
    inf.cor <- cor(expr.mat)
    inf.pcor <- cor2pcor(inf.cor)

    test.results <- ggm.test.edges(inf.pcor, plot=FALSE)
    net.edges <- test.results[test.results$pval < p.cutoff, ]

    nodes1 <- colnames(inf.cor)[net.edges$node1]
    nodes2 <- colnames(inf.cor)[net.edges$node2]

    net.edges.df <- data.frame(nodes1, nodes2, weight=net.edges$pcor)

    g <- graph.data.frame(net.edges.df)

    adj <- as.matrix(get.adjacency(g, attr="weight"))
    genes <- which(!(colnames(adj) %in% tfs))
    badj.mat <- adj[tfs, genes]
    g.bp <- graph.incidence(badj.mat, weighted="weight")

    return(g.bp)
}

buildNetMI <- function(g.MI, tfs, annotation) {
    adj <- as.matrix(get.adjacency(g.MI, attr="weight"))
    pids <- annotation$probeID[which(annotation$probeID %in% rownames(adj))]
    tfs <- tfs[which(tfs %in% pids)]
    genes <- which(!(pids %in% tfs))
    badj.mat <- adj[tfs, genes]
    g <- graph.incidence(badj.mat, weighted="weight")

    return(g)
}


getBadj  <- function(g, weight) {
    if(weight == TRUE) {
      adj.mat  <- as.matrix(get.adjacency(g, attr="weight"))
    } else {
      adj.mat  <- as.matrix(get.adjacency(g))
    }
    tfs <- which(get.vertex.attribute(g, "type") == FALSE)
    
    badj.mat  <- adj.mat[tfs, -tfs]
    
    return(badj.mat)
}

olapAnalysis <- function(badj.mat, vis=FALSE) {
    olaps <- badj.mat %*% t(badj.mat)
    l.regs <- sort.int(rowSums(olaps), decreasing=T, index.return=T)$ix[1:10]
    l.olaps <- olaps[l.regs, l.regs]
    g.olaps <- graph.adjacency(l.olaps, mode="undirected", weighted="weight", diag=FALSE)

    if (vis) {
        rdp <- RedPort()
        calld(rdp)
        addGraph(rdp, g.olaps)
    }

    return(g.olaps)
}


#write the annotated gene IDs to be used so can do the reduction in python
load("data/annotation.RData")
load("data/g.MI.RData")
write.table(annotation$probeID, file="data/annotation.dat", quote=FALSE, row.names=FALSE,
            col.names=FALSE)
expr.data <- read.table(file="data/disc_set/discovery_ExpressionMatrix_red.txt",
                            header=T, comment.char="", row.names=1)
tfs <- annotation[annotation$is.TF == TRUE, 1]
