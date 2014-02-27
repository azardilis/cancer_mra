ndReg <- function(cor, beta=0.5, alpha=0.1) {

    n.tf <- nrow(cor)
    n <- ncol(cor)

    #removing self-loops
    diag(cor) <- 0

    #making TF-TF symmetric
    tf.net <- cor[1:n.tf, 1:n.tf]
    xy <- which((tf.net != t(tf.net)) != 0, arr.ind=T)

    tf.net.final <- tf.net
    if (dim(xy)[1] != 0) {
      for (i in 1:nrow(xy)) {
        x <- xy[i, ][1]
        y <- xy[i, ][2]
        if (tf.net[x, y] != 0 & tf.net[x, y] != 0) {
          tf.net.final[x, y] <- (tf.net[x,y] + tf.net[y, x]) / 2
          tf.net.final[y, x] <- tf.net.final[x, y]
        } else if (tf.net[x, y] == 0) {
          tf.net.final[x, y] <- tf.net[y, x]
          tf.net.final[y, x] <- tf.net.final[x, y]
        } else if (tf_net[y,x] == 0) {
          tf.net.final[x, y] -> tf.net[x, y]
          tf.net.final[y, x] -> tf.net.final[x, y]
        }
        
      }  
    }
    
    cor[1:n.tf, 1:n.tf] <- tf.net.final

    y <- quantile(c(cor), 1-alpha)
    mat.th <- cor * (cor >= y)

    temp.net <- (mat.th>0)*1.0
    temp.net.remain <- (mat.th==0)*1.0
    mat.th.remain = cor*temp.net.remain
    m11=max(mat.th.remain)

    eig.res <- NULL

    #if (rcond(eig.res$vectors) < 10^-10) {
    rp <- 0.001
    rand.tf <- rp * matrix(runif(n.tf*n.tf), nrow=n.tf)
    rand.tf <- (rand.tf + t(rand.tf)) / 2
    diag(rand.tf) <- 0
    rand.target <- rp * matrix(runif(n.tf*(n-n.tf)), nrow=n.tf)
    mat.rand <- cbind(rand.tf, rand.target)
    mat.th <- mat.th + mat.rand
    mat.th <- rbind(mat.th, matrix(rep(0, (n-n.tf)*n), nrow=(n-n.tf)))
    
    eig.res <- eigen(mat.th)
    
    U <- eig.res$vectors
    D <- matrix(rep(0, nrow(mat.th)*nrow(mat.th)), nrow=nrow(mat.th))
    diag(D) <- eig.res$values

    lam.n <- abs(min(min(diag(D)), 0))
    lam.p <- abs(max(max(diag(D)), 0))

    m1 <- lam.p * (1-beta)/beta
    m2 <- lam.n * (1+beta)/beta
    scale.eigen <- max(m1, m2)

    D <- D / (scale.eigen + D)

    nmat <- U %*% D %*% solve(U)

    net.new2 <- nmat[1:n.tf, ]
    m2 <- min(net.new2)
    net.new3 <- (net.new2 + max(m11-m2,0)) * temp.net;
    mat.nd <- net.new3 + mat.th.remain

    return(mat.nd)
}

edges  <- read.table("dream/net1_expression_data_Correlation_predictions.txt")
colnames(edges)[3]  <- "weight"
g.df  <- graph.data.frame(edges)
adj.mat  <- abs(as.matrix(get.adjacency(g.df, attr="weight")))
tfs  <- as.character(read.table("dream/net1_transcription_factors.tsv")$V1)

badj.mat  <- adj.mat[tfs, ]
badj.mat.nd  <- ndReg(badj.mat)

edges.out  <- read.table("dream/net1_expression_data_Pearson-ND_predictions.txt")
colnames(edges.out)[3]  <- "weight"
g.out  <- graph.data.frame(edges.out)
adj.out  <- abs(as.matrix(get.adjacency(g.out, attr="weight")))

badj.out  <- adj.out[tfs, ]

#order it the same way as the other output so can be comparable
badj.out  <- badj.out(rownames(badj.mat.nd), colnames(badj.mat.nd))





