mraAnalysis <- function(g, exp="Exp1", cont="E2FGF10", p.cutoff=0.05) {
  adj.mat  <- as.matrix(get.adjacency(g))
  tfs <- which(get.vertex.attribute(g.mi, "type") == FALSE)
  genes  <- setdiff(1:ncol(adj.mat), tfs)
  badj.mat  <- adj.mat[tfs, genes]
  universe  <- colnames(adj.mat)
  sigt  <- Fletcher2013pipeline.deg(what=exp)
  hits  <- get(cont, sigt)
  hits  <- hits[which(hits %in% universe)]
  
  gene.sets <- apply(badj.mat, 1, function(x) { names(which(x == 1))})
  hgt.results  <- as.data.frame(multiHyperGeoTest(gene.sets, universe, hits))
  sign.tfs  <- hgt.results[which(hgt.results$Pvalue < p.cutoff), ]
    
  return(sign.tfs)
}

#correlation network
multTest  <- function(g, p.cutoff) {
  res <- vector(mode = "list", length = 3)
  contrasts  <- c("E2FGF10", "E2AP20187", "TetE2FGF10")
  for(i in 1:3) {
    print(i)
    hgt.result <- mraAnalysis(g, exp=paste("Exp", i, sep=""), cont=contrasts[i],
                              p.cutoff=p.cutoff)
    res[[i]]  <- hgt.result
  }
  
  return(res)
}

getOlapsExp  <- function(res) {
  fg.ids  <- lapply(res, rownames)
  tfs.s  <- Reduce(intersect, fg.ids)
  
  return(tfs.s)
}

annotation[which(annotation$probeID %in% tfs.s)]


all.tfs  <- unique(c(tfs.cor, tfs.mi.dpi, tfs.mi, tfs.pcor))
tfs.list  <- list(tfs.cor, tfs.pcor, tfs.mi, tfs.mi.dpi)

m  <- t(sapply(all.tfs, function(tf) {sapply(tfs.list, function(l){tf %in% l})}))
sm  <- apply(m, 1, function(x){table(x)["TRUE"]})
sm  <- sm[sm >= 2]
tfs.s  <- names(sm)
tfs.symbol  <- annotation[annotation$probeID %in% tfs.s, 2]
tfs.table  <- cbind(tfs.symbol, m[tfs.s, ])
colnames(tfs.table) <- c("Name", "Cor", "PCor", "MI", "MI+DPI")






