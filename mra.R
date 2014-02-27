mraAnalysis <- function(g, exp="Exp1", cont="E2FGF10", p.cutoff=0.05) {
  adj.mat  <- as.matrix(get.adjacency(g))
  tfs <- which(get.vertex.attribute(g.mi, "type") == FALSE)
  badj.mat  <- adj.mat[tfs, ]
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
res <- vector(mode = "list", length = 3)
contrasts  <- c("E2FGF10", "E2AP20187", "TetE2FGF10")
for(i in 1:3) {
  print(i)
  hgt.result <- mraAnalysis(g.pcor, exp=paste("Exp", i, sep=""), cont=contrasts[i])
  res[[i]]  <- hgt.result
}

fg.ids  <- lapply(res, rownames)
tfs.s  <- Reduce(intersect, fg.ids)
annotation[which(annotation$probeID %in% tfs.s)]







