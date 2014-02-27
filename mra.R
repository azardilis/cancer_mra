library(HTSanalyzeR)

overlapTest  <- function(badj.mat, hits) {
  univ  <- colnames(badj.mat)
  gene.sets <- apply(badj.mat, 1, function(x) { names(which(x == 1))})
  
  hgt.results  <- multiHyperGeoTest(gene.sets, univ, hits)
  
  return(hgt.results)
}

sigt  <- Fletcher2013pipeline.deg(what="Exp1")
hits  <- sigt$E2FGF10

hgt.results  <- overlapTest(badj.mat, hits)
top.tfs  <- names(sort(hgt.results[, 'Observed Hits'], decreasing=T)[1:5])
#get the names from annotation
