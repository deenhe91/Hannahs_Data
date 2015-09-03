

##must be unique sequences for mast

dscru <- unique(dscr)

## how do i find corresponding seqs?

matchlist_mast <- 0

for (i in 1:17562){
  
  m = which(dscr == dscru[i])
  matchlist_mast <- c(matchlist_mast, m[1])

}

mastseqs <-  seqs[matchlist_mast, 2]


##

j = 0

sink ("mast.fasta")

for (i in 2:17562) { #matchlist[1] is zero
  
  j <- j+1
  cat(">")
  cat(dscru[i])
  cat("\n")
  cat(as.character(mastseqs[j]))
  cat("\n")
  
}

sink()


## pull out genes from mast output

mast_output <- mast_output[8:3110,]

mast_genes <- mast_genes[2:3102]  #get rid of first z and last bracket

mast_genes <- "z"

for (i in 2:3103) {
  tmp <- strsplit(as.character(mast_output[i]), " ")
  tmp2 <- tmp[[1]][1]
  mast_genes <- c(mast_genes, tmp2)
}


mastagenes <- "z"

for (i in 1:17561) {
  tmp <- strsplit(as.character(mast_tmp[i]), " ")
  tmp2 <- tmp[[1]][1]
  mastagenes <- c(mastagenes, tmp2)
}


##overlap with DEGs?

common <- Reduce(intersect, list(udegs, mast_genes))


##

common_tau <- Reduce(intersect, list(g.tau[,2], mast_genes))
                     
matchlist_tau <- 0

  for (i in 1:1146){
  
  m <- which(gene_names == g.tau[i,2])
  
  matchlist_tau <- c(matchlist_tau, m)
  
}
  
nontau <- gene_names[-(matchlist_tau)]


# length decreases every time loop runs. 
#
tmp_match <- 0

for (i in 1:206) {
  
 m <- which(g.het[i,2] == common_het)
 tmp_match <- c(tmp_match, m)
 
}
