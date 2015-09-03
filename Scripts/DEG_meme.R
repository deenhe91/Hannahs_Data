getwd()
g.het <- read.table("g.het.txt", sep = "\t")
dim(g.het)

##

g.het$newcol <- 1

for (i in 1:223) {
  g.het[i,2] <- as.character(g.het[i,1])
}

g.het[,1] <- 1


hetdeg <- Reduce(intersect, list(gene_names, as.character(g.het[,2])))

matchlist_het <- 0

for (i in 1:168) {
  
  m <- which(g.het[,2] == hetdeg[i])
  
  matchlist_het <- c(matchlist_het, m)
  
}

##
g.ho <- read.table("g.ho.txt", sep = "\t")
dim(g.ho)

g.ho$newcol <- 1

for (i in 1:757) {
  g.ho[i,2] <- as.character(g.ho[i,1])
}

g.ho[,1] <- 2

##
g.tau <- read.table("g.tau.txt", sep = "\t")
dim(g.tau)

g.tau$newcol <- 1

for (i in 1:1146) {
  g.tau[i,2] <- as.character(g.tau[i,1])
}

g.tau[,1] <- 3

##

g.tas <- read.table("g.tas.txt", sep = "\t")
dim(g.tas)

g.tas$newcol <- 1

for (i in 1:47) {
  g.tas[i,2] <- as.character(g.tas[i,1])
}

g.tas[,1] <- 4

##

g.tpm <- read.table("g.tpm.txt", sep = "\t")
dim(g.tpm)

g.tas$newcol <- 1

for (i in 1:303) {
  g.tpm[i,2] <- as.character(g.tpm[i,1])
}

g.tpm[,1] <- 5


##concatenate tables?

degs <- c(g.het[,2], g.ho[,2], g.tau[,2], g.tas[,2], g.tpm[,2])
udegs <- unique(degs)

commondegs <- Reduce(intersect, list(gene_names, udegs))

matchlist_degs <- 0

for (i in 1:1284) {
  
  m <- which(df[,2] == commondegs[i])
  matchlist_degs <- c(matchlist_degs, m[1])
}

#probably easier way to do this but 'easier' way looks affa complex
#so...

df <- data.frame(mat.or.vec(2029,2))

df[1:223,2] <- g.het[,2]
df[1:223,1] <- 1

df[224:980,2] <- g.ho[,2]
df[224:980,1] <- 2

df[981:1679,2] <- g.tau[rand,2]
df[981:1679,1] <- 3

df[1680:1726,2] <- g.tas[,2]
df[1680:1726,1] <- 4

df[1727:2029,2] <- g.tpm[,2]
df[1727:2029,1] <- 5

write.table(df, file = "degmeme.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



##pick out appropriate fasta files

matchlist_fasta <- 0

for (i in 1:1284) {
  
  m <- which((gene_names) == commondegs[i])
  
  matchlist_fasta <- c(matchlist_fasta, m)
}

genes_fasta_deg <- gene_names[matchlist_fasta]

#create fasta files

sink ("meme_deg.fasta")

for (i in 2:1697) {
  
  cat(">")
  cat(dscr[matchlist_fasta[i]+1])
  cat("\n")
  cat(as.character(seqs[matchlist_fasta[i],2]))
  cat("\n")
  
}

sink()

#non DEGs for random sample in memelab to check whether sp1 comes up always, or is significant in DEGs

matchlist_genes <- 0

for (i in 1:1789){
  
  m = which(genes == udegs[i])
  matchlist_genes <- c(matchlist_genes, m)
  
}

non_degs <-  genes[-(matchlist_genes)]

non_deg_rand <- sample(10832,500)

df_rand <- data.frame(mat.or.vec(2000,2))
df_rand[1:500, 2] <- genes[non_deg_rand]
df_rand[1:500, 1] <- 1

df_rand[501:1000, 2] <- genes[non_deg_rand] 
df_rand[501:1000, 1] <- 2

df_rand[1001:1500, 2] <- genes[non_deg_rand]
df_rand[1001:1500, 1] <- 3

df_rand[1501:2000, 2] <- genes[non_deg_rand]
df_rand[1501:2000, 1] <- 4
