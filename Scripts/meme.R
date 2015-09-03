
library(CHNOSZ)
file <- system.file("extdata/fasta/EF-Tu.aln", package="CHNOSZ")

# Function
ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

#

setwd("/Users/hannahdeen/Desktop")

#

seqs<-ReadFasta("mm9_feP4k.fasta")

descs <- as.character(seqs[,1])

gene_names <- mat.or.vec(dim(seqs)[1],1) #gene IDs from whole genome fasta

seq_descr <- "y" #geneID and descriptive strip for fasta file

for (i in 1:dim(seqs)[1]){
  
  tmo <- strsplit(descs[i], " ")[[1]][3:30]
  
  tmp <- strsplit(descs[i], " ")[[1]][2]
  
  tmq <- strsplit(tmp, "_")[[1]][1]
  
  tmr <- c(tmq, tmo)
  
  tms <- paste0(tmr[1:30], collapse = " ")
  
  tmp_length <- nchar(tmp)
  
  seq_descr <- c(seq_descr, tms)
  
  gene_names[i] <- substr(tmp,1,(tmp_length -2 ))
  
}

gene_names <- toupper(gene_names)

##

gene_names <- "x"

for (i in 1:1696) {
  
  m <- strsplit(as.character(memes[i,1]), " ")
  gene_names <- c(gene_names, m[[1]][1])
}


dscr <- "z"

for (i in 2:21240) {
  
  klm <- strsplit(seq_descr[i], "NA")[[1]][1]
  
  dscr <- c(dscr, klm)
}

dscr <- toupper(dscr) ##seq_descr without NAs and capitalised to match run3 genes

  
##all genes in all clusters

run3 <- read.table("run3.txt", sep="\t")

run3genes <- as.character(unique(run3[,2]))


##find genes that are NOT in the fasta file. commongenes_50 and 

commongenes <- Reduce(intersect, list(gene_names, run3genes))
##generates 493 genes, no duplicates

matchlist <- 0

for (i in 1:994) {
  
  m <- which((gene_names) == run3genes[i])
  
  matchlist <- c(matchlist, m[1])
}

gene_names[matchlist]

## matchlist generates 631

meme_seqs <- seqs[matchlist,1:2]

##matching first 50 genes


#find out positions of common genes in meme50genes and then delete these.
matchlist_50 <- 0

for (i in 1:656) {
  
  m <- which((gene_names) == meme50genes[i])
  
  matchlist_50 <- c(matchlist_50, m+1)
}

meme_seqs_50 <- seqs[matchlist_50,]

##
##meme50 <- read.table("meme50.txt")
##meme50genes <- unique(as.character(meme50[,2]))
##commongenes_50 <- Reduce(intersect, list(gene_names, meme50genes))


j = 0

sink("wigwam.fasta")

for (i in 2:495) {
  
  j <- j+1
  cat(">")
  cat(dscr[matchlist[i]+1])
  cat("\n")
  cat(as.character(seqs[matchlist[i]+1,2]))
  cat("\n")
  
}

sink()

#### for first 50 genes
j = 0

sink ("meme_50.fasta")

for (i in 2:420) { #matchlist[1] is zero
  
 j <- j+1
  cat(">")
  cat(dscr[matchlist_50[i]+1])
  cat("\n")
  cat(as.character(meme_seqs_50[(j),2]))
  cat("\n")
  
}

sink()

##modmeme <- strsplit(as.character(memegenes[,1]), "_")

meme.table <- mat.or.vec(5542,2)

meme.table <- data.frame(meme.table)

vecs <- "x"

for (i in 1:5542) {
  
  vec <- strsplit(as.character(run3[i,1]), "_")
  vecb <- strsplit(vec[[1]][1], "e")
  meme.table[i,1] <- vecb[[1]][2]
  meme.table[i,2] <- as.character(run3[i,2])
  vecs <- c(vecs, vecb[[1]][2])
}


j = 0

sink ("mast.fasta")

for (i in 2:21239) { #matchlist[1] is zero
  
  j <- j+1
  cat(">")
  cat(dscr[i])
  cat("\n")
  cat(as.character(seqs[j,2]))
  cat("\n")
  
}

sink()