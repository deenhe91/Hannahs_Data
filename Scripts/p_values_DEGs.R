
#load mouseac data

getwd()

setwd("/Users/hannahdeen/Dropbox/Hannah Deen")

list.files()

setwd("Data")

setwd("mouseac")

mouseac <- read.csv("mouseac.csv")




#rearrange order for comparisons

mouseac[,"MODELDIS"]

mouseac[,"MODELDIS"] <- factor(mouseac[, "MODELDIS"], c("WILD", "HET_TASTPM", "HO_TASTPM", "TAS10", "TAU", "TPM"))



#regression analysis of first gene

fit <- lm(X0610005K03RIK ~ factor(AGE_cat) * REGION * MODELDIS , data=mouseac)

summary(fit)

p.val<-summary(fit)$coef[,4]




#making p value table (p.val.table)

p.val.table <- mat.or.vec(12570,72)

colnames(p.val.table) <- names(p.val)

p.val.table <- data.frame(p.val.table)

colnames(mouseac)[8:12577]

genes <- colnames(mouseac)[8:12577]

rownames(p.val.table) <- genes



# saving p.val.table
#first changed setwd to appropriate folder, then:
# write.table(p.val.table, file = "p.val.table.csv" , sep="," )

#q value table (q.val.table)
#how can loop be used to create q values?




q.val.table <- mat.or.vec(12570,72)

q.val.table <- data.frame(q.val.table)

colnames(q.val.table) <- names(p.val)

colnames(mouseac)[8:12577]

genes <- colnames(mouseac)[8:12577]

rownames(q.val.table) <- genes



sig.val.table <- mat.or.vec(12570,72)

sig.val.table <- data.frame(sig.val.table)

colnames(sig.val.table) <- names(p.val)

rownames(sig.val.table) <- genes


# filling tables with: 
#p values from regression analysis (p.val.table)
#q.vals (BH multiple correction) for q.val.table
#sig.vals (q.val < 0.1) for sig.val.table



row <- 0

for (gene.expr in mouseac[,8:12577]) {  
  
  fit <- lm(gene.expr ~ factor(mouseac[,4]) * mouseac[,3] * mouseac[,5])
  p.val<-summary(fit)$coef[,4]
  
  row <- row + 1
  
  p.val.table[row,] <- p.val
  
  q.val <- p.adjust(p.val, method = "BH")
  
  q.val.table[row,] <- q.val
  
  sig.vals <- q.val < 0.1
  
  sig.val.table[row,]<-sig.vals
  
}

setwd("/Users/hannahdeen/Dropbox/Hannah Deen/excel tables")
TPM.table <- mat.or.vec(12570,8)
TPM.table <- data.frame(TPM.table)
TPM.cols <- c(11, 30, 31, 32, 42, 70, 71,72)
colnames(TPM.table) <- names(p.val.table[,TPM.cols])
rownames(TPM.table) <- genes

for (i in 1:8) {
TPM.table[,i] <- sig.val.table[,TPM.cols[i]]
}
TPM.table$newcol <- apply(TPM.table, 1, any)
tmp <- which(TPM.table[,13] == TRUE)
sig.genes_tpm <- rownames(TPM.table[tmp,])
length(sig.genes_tpm)

write.table(p.val.table, file = "p.val.table.csv", sep=",", col.names=NA)
write.table(q.val.table, file = "q.val.table.csv", sep=",", col.names=NA)
write.table(sig.val.table, file = "sig.val.table.csv", sep=",", col.names=NA)
#col.names=NA

setwd("Users/hannahdeen/Dropbox/Hannah Deen/Data/mouseac")


#retrieve column numbers that contain genetic background names. e.g. "TAU"
#using grep function, e.g., grep("TAU", colnames(p.val.table))


#inserting sig.vals into TAU.table

TAU.table <- mat.or.vec(12570,12)

TAU.table <- data.frame(TAU.table)

TAU.cols <- grep("TAU", colnames(p.val.table))

colnames(TAU.table) <- names(p.val.table[,TAU.cols])
rownames(TAU.table) <- genes


for (i in 1:12) {
  
  TAU.table[,i] <- sig.val.table[,TAU.cols[i]]
}


#use TAU.table$newcol <- apply(TAU.table, 1, function(row) any(TAU.table[1:12,1]==1)) 
#but in a loop?


tmp <- c(1,2,3,4,6,10,11,12) 
TAU.table$newcol <- apply(TAU.table[,tmp], 1, any)
tmp <- which(TAU.table[,13] == TRUE)
sig.genes_tau <- rownames(TAU.table[tmp,])



#generate DEG lists for each age within a genetic background


tmp <- c(1,2,6,10)
TAU.table$tau16 <- apply(TAU.table[,tmp], 1, any)
tau16 <- which(TAU.table[,14] == TRUE)
tau16 <- rownames(TAU.table[tau16,])

tmp <- c(1,3,6,11)
TAU.table$tau32 <- apply(TAU.table[,tmp], 1, any)
tau32 <- which(TAU.table[,15] == TRUE)
tau32 <- rownames(TAU.table[tau32,])

tmp <- c(1,4,6,12)
TAU.table$tau72 <- apply(TAU.table[,tmp], 1, any)
tau72 <- which(TAU.table[,16] == TRUE)
tau72 <- rownames(TAU.table[tau72,])

write(sig.genes_tau, file = "DEG_TAU.txt", sep = "\t")

  
write.table(TAU.table, file = "TAU.table4.csv", sep=",", col.names=NA)


#other genetic backgrounds

HET_TASTPM.table <- mat.or.vec(12570,12)

HET_TASTPM.table <- data.frame(HET_TASTPM.table)

HET_TASTPM.cols <- grep("HET_TASTPM", colnames(p.val.table))

colnames(HET_TASTPM.table) <- names(p.val.table[,HET_TASTPM.cols])
rownames(HET_TASTPM.table) <- genes


for (i in 1:12) {
  
  HET_TASTPM.table[,i] <- sig.val.table[,HET_TASTPM.cols[i]]
}

tmp <- c(1, 2, 3, 4, 6, 10, 11, 12)
HET_TASTPM.table$hipcol <- apply(HET_TASTPM.table[,tmp], 1, any)
tmp <- which(HET_TASTPM.table[,13] == TRUE)
sig.genes_het <- rownames(HET_TASTPM.table[tmp,])
length(sig.genes_het)

tmp <- c(1,2,6,10)
HET_TASTPM.table$het16 <- apply(HET_TASTPM.table[,tmp], 1, any)
het16 <- which(HET_TASTPM.table[,14] == TRUE)
het16 <- rownames(HET_TASTPM.table[het16,])
length(het16)

tmp <- c(1,3,6,11)
HET_TASTPM.table$het32 <- apply(HET_TASTPM.table[,tmp], 1, any)
het32 <- which(HET_TASTPM.table[,15] == TRUE)
het32 <- rownames(HET_TASTPM.table[het32,])
length(het32)

tmp <- c(1,4,6,12)
HET_TASTPM.table$het72 <- apply(HET_TASTPM.table[,tmp], 1, any)
het72 <- which(HET_TASTPM.table[,16] == TRUE)
het72 <- rownames(HET_TASTPM.table[het72,])
length(het72)


write.table(HET_TASTPM.table, file = "HET_TASTPM.table.csv", sep=",", col.names=NA)

##

HO_TASTPM.table <- mat.or.vec(12570,12)

HO_TASTPM.table <- data.frame(HO_TASTPM.table)

HO_TASTPM.cols <- grep("HO_TASTPM", colnames(p.val.table))

colnames(HO_TASTPM.table) <- names(p.val.table[,HO_TASTPM.cols])
rownames(HO_TASTPM.table) <- genes

for (i in 1:12) {
  
  HO_TASTPM.table[,i] <- sig.val.table[,HO_TASTPM.cols[i]]
}

tmp <- c(1,2,3,4,6,10,11,12) 
HO_TASTPM.table$newcol <- apply(HO_TASTPM.table[,tmp], 1, any)
tmp <- which(HO_TASTPM.table[,13] == TRUE)
sig.genes_ho <- rownames(HO_TASTPM.table[tmp,])
length(sig.genes_ho)


tmp <- c(1,2,6,10)
HO_TASTPM.table$ho16 <- apply(HO_TASTPM.table[,tmp], 1, any)
ho16 <- which(HO_TASTPM.table[,14] == TRUE)
ho16 <- rownames(HO_TASTPM.table[ho16,])
length(ho16)

tmp <- c(1,3,6,11)
HO_TASTPM.table$ho32 <- apply(HO_TASTPM.table[,tmp], 1, any)
ho32 <- which(HO_TASTPM.table[,15] == TRUE)
ho32 <- rownames(HO_TASTPM.table[ho32,])
length(ho32)

tmp <- c(1,4,6,12)
HO_TASTPM.table$ho72 <- apply(HO_TASTPM.table[,tmp], 1, any)
ho72 <- which(HO_TASTPM.table[,16] == TRUE)
ho72 <- rownames(HO_TASTPM.table[ho72,])
length(ho72)


write.table(HO_TASTPM.table, file = "HO_TASTPM.table.csv", sep=",", col.names=NA)

##

TAS10.table <- mat.or.vec(12570,12)

TAS10.table <- data.frame(TAS10.table)

TAS10.cols <- grep("TAS10", colnames(p.val.table))

colnames(TAS10.table) <- names(p.val.table[,TAS10.cols])
rownames(TAS10.table) <- genes

for (i in 1:12) {
  
  TAS10.table[,i] <- sig.val.table[,TAS10.cols[i]]
}

tmp <- c(1,2,3,4,6,10,11,12) 
TAS10.table$newcol <- apply(TAS10.table[,tmp], 1, any)
tmp <- which(TAS10.table[,13] == TRUE)
sig.genes_tas <- rownames(TAS10.table[tmp,])
length(sig.genes_tas)


tmp <- c(1,2,6,10)
TAS10.table$tas16 <- apply(TAS10.table[,tmp], 1, any)
tas16 <- which(TAS10.table[,14] == TRUE)
tas16 <- rownames(TAS10.table[tas16,])
length(tas16)

tmp <- c(1,3,6,11)
TAS10.table$tas32 <- apply(TAS10.table[,tmp], 1, any)
tas32 <- which(TAS10.table[,15] == TRUE)
tas32 <- rownames(TAS10.table[tas32,])
length(tas32)

tmp <- c(1,4,6,12)
TAS10.table$tas72 <- apply(TAS10.table[,tmp], 1, any)
tas72 <- which(TAS10.table[,16] == TRUE)
tas72 <- rownames(TAS10.table[tas72,])
length(tas72)


write.table(TAS10.table, file = "TAS10.table.csv", sep=",", col.names=NA)

##

TPM.table <- mat.or.vec(12570,12)

TPM.table <- data.frame(TPM.table)

TPM.cols <- c(11,30,31,32,41,42,67,68,69,70,71,72)

colnames(TPM.table) <- names(p.val.table[,TPM.cols])
rownames(TPM.table) <- genes


for (i in 1:12) {
  
  TPM.table[,i] <- sig.val.table[,TPM.cols[i]]
}

tmp <- c(1,2,3,4,6,10,11,12)
TPM.table$newcol <- apply(TPM.table[,tmp], 1, any)
tmp <- which(TPM.table[,13] == TRUE)
sig.genes_tpm <- rownames(TPM.table[tmp,])
length(sig.genes_tpm)


tmp <- c(1,2,6,10)
TPM.table$tpm16 <- apply(TPM.table[,tmp], 1, any)
tpm16 <- which(TPM.table[,14] == TRUE)
tpm16 <- rownames(TPM.table[tpm16,])
length(tpm16)

tmp <- c(1,3,6,11)
TPM.table$tpm32 <- apply(TPM.table[,tmp], 1, any)
tpm32 <- which(TPM.table[,15] == TRUE)
tpm32 <- rownames(TPM.table[tpm32,])
length(tpm32)

tmp <- c(1,4,6,12)
TPM.table$tpm72 <- apply(TPM.table[,tmp], 1, any)
tpm72 <- which(TPM.table[,16] == TRUE)
tpm72 <- rownames(TPM.table[tpm72,])
length(tpm72)


write.table(TPM.table, file = "TPM.table.csv", sep=",", col.names=NA)

###
#extract significant differentially expressed genes from tables?











list = 0

for (i in 1:279) {
  
  m = which(het32 == same_het[i])
  list <- c(list, m)
  
}


clist <- 0

for (i in 1:10) {
  
  c <- which(commons == commoncommons[i])
  clist <- c(clist, c)
}

