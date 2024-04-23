###################################
######data1 is the phenotypic data for the target population, including outcome and covariates.
######data1 is the phenotypic data for the auxiliary population, including outcome and covariates.
######G1 is the target population genotype data (matched to 1000 Genomes Project).
######G2 is the auxiliary population genotype data (matched to 1000 Genomes Project).

Individual_level_transPGS <- function(data1, data2,G1,G2) {
###############################auxiliary population
ph1 <- as.matrix(data1[,1])
cov1 <- as.matrix(data1[,-1])
G1 <- as.matrix(G1)
fit1 <- lmm_PXEM(y=ph1, X=cov1, G=G1, PXEM=TRUE, maxIter=1000)
beta1<-as.data.frame(fit1$mub)

###############################target population
ph2 <- as.matrix(data2[,1])
cov2 <- as.matrix(data2[,-1])
G2 <- as.matrix(G2)
fit2 <- lmm_PXEM(y=ph2, X=cov2, G=G2, PXEM=TRUE, maxIter=1000)
beta2<-as.data.frame(fit2$mub)

##############transfer learning
ph2 <- as.matrix(data2[,1])
cov2 <- as.matrix(data2[,-1])
G2 <- as.matrix(G2)

beta1 <- as.matrix(beta1)
PGS <- G2 %*% beta1

G3 <- cbind(PGS,G2)

fit3 <- lmm_PXEM(y=ph2, X=cov2, G=G3, PXEM=TRUE, maxIter=1000)
coeff<-as.data.frame(fit3$mub)

O <- NULL
P <- NA
w=coeff[1,]
for (i in 1:dim(beta1)[1]) {
  O <- (w * beta1[i, 1] + coeff[i + 1, "V1"])
  P <- rbind(P, O)
}
P <- as.data.frame(P)
betaw <- as.matrix(P[-1,])
result <- cbind(beta2, betaw)
colnames(result) <- c("original_beta","tl_beta")
return(result)
}
