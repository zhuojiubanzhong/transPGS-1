###################################
######T is the GWAS summary statistics for the target and auxiliary populations, including marginal effects as well as standard errors.
######G1 is the target population genotype data (matched to 1000 Genomes Project).
######G2 is the auxiliary population genotype data (matched to 1000 Genomes Project).

Summary_level_transPGS <- function(T, G1, G2) {
  lambda=0.95
  m2 <- dim(G1)[2]
  r3_3 <- as.matrix(cbind(T$tar_beta, T$tar_se))
  R2 <- lambda * as.matrix(cor(G1)) + diag(1 - lambda , nrow = m2, ncol = m2)
  BETA2 <- r3_3[, 1]
  SE2 <- r3_3[, 2]
  S2 <- diag(SE2)
  SRS2 <- S2 %*% R2 %*% S2
  eig2 <- eigen(SRS2)
  d2 <- eig2$values
  U2 <- eig2$vectors
  indexd2 <- which(d2 < 0)
  if (length(indexd2) > 0) d2[indexd2] <- 1
  d22 <- 1 / sqrt(d2)
  if (length(indexd2) > 0) d22[indexd2] <- 0
  SRS122 <- U2 %*% diag(d22) %*% t(U2)
  y2 <- c(SRS122 %*% BETA2)
  Z2 <- as.matrix(SRS122 %*% (S2 %*% R2 %*% diag(1 / SE2)))
  y2 <- as.matrix(y2)
  Z2 <- as.matrix(Z2)
  G22 <- matrix(rep(1, dim(y2)[1], ncol = 1))
  fit2 <- lmm_PXEM(y = y2, X = Z2, G = G22, PXEM = TRUE, maxIter = 1000)
  betax2 <- as.data.frame(fit2$alpha)######target population joint beta

  m2 <- dim(G2)[2]
  R2 <- lambda  * as.matrix(cor(G2)) + diag(1 - lambda , nrow = m2, ncol = m2)
  r3_3 <- as.matrix(cbind(T$aux_beta, T$aux_se))
  BETA2 <- r3_3[, 1]
  SE2 <- r3_3[, 2]
  S2 <- diag(SE2)
  SRS2 <- S2 %*% R2 %*% S2
  eig2 <- eigen(SRS2)
d2 <- eig2$values
U2 <- eig2$vectors
indexd2 <- which(d2 < 0)
if (length(indexd2) > 0) d2[indexd2] <- 1
d22 <- 1 / sqrt(d2)
if (length(indexd2) > 0) d22[indexd2] <- 0
SRS122 <- U2 %*% diag(d22) %*% t(U2)
y1 <- c(SRS122 %*% BETA2)
Z1 <- as.matrix(SRS122 %*% (S2 %*% R2 %*% diag(1 / SE2)))
y1 <- as.matrix(y1)
Z1 <- as.matrix(Z1)
G11 <- matrix(rep(1, dim(y1)[1], ncol = 1))
fit1 <- lmm_PXEM(y = y1, X = Z1, G = G11, PXEM = TRUE, maxIter = 1000)
betax1 <- as.data.frame(fit1$alpha)
betax1 <- as.matrix(betax1)######auxiliary population joint beta


X <- Z2 %*% betax1
X1 <- cbind(X, Z2)
fitw <- cv.glmnet(X1, y2, intercept = FALSE, standardize = FALSE, alpha = 0, nlambda = 100, family = "gaussian", nfolds = 10)
coeff <- as.matrix(coef(fitw, s = "lambda.min"))
coeff <- as.matrix(coeff)
w <- coeff[2]######beta  after transfer learning 
O <- NULL
P <- NA
for (i in 1:dim(betax1)[1]) {
  O <- (w * betax1[i, 1] + coeff[i + 2, "1"])
  P <- rbind(P, O)
}
P <- as.data.frame(P)
beta1 <- as.matrix(P[-1,])
result <- cbind(betax2, beta1)
colnames(result) <- c("original_beta","tl_beta")
return(result)
}