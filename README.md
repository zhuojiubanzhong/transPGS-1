# transPGS：Polygenic prediction for underrepresented populations through transfer learning by utilizing shared genetic similarity shared with European populations

# Introduction
Current GWASs have been predominantly conducted in individuals of European (EUR) ancestry, with 94.6% in the EUR population and only 3.7% in the East Asian (EAS) population and 0.2% in the African (AFR) population.
Due to this underrepresentation, the performance of PGS behaves poorly in non-EUR populations, particularly in populations of AFR ancestry. For example, the PGS accuracy reduced approximately 78% across multiple traits in individuals of AFR ancestry
relative to those of EUR ancestry; similarly, the accuracy of PGS across traits was on average 40% lower in individuals of South Asian ancestry and 5% lower in individuals of EAS ancestry compared to that in those of EUR ancestry. The poor transferability 
of PGS derived from EUR ancestry data to non-EUR populations leads to great concern in health disparities . Therefore, there is an urgent need to develop novel PGS methods which can exploit data across diverse populations to better perform genetic risk prediction.

Increasing sample sizes in non-EUR GWASs for the understanding of genetic architecture underlying complex phenotypes is a necessary road for understudied populations such as EAS and AFR; but this requires plenty of expense and time. 
Alternatively, integrating existing knowledge available from EURs into non-EURs by novel approaches may be another promising strategy to improve the portability of PGS. Actually, there is a deal of evidence that significant genetic
similarity exists between the EUR and non-EUR populations at both SNP and gene levels. Such genetic similarity provides theoretical and biological support for trans-ethnic leveraging of EUR information into non-EUR studies.
Currently, there are a range of trans-ethnic statistical methods that help enhance the transferability of PGS across distinct ancestral groups; however, how to optimally integrate EUR information into non-EUR genetic research remains unknown.

In the machine learning filed, transfer learning is recognized as a novel technique that enables the utilization of knowledge acquired from auxiliary samples to enhance learning capability in target samples, which are distinct yet 
related to the former (Bastani, 2021; Li, et al., 2023; Lu, et al., 2023; Pan and Yang, 2009; Tian, et al., 2022; Zhao, et al., 2022). By borrowing the idea of transfer learning,we propose transPGS to boost the genetic prediction accuracy 
in non-EUR populations by leveraging trans-ethnic genetic similarity shared with the EUR population. Theoretically, transPGS helps capture genetic information across diverse ancestral populations and renders the prediction more efficiently and accurately.

# Example
```ruby
library(data.table)
library(SKAT)
library(glmnet)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
sourceCpp("lmm_PXEM.cpp")
source("Individual_level_transPGS.R")
source("Summary_level_transPGS.R")

################ Individual-level transPGS #######################

###### data1 is the phenotypic data for the target population, including outcome and covariates.
###### data2 is the phenotypic data for the auxiliary population, including outcome and covariates.
###### G1 is the target population genotype data (matched to 1000 Genomes Project).
###### G2 is the auxiliary population genotype data (matched to 1000 Genomes Project).

data1 <- data.frame(fread("data1.txt"))
data2 <- data.frame(fread("data2.txt"))
G1 <- data.frame(fread("geno1.txt"))
G2 <- data.frame(fread("geno2.txt"))
a1 <- Individual_level_transPGS(data1, data2,G1,G2)
head(a1)
original_beta  tl_beta
1      2.106475 2.082117
2      1.961041 1.909088
3      1.160025 1.239434
4      2.407902 2.368119
5      1.800083 1.832781
6      2.318421 2.314352

#### original_beta is the effect before transfer learning
#### tl_beta is the effect after transfer learning  

################ Summary_level transPGS #######################

###### T is the GWAS summary statistics for the target and auxiliary populations, including marginal effects as well as standard errors.
###### G1 is the target population genotype data (matched to 1000 Genomes Project).
###### G2 is the auxiliary population genotype data (matched to 1000 Genomes Project).

T <- data.frame(fread("data.txt"))
G1 <- data.frame(fread("target_geno.txt"))
G2 <- data.frame(fread("auxiliary_geno.txt"))
a2 <- Summary_level_transPGS(T,G1,G2)
head(a2)
original_beta      tl_beta
1   -0.05005137  0.015091193
2    0.09789231  0.021784294
3   -0.08230863 -0.009494420
4   -0.05116747 -0.006389153
5    0.26905642  0.007526012
6    0.02256624 -0.008591966

#### original_beta is the joint effect before transfer learning
#### tl_beta is the joint effect after transfer learning      
```
  
# Cite
Yiyang Zhu<sup>$</sup>, Wenying Chen<sup>$</sup> , Kexuan Zhu and Ping Zeng<sup>#</sup> (2024). Polygenic prediction for underrepresented populations through transfer learning by utilizing shared genetic similarity shared with European populations.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.
