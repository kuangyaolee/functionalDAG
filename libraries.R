############################################################### 
#                        Run libraries                        #
###############################################################
library(graph)
library(Rgraphviz)
library(pcalg)
library(e1071)
library(mvtnorm)
library(kernlab)
library(MASS)
library(igraph)
library(corpcor)
library(CAM)
library(fda)

###############################################################
# Calculate the gamma parameter for Gaussian Kernel
###############################################################
# Input :
#         A: Input matrix (n by p)

# Output:
#         gam : the value of gamma
###############################################################
CalGam=function(A){
  if(is.vector(A)){
    n = length(A)}
  else n = dim(A)[1]
  tmp=rowSums(as.matrix(A*A))%*%t(rep(1,n))
  K=tmp+t(tmp)-2*A%*%t(A)
  K=K*(K>=0)
  tou=as.double(sum(sqrt(K))/(n*(n-1)))
  gam=1/(2*tou^2)
  return(gam)
}


###############################################################
# Calculate the Gaussian Kernel matrix
###############################################################
# Input :
#         A: Input matrix (n by p)
#         B: Bases matrix (m by p, m can be different from n)

# Output:
#         K : Gaussian kernel matrix (m by n) 
###############################################################
KGaussian=function(gamma,A,B){
  if(is.vector(A)){
    n = length(A)}
  else n = dim(A)[1]
  if(is.vector(B)){
    m = length(B)}
  else m = dim(B)[1]
  tmp_1=rowSums(as.matrix(A*A))%*%matrix(1,1,m)
  tmp_2=rowSums(as.matrix(B*B))%*%matrix(1,1,n)
  K=tmp_1+t(tmp_2)-2*A%*%t(B)
  K=exp(-K*gamma)
  return(K)
}


###############################################################
# Calculate the Brownian Kernel matrix
###############################################################
# Brownian kernel, i.e. k(a,b) = min(a,b)
# Inputs %
# A: m x 1 %
# B: n x 1 %

# Outputs %
# K: m x n Gram matrix %
# by Kuang-Yao Lee
###############################################################

KBrownian = function(A, B)
{
  m = length(A);
  n = length(B);
  
  AA = A %*% t(rep(1,n))
  BB = t(B %*% t(rep(1,m)))
  
  K = (AA-BB>0)*BB + (AA-BB<=0)*AA
  return(K)
}



###############################################################
# %% find the lower dimension approximation of K, i.e. reduced kernel
# % take \lambda_k if \lambda_k/\lambda_1>C
# % Status =1 if reducing the space in L_2; =0 if in RKHS.
###############################################################
RedKer3<-function(K, NumOfDim)
{
  eig_K=eigen(K);d = eig_K$values;U=eig_K$vectors
  COEF = U[,1:NumOfDim]
  SCORE = COEF %*% diag(sqrt(d[1:NumOfDim]));d_max=d[1];
  VALUES = d[1:NumOfDim]
  return(list(SCORE=SCORE,VALUES =VALUES))
}

############################################################### 
# function: Q = I - J/n                                   
############################################################### 
qmat = function(n) return(diag(n)-rep(1,n)%*%t(rep(1,n))/n)

###############################################################
#          function: power of a matrix 
###############################################################
matpower = function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}



## Purpose: Compute Structural Hamming Distance between graphs g1 and g2
## ----------------------------------------------------------------------
## Arguments:
## - am1, am2: Input graphs
## (graph objects;connectivity matrix where m[x,y]=1 iff x->1
## and m[x,y]=m[y,x]=1 iff x-y; pcAlgo-objects)
## ----------------------------------------------------------------------
## Author: Markus Kalisch, Date:  1 Dec 2006, 17:21
## inputs are adjacency matrices; modified by Kuang-Yao Lee 05-12-18

## Idea: Transform g1 into g2
shd_am <- function(am1,am2)
{
  am1[am1 != 0] <- 1
  
  am2[am2 != 0] <- 1
  
  shd <- 0
  ## Remove superfluous edges from g1
  s1 <- am1 + t(am1)
  s2 <- am2 + t(am2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  am1[ind] <- 0
  shd <- shd + length(ind)/2
  ## Add missing edges to g1
  ind <- which(ds < 0)
  am1[ind] <- am2[ind]
  shd <- shd + length(ind)/2
  ## Compare Orientation
  d <- abs(am1-am2)
  ## return
  shd + sum((d + t(d)) > 0)/2
}

## Purpose: Return TPR, FPR and TDR of comparison of two undirected graphs
## ----------------------------------------------------------------------
## Arguments:
## - gl: Estimated graph (may be directed, but the direction will
##       be dropped internally)
## - gt: True graph (may be directed, but the direction will
##       be dropped internally)
## ----------------------------------------------------------------------
## Author: Markus Kalisch, Date: 26 Jan 2006, 17:35
## inputs are adjacency matrices; modified by Kuang-Yao Lee 05-12-18

## When 'graph' returns the 'weight matrix' again:
##   ml <- as(ugraph(gl), "matrix")
##   mt <- as(ugraph(gt), "matrix")
##   p <- dim(ml)[1]
compareGraphs_am <- function(ml,mt) {
  p <- dim(ml)[2]
  
  mt[mt != 0] <- rep(1,sum(mt != 0))
  ml[ml != 0] <- rep(1,sum(ml != 0)) ## inserted to fix bug
  
  ## FPR :=  #{misplaced edges} / #{true gaps}
  diffm <- ml-mt
  nmbTrueGaps <- (sum(mt == 0)-p)/2
  ##  print(list(p=p,sum=sum(mt==0),mt=mt,ml=ml,nmbTrueGaps=nmbTrueGaps,diffm=diffm))
  fpr <- if (nmbTrueGaps == 0) 1 else (sum(diffm > 0)/2)/nmbTrueGaps
  
  ## TPR := #{correctly found edges} / #{true edges}
  diffm2 <- mt-ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  tpr <- if (nmbTrueEdges == 0) 0 else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
  
  ## TDR := #{correctly found edges} / #{found edges}
  trueEstEdges <- (nmbTrueEdges-sum(diffm2 > 0)/2) ## #{true edges} - #{not detected}
  tdr <-
    if (sum(ml == 1) == 0) { ## no edges detected
      if (trueEstEdges == 0) 1 ## no edges in true graph
      else 0
    } else trueEstEdges/(sum(ml == 1)/2)
  
  ## return named vector:
  c(tpr = tpr, fpr = fpr, tdr = tdr)
}

############################################################### 
# z-transformation
############################################################### 
z_trans <- function(r) .5 * log((1+r)/(1-r))





################################################################
# The mercer decomposition of the Brownian kernel
################################################################
Brownian2 <- function(tt){
  K_T = matrix(0,length(tt),100)
  J=1:100
  for (i in 1:length(tt)){
    K_T[i,]=   sqrt(2) * sin((J-1/2)*pi*tt[i])
  }
  return(K_T)
}

################################################################
# The (firt five) fourier transformation function
################################################################
fourier <- function(tt){
  K_T = matrix(0,length(tt),5)
  J=1:5
  for (i in 1:length(tt)){
    K_T[i,]= c(1, sqrt(2) * sin(2 * pi * tt[i]), sqrt(2) * cos(2 * pi * tt[i]), sqrt(2) * sin(4 * pi * tt[i]), sqrt(2) * cos(2 * pi * tt[i]) )
  }
  return(K_T)
}



###############################################################
# Compute the coordinates of functions (version 2, updated 07-23-18)
# added compute gamma_T when using Gaussian kernel, and unbalanced design
###############################################################
# % %
#   % Inputs %
#   % f_t: N x 1 f(T)%
#  % t: N x 1 time points T %
#   % ind_ind: N x 1 indicator for individuals
# % NumofBasis: the dimension of \Omega \hi n
# eps_T: eidge parameter (only needed for balanced design)
# 
# % %
#   % Outputs %
# % XX:  n x NumofBasis; the coordinates for each function
# % by Kuang-Yao Lee
###############################################################

CalCoor2<-function(f_t, ind_ind,tt,kernel_T, NumofBasis,eps_T,UNB,K_T_SCORE)
{   
  n = length( unique(ind_ind) )
  m = length(tt)/n
  f_t_tmp = t(array(f_t, c(m,n)))
  if(UNB==0){
    tt = unique(tt)
    if (kernel_T ==1){ 
      gamma_T=CalGam(tt)
      K_T = KGaussian (gamma_T, tt,tt);
    }else{
      K_T = KBrownian (tt, tt);
    }
    
    SCORE = RedKer3(K_T, NumofBasis)
    
    KK_T = SCORE$SCORE #dim(K_T)[1] x NumofBasis
    
    XX = t( diag( 1/(SCORE$VALUES + eps_T * SCORE$VALUES[1]) ) %*% t(KK_T) %*% t(f_t_tmp) )
  }else{
    XX = matrix(0,n,NumofBasis)
    
    for(i in 1:n){
      XX_tmp = matpower(t(K_T_SCORE$SCORE[ind_ind==i,]) %*% K_T_SCORE$SCORE[ind_ind==i,] + eps_T * K_T_SCORE$VALUES[1] * diag(NumofBasis),-1) %*% t(K_T_SCORE$SCORE[ind_ind==i,]) %*% f_t_tmp[i,]
      XX[i,]=XX_tmp
    }
    
  }
  
  
  return(list(XX=XX))
}




############################################################### 
# test CI using the H-S norm of fPCO 
# will report the pvalue based on normal approximation 
# 07-23-18 by Kuang-Yao Lee
############################################################### 
CalfPCO<- function(I,J,S,suffStat)
{
  GX_orth =suffStat$GX_orth; epsilon_1 = suffStat$epsilon_1; 
  epsilon_2 = suffStat$epsilon_2; NumOfDim = suffStat$NumOfDim
  
  n=dim(GX_orth)[1]; ptimesNumOfDim=dim(GX_orth)[2];p=ptimesNumOfDim/NumOfDim
  
  Xu=GX_orth[,(1+(I-1)*NumOfDim):(I*NumOfDim) ];Xv=GX_orth[,(1+(J-1)*NumOfDim):(J*NumOfDim) ]
  
  S_ind_tmp =array(1:ptimesNumOfDim,c(NumOfDim,p))
  S_ind = array(S_ind_tmp[,S],c(1,length(S)*NumOfDim))
  Gw = GX_orth[,S_ind]
  
  if (length(S)==0){
    Pw = matrix(0,n,n)
  }else{
    Eig= eigen(Gw %*% t(Gw)); temp_U = Eig$vectors; temp_D=Eig$values;eig_1=max(temp_D);
    Pw=temp_U %*% diag( (temp_D /(temp_D+eig_1*epsilon_1))^2 ) %*% t(temp_U)
  }
  
  SigUUW= t(Xu) %*% ( diag(n) - Pw ) %*% Xu
  SigUVW= t(Xu) %*% ( diag(n) - Pw ) %*% Xv
  SigVVW= t(Xv) %*% ( diag(n) - Pw ) %*% Xv
  
  eig_2 = max( c( as.double(eigen(SigUUW, only.values = TRUE)$values[1]), as.double(eigen(SigVVW, only.values = TRUE)$values[1]) ) )
  
  if(NumOfDim==1){
    fPCO = SigUUW^(-.5) * SigUVW * SigVVW^(-.5)
  }else{
    fPCO = matpower(SigUUW + epsilon_2*eig_2 * diag(NumOfDim), -.5) %*% SigUVW %*% matpower(SigVVW + epsilon_2*eig_2 * diag(NumOfDim), -.5)
  }
  # fPCOnorm = as.double((eigen(fPCO %*% t(fPCO))$values[1])^.5) #operator norm
  fPCOnorm = as.double(sum (fPCO^2)^.5) #H-S norm
  
  # based on normal approximation
  test_stat = sqrt(n- length(S) - 3) * fPCOnorm
  # test_stat_tmp = sqrt(n- length(S_ind) - 3) * z_trans(fPCO) # z-transformation
  # test_stat = sum(test_stat_tmp^2)
  
  if (is.na(test_stat)) 0 else test_stat
  p_val = 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  # p_val = pchisq(test_stat, df = NumOfDim^2, lower.tail = FALSE)
  
  return(p_val)
}

