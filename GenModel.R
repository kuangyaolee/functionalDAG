###############################################################
# Generate DAG
###############################################################
# Input :
#         p: number of nodes
#         prob: sparsity rate
#         hub: 0 for random graph; 1 for 10 independent hubs

# Output:
#         pxp adjacent matrix
# by Kuang-Yao Lee
###############################################################
genDAG <- function (p, rate, hub)
{
  DAG=matrix(0,p,p);
  #hub
  if(hub==1){
    HUB_size=p/10;
    NumOfHubs=10;
    for (i in 1:NumOfHubs){
      DAG[(1+(i-1)*HUB_size),(2+(i-1)*HUB_size):(i*HUB_size)]=rep(1,HUB_size-1)
    } }
  else{
    # generate random graph
    ComG = upper.tri(matrix(1,p,p),diag=FALSE)*1
    IND=which(ComG==1)
    IND_DAG=sample(IND,round(p*(p-1)*rate/2),replace=FALSE);
    DAG[IND_DAG]=1
  }
  
  # moralized graph
  mDAG=DAG
  id_kids=which(colSums(DAG)>=2) #find v-structure, i.e. x->z<-y
  for (i in 1:length(id_kids) ){
    id_parents=which(DAG[,id_kids[i]]!=0)
    for (j in 2:length(id_parents) ){
      for (k in 1:(j-1)){
        mDAG[id_parents[k],id_parents[j]]=1 # connecting parents 
      }
    }
  }
  mDAG = ((mDAG + t(mDAG))>0)*1
  list(DAG = DAG, mDAG = mDAG)
}


###############################################################
# Generate random error function (for both balanced/unbalanced designs)
###############################################################
# % inputs:
#   % n: sample size
# % t_obs: (N=m*n) x 1 observed time points
# % ind_ind: (N=m*n) x 1 indicator for individuals
# % NumOfCom: #(mixture components)
#   % kernel_type: 1 for Gaussian; 2 for Brownian motion
# % kernel_par: parameters for given kernels, e.g. width for Gaussian
# % loc_type: distributions of \ti_i; 1 for Uniform(0, 1)
# % var_type: distributions of \xi_i: 1 for standard normal
# % outputs:
#   % E_t: (N=m*n) x 1 observed f_t
# % by Kuang-Yao Lee
###############################################################

GenRanFun<- function(n, t_obs,ind_ind,NumOfCom,kernel_type,kernel_par,loc_type,var_type)
{
  N = length(t_obs);
  E_t = matrix(0,N,1);
  
  
  for (i in 1:n){
    if (kernel_type ==1) {
      Xi = var_type *rnorm(NumOfCom)
      ti = loc_type *runif(NumOfCom)
      #RBF
      K_T = KGaussian (kernel_par, t_obs[ind_ind==i],ti)
      E_t[ind_ind==i] =  K_T %*% Xi
    }else
      if (kernel_type ==2){ 
        Xi = var_type *rnorm(NumOfCom)
        # ti = loc_type *runif(NumOfCom)
        ti = seq(1/NumOfCom,1,by=1/NumOfCom) # fix ti to generate Gausian process
        #Brownian
        K_T = KBrownian (t_obs[ind_ind==i], ti)
        E_t[ind_ind==i] = K_T %*% Xi
      }else
        if (kernel_type==3){
          Xi = var_type *rnorm(100)
          #Brownian eigenfunction
          K_T = Brownian2 (t_obs[ind_ind==i])
          E_t[ind_ind==i] = K_T %*% Xi
        }else
          if (kernel_type==4){
            Xi = var_type *rnorm(5)
            #five Fourier basis functions
            K_T = fourier (t_obs[ind_ind==i])
            E_t[ind_ind==i] = K_T %*% Xi
          }
  }
  
  
  
  return(E_t)
}


###############################################################
# generate synthetic models (for both balanced/unbalanced designs)
###############################################################
# % generate synthetic models
# % Inputs %
#   % n: sample size
# % p: dimensionality
#   m: number of replicates for each individual
# % model: model indicator
# % %
#   % Outputs %
#   % dat: (N=n*m) x p, f(T) for n individuals, each of which observes p nodes and m_i time points
# % t: (N=n*m) x 1 time points
# % ind_ind: (N=n*m) x 1 indicator for individuals
# % by Kuang-Yao Lee
# ###############################################################
GenMod<-function(n,p,m,model,UNB,STD,HUB=0,rate=1,kernel_type,kernel_par,loc_type,var_type){
  N = sum(m * n)
  
  if(UNB==0){# balanced 
    tt = array( seq(1/m,1,by=1/m) %*% t(rep(1,n)), c(N,1))
  }else{ # unbalanced
    # draw t without replacement
    tt =sample(seq(.01,1, by=.01),m,replace = FALSE)
    for (i in 2:n) tt =c(tt,sample(seq(.01,1, by=.01),m,replace = FALSE))
  }
  ind_ind = array(rep(1,m)%*% t(1:n),c(N,1))
  
  dat=matrix(0,N,p)	
  
  DAG =genDAG(p,2/(p-1)*.7*rate,HUB)
  true_DAG=DAG$DAG
  
  NumOfCom = 10; 
  # kernel_type =2; kernel_par =7; loc_type =1;var_type=1;
  
  if (model==1){ # linear model
    f_t1 = GenRanFun(n, tt,ind_ind,NumOfCom,kernel_type,kernel_par,loc_type,var_type)
    
    dat[,1] = f_t1
    
    for (i in 2:p){
      f_ti=  as.matrix(dat[,1:(i-1)]) %*% true_DAG[1:(i-1),i] + GenRanFun(n, tt,ind_ind,NumOfCom,kernel_type,kernel_par,loc_type,var_type)
      
      dat[,i] = f_ti 
    }
  }
  return(list(dat = dat, tt =tt, ind_ind = ind_ind,true_DAG = true_DAG))
}


################################################################
# Eigenfunctions of Brownian motion kernel based on Mercer's theorem
# tt: 1x1 time; j: the j-th eigenfunction
# Kuang-Yao Lee
# 01-13-20
################################################################
BrownEigen <- function(tt, j){sqrt(2) * sin(pi * (j-1/2)%*% t(tt) )}

