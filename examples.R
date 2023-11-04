###############################################################################
###############################################################################
###############################################################################
# Simulations

# primary parameters
n=50 # sample size
p=50 # dimensionality, i.e. the number of random functions
m=50 # number of obserational time points

# other parameters
rate = 1.5 #sparsity, "rate =1 (1.5)" corresponds to "q=0.7 (1.05)"
sig=0.05
model=1
HUB=0
UNB=0
NumofBasis=10
NumOfDim=2
eps_1=.1;
eps_2=.5;
STD=0
kernel_T=2;
kernel_type=2;kernel_par=7;loc_type=1;
var_type<- if(UNB==1) 2 else 5
eps_T<- if(UNB==1) .01 else .01

# generate multivariate functional data
X=GenMod(n,p,m,model,UNB,STD,HUB=0,rate,kernel_type,kernel_par,loc_type,var_type) #old model

# fucntional PC-algorithm
Q = qmat(NumofBasis)
Qn = qmat(n)

A = matrix(0, n,p * NumOfDim)
for (i in 1:p){
  XX<-CalCoor2(X$dat[,i], X$ind_ind,X$tt,kernel_T, NumofBasis,eps_T,UNB,K_T_SCORE)
  
  XX = XX$XX
  # K-L expansion
  eig_K=eigen(Q %*% t(XX) %*% Qn %*% XX %*% Q);d = eig_K$values;U=eig_K$vectors
  
  # calculate fPC scores
  A[,(1+(i-1)*NumOfDim):(i*NumOfDim)] = ( Qn %*% XX %*% Q ) %*% U[,1:NumOfDim]
}

suffStat <- list(GX_orth = A, epsilon_1 = eps_1, epsilon_2 = eps_2, NumOfDim = NumOfDim)
pc.fit <- pc(suffStat=suffStat, indepTest=CalfPCO, p = p,alpha = sig) #ridge-type regularization


# compute the structural Hamming distance (SHD) between the estiamted
# graph and the truth
shd_am(wgtMatrix(pc.fit),dag2cpdag(X$true_DAG))

# compute the true discovery rate (TDR) between the estiamted
# graph and the truth
compareGraphs_am(wgtMatrix(pc.fit),dag2cpdag(X$true_DAG))[3]


###############################################################################
###############################################################################
###############################################################################
# Time-course data analysis
# load in data
dat <- read.csv(file="insilico.csv", header=TRUE, sep=",")
trueG <- read.csv(file="trueGraph.csv", header=FALSE, sep=",")
NODE <- read.csv(file="nodeNames.csv", header=FALSE, sep=",")
diag(trueG)=0
trueG=as.matrix(trueG)

# data pre-processing
IND = matrix(0,660,20)
IND[,1]= (dat$Inhibitor=="" & dat$Stimulus=="loLIG1")
IND[,2]= (dat$Inhibitor=="INH1" & dat$Stimulus=="loLIG1")
IND[,3]= (dat$Inhibitor=="INH2" & dat$Stimulus=="loLIG1")
IND[,4]= (dat$Inhibitor=="INH3" & dat$Stimulus=="loLIG1")
IND[,5]= (dat$Inhibitor=="" & dat$Stimulus=="hiLIG1")
IND[,6]= (dat$Inhibitor=="INH1" & dat$Stimulus=="hiLIG1")
IND[,7]= (dat$Inhibitor=="INH2" & dat$Stimulus=="hiLIG1")
IND[,8]= (dat$Inhibitor=="INH3" & dat$Stimulus=="hiLIG1")

IND[,9]= (dat$Inhibitor=="" & dat$Stimulus=="loLIG2")
IND[,10]= (dat$Inhibitor=="INH1" & dat$Stimulus=="loLIG2")
IND[,11]= (dat$Inhibitor=="INH2" & dat$Stimulus=="loLIG2")
IND[,12]= (dat$Inhibitor=="INH3" & dat$Stimulus=="loLIG2")
IND[,13]= (dat$Inhibitor=="" & dat$Stimulus=="hiLIG2")
IND[,14]= (dat$Inhibitor=="INH1" & dat$Stimulus=="hiLIG2")
IND[,15]= (dat$Inhibitor=="INH2" & dat$Stimulus=="hiLIG2")
IND[,16]= (dat$Inhibitor=="INH3" & dat$Stimulus=="hiLIG2")

IND[,17]= (dat$Inhibitor=="" & dat$Stimulus=="loLIG1+loLIG2")
IND[,18]= (dat$Inhibitor=="" & dat$Stimulus=="loLIG1+hiLIG2")
IND[,19]= (dat$Inhibitor=="" & dat$Stimulus=="hiLIG1+loLIG2")
IND[,20]= (dat$Inhibitor=="" & dat$Stimulus=="hiLIG1+hiLIG2")

#remove IND5,7,19,20 
dat<-dat[((IND[,5]==1) + (IND[,7]==1) + (IND[,19]==1) + (IND[,20]==1))==0,]

m=11;n <- dim(dat)[1]/m;p=20
ind_ind=rep(1:n,11)
tt = c(0, 1,2,4,6,10,15,30,45,60,120)
ttt = array(t(rep(1,n) %*% t(tt)),c(n*m,1))

X<-matrix(0,n*m,p)
for (i in 1:n){
  X[((1+(i-1)*m):(i*m)),] = as.matrix(dat[(ind_ind==i),-c(1:4)]) }

# correct ind_ind indexing
ind_ind=array(rep(1,m) %*% t(1:n),c(n*m,1))


# parameters 
kernel_T=2
UNB=0;
NumOfDim=2 
eps_T=.01
NumofBasis = 10
Q = qmat(NumofBasis)
Q_n = qmat(n)


# functional PC-algorithm
A = matrix(0, n,p * NumOfDim)
for (i in 1:p){
  XX<-CalCoor2(X[,i], ind_ind,ttt,kernel_T=2, NumofBasis,eps_T,UNB,K_T_SCORE)
  
  XX = XX$XX
  # K-L expansion
  eig_K=eigen(Q %*% t(XX) %*% Q_n %*% XX %*% Q);d = eig_K$values;U=eig_K$vectors
  
  # calculate fPC scores
  A[,(1+(i-1)*NumOfDim):(i*NumOfDim)] = (Q_n %*% XX %*% Q ) %*% U[,1:NumOfDim]
}

eps_1=.1;
eps_2=.5;
sig=.05

suffStat <- list(GX_orth = A, epsilon_1 = eps_1, epsilon_2 = eps_2, NumOfDim = NumOfDim)
pc.fit <- pc(suffStat, indepTest=CalfPCO, p = p,alpha = sig) #ridge-type regularization

# compute the structural Hamming distance (SHD) between the estiamted
# graph and the truth
shd_am(wgtMatrix(pc.fit),dag2cpdag(trueG))

# compute the true discovery rate (TDR) between the estiamted
# graph and the truth
compareGraphs_am(wgtMatrix(pc.fit),dag2cpdag(trueG))[3]

# generate the estimated graph
g_fpc=graph.adjacency(wgtMatrix(pc.fit),mode="directed")
coords <- layout_in_circle(g_fpc)
par(mar=c(.1,.1,.1,.1))
plot(g_fpc,vertex.label=NODE$V1, vertex.shape="none",
     vertex.color="none",
     vertex.label.dist=0,
     edge.color="black",vertex.label.cex = 1.2,layout=coords,edge.width=.3,
     vertex.label.color="blue", edge.arrow.size=0.9, edge.arrow.width=1,frame=F)



