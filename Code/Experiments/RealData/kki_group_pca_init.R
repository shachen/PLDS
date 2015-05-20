# run group pca on kirby 21 motor data to get initialization
# pca
files = list.files(path="../../../Data/Base",full.names=TRUE,pattern="KKI*")
data11 = read.csv(files[1],header=TRUE)
data12 = read.csv(files[4],header=TRUE)
data21 = read.csv(files[2],header=TRUE)
data22 = read.csv(files[3],header=TRUE)

require(fastICA)
m = 11
f11=fastICA(t(data11),n.comp=m)
f12=fastICA(t(data12),n.comp=m)
f21=fastICA(t(data21),n.comp=m)
f22=fastICA(t(data22),n.comp=m)

X11 = (f11$X)%*%(f11$K)
X12 = (f12$X)%*%(f12$K)
X21 = (f21$X)%*%(f21$K)
X22 = (f22$X)%*%(f22$K)

X = cbind(X11,X12,X21,X22)
XtX = t(X)%*%X

sv=svd(XtX)
Sigma=diag(sqrt(sv$d))
SigmaInv=diag(1/sqrt(sv$d))
U=sv$u
V=X%*%U%*%SigmaInv

V.l=V[,1:m]
Sigma.l=Sigma[1:m,1:m]
U.l=t(U)[1:m,]
X.app=V.l%*%Sigma.l%*%U.l

W0=c()
for(i in 1:4){
	W0=c(W0,list(t(solve(Sigma[1:m,1:m]%*%U.l[,((i-1)*m+1):(i*m)]))))
}

W11 = W0[[1]]
W12 = W0[[2]]
W21 = W0[[3]]
W22 = W0[[4]]

require(R.matlab)
writeMat(con="../../../Data/Results/RealData/kki_motor_init.mat",a11=W11,a12=W12,a21=W21,a22=W22,c11=X11,c12=X12,c21=X21,c22=X22)
