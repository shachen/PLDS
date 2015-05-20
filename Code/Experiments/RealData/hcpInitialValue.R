# calculate initial values for full hcpdata
data = read.csv("../../../Data/Base/hcpmotordata_smoothed.csv",header=T)
data = scale(data,scale=F)
require(fastICA)
d = 149 # number of latents states picked with profile-likelihood method
fica = fastICA(t(data),n.comp=d) # use fastICA for dimension reduction
c = (fica$X) %*% (fica$K) # project to first d principal directions

require(vars)
y = fica$K
colnames(y) = paste("y",1:d,sep='')
var = VAR(y,p=1,type="none") # use vector AR model to estimate A
a = NULL
for(i in 1:d){
	a = rbind(a,var$varresult[[i]]$coefficients)
}

require(R.matlab)
writeMat(con="../../../Data/Results/RealData/hcp_init_val.mat",a=a,c=c)

# calculate initial values with partial hcpdata
N = 1000 # the number of data points to use

dataAll = read.csv("../../../Data/Base/hcpmotordata_smoothed.csv",header=T)
data = scale(dataAll,scale=F)[1:N,]
require(fastICA)
d = 149 # number of latents states picked with profile-likelihood method
fica = fastICA(t(data),n.comp=d) # use fastICA for dimension reduction
c = (fica$X) %*% (fica$K) # project to first d principal directions

require(vars)
y = fica$K
colnames(y) = paste("y",1:d,sep='')
var = VAR(y,p=1,type="none") # use vector AR model to estimate A
a = NULL
for(i in 1:d){
	a = rbind(a,var$varresult[[i]]$coefficients)
}

require(R.matlab)
writeMat(con="../../../Data/Results/RealData/hcp_partial_init_val.mat",a=a,c=c,pcax=y[N,])

