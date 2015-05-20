# this code picks the optimal number of latents states

source("../../Functions/getElbows.R")
# step 1: calculate eigen values
files = list.files(path="../../../Data/Base",full.names=TRUE,pattern="KKI*")
subj = 1
data = read.csv(files[subj],header=TRUE)


cdata = scale(data,scale=F)
covmat = cov(t(cdata))
eigvals = eigen(covmat)$values
elbows = getElbows(eigvals,n=5)
var_pct = cumsum(eigvals)/sum(eigvals)
var_pct[elbows]
elbows


pdf("../../../Figures/pdf/EigenValues.pdf")
plot(cumsum(eigvals)/sum(eigvals),type="l",lwd=3,ylab="percentage",xlab="N",main="Variances Interpreted by Eigenvals",xlim=c(1,210))
grid(nx = 20, ny = 10, col = "lightgray", lty = "dotted")
dev.off()

t = dim(data)[1]
proflik = rep(0,t-1)
for(i in 4:(t-1)){
  mu1 = mean(eigvals[4:i])
  mu2 = mean(eigvals[(i+1):t])
  sigma = sqrt(((i-4)*var(eigvals[4:i])+(t-i-1)*var(eigvals[(i+1):t]))/(t-5))
  proflik[i] = sum(log(dnorm(eigvals[4:i],mean = mu1,sd = sigma))) + sum(log(dnorm(eigvals[(i+1):t],mean=mu2,sd=sigma)))
}
pdf("../../../Figures/pdf/profileLikelihood.pdf")
plot(proflik[-(1:3)],type="l",main="Profile Likelihood Plot",xlab="N",ylab="profile likelihood",lwd=3)
segments(11-3,min(proflik[!is.na(proflik)])-5,11-3,max(proflik[!is.na(proflik)]),lwd=3,col="Blue")
dev.off()
