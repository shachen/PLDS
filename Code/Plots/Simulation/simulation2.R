require(R.matlab)
source('../../Functions/canon_cor.R')
cancor = NULL

data = readMat('../../../Data/Results/Simulation/p-300-d-10-T-100-sim2-result.mat')
p = as.numeric(data$p)
d = as.numeric(data$d)
T = as.numeric(data$T)
penaltyC = as.vector(data$penaltyC)
penaltyA = as.vector(data$penaltyA)
cormat = data$cormat # first half of cormat are correlation <aap,A>; second half <ccp,C>
npenal = length(penaltyC)

A = data$A;
C = data$C;
aa = data$aa;
cc = data$cc;
paa = data$paa;
pcc = data$pcc;
paacor = cor(paa,A);
pcccor = cor(pcc,C);
I=diag(rep(1,dim(paacor)[1]))
tmp = pMatrix.min(paacor,I)
aap = paa[,tmp$pvec]
tmp = pMatrix.min(pcccor,I)
ccp = pcc[,tmp$pvec]
C = scale(C)
cc = scale(cc)
ccp = scale(ccp)
A = scale(A)
aa = scale(aa)
aap = scale(aap)
C = apply(C,2,rev)

for(i in 1:d){
	avg1 = mean(ccp[1:(p/2),i]);
	avg2 = mean(ccp[(p/2+1):p,i]);
	if(avg1 < avg2){
		ccp[,i] = -ccp[,i];
	}
}

for(i in 1:d){
	avg1 = mean(cc[1:(p/2),i]);
	avg2 = mean(cc[(p/2+1):p,i]);
	if(avg1 < avg2){
		cc[,i] = -cc[,i];
	}
}

for(i in 1:npenal){
	aacor = cormat[[1]]
	cccor = cormat[[1 + npenal]]
	aapcor = cormat[[i]]
	ccpcor = cormat[[npenal + i]]

	if(sum(is.na(aapcor)) | sum(is.na(ccpcor))){
		cat(paste("penalty ", i," bad result!\n"));
		next;
	}

	I=diag(rep(1,dim(aacor)[1]))
	xaa = mean(diag(pMatrix.min(aacor,I)$A))
	xaap = mean(diag(pMatrix.min(aapcor,I)$A))
	xcc = mean(diag(pMatrix.min(cccor,I)$A))
	xccp = mean(diag(pMatrix.min(ccpcor,I)$A))
	newrow = c(xaa,xaap,xcc,xccp)
	cancor = rbind(cancor,newrow)
}

plotname = paste("../../../Figures/pdf/","p-",p,"-d-",d,"-T-",T,"-sim2-plot.pdf",sep='')
pdf(plotname)
matplot(cancor,type="l")
dev.off()

require(fields)
require(colorRamps)
require(RColorBrewer)
pdf(paste("../../../Figures/pdf/","p-",p,"-d-",d,"-T-",T,"-sim2-heatmap.pdf",sep=''),width=10,height=7)
par(mfrow=c(2,3))
col=tim.colors(n = 64);
#col = rainbow(n=12)
#col = matlab.like(n=12)
#col = matlab.like2(n=64)
colA = col
colC = col
image.plot(C,main="C",axes=FALSE,col=colC,cex = 2)
image.plot(cc,main= expression(hat(C)[-infinity]),axes=FALSE,col=colC,cex = 2)
image.plot(ccp,main=expression(hat(C)[lambda[m]]),axes=FALSE,col=colC,cex = 2)
image.plot(t(A[d:1,]),main="A",axes=FALSE,col=colA,cex = 2)
image.plot(t(aa[d:1,]),main=expression(hat(A)[-infinity]),axes=FALSE,col=colA,cex = 2)
image.plot(t(aap[d:1,]),main=expression(hat(A)[lambda[m]]),axes=FALSE,col=colA,cex = 2)
dev.off()
