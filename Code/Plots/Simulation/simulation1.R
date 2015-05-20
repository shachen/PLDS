require(R.matlab)
source('../../Functions/canon_cor.R')
cancor = NULL

data = readMat('../../../Data/Results/Simulation/p-300-d-10-T-100-sim1-result.mat')
p = data$p
d = data$d
T = data$T
penaltyC = as.vector(data$penaltyC)
penaltyA = as.vector(data$penaltyA)
cormat = data$cormat # first half of cormat are correlation <aap,A>; second half <ccp,C>
npenal = length(penaltyC)

for(i in 1:npenal){
	aacor = cormat[[1]]
	cccor = cormat[[1 + npenal]]
	aapcor = cormat[[i]]
	ccpcor = cormat[[npenal + i]]

	if(sum(is.na(aapcor)) | sum(is.na(ccpcor))){
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

plotname = paste("../../../Figures/pdf/","p-",p,"-d-",d,"-T-",T,"-sim1-plot.pdf",sep='')
pdf(plotname)
matplot(cancor,type="l")
dev.off()

