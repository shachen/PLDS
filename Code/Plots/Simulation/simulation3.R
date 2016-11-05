# estimate & prediction accuracy
# estimate accuracy
require(R.matlab)
source('../../Functions/canon_cor.R')
cancor = NULL

data = readMat('../../../Data/Results/Simulation/p-300-d-10-T-100-sim1-result.mat')
p = data$p
d = data$d
T = data$T
penaltyCe = as.vector(data$penaltyC)
penaltyAe = as.vector(data$penaltyA)
cormat = data$cormat # first half of cormat are correlation <aap,A>; second half <ccp,C>
npenal = length(penaltyCe)

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

# prediction accuracy
data = readMat('../../../Data/Results/Simulation/p-300-d-10-T-100-sim3-result.mat')
T_future = as.numeric(data$T.future)
penaltyCp = as.vector(data$penaltyC)
penaltyAp = as.vector(data$penaltyA)
pred_accuracy = as.vector(data$accuracy)

pdf("../../../Figures/pdf/p-300-d-10-T-100-sim3-accuracy.pdf")
# plot both accuracies in one graph
par(mar=c(5,4,4,4)+0.1)
# plot estimation accuracy
plot(log(penaltyCe[-1],4),smooth.spline(cancor[-1,4], df = 7)$y,axes=F,ylim=c(0,0.4),xlab="",ylab="",
     type="l",lwd=3,col="orange",main="Prediction/Estimation Accuracy vs Penalty", lty = 1)
axis(2,ylim=c(0,0.4),col="black")
mtext("Estimation Accuracy",side=2,line=2.5,adj=1)

# allow a second plot on same graph
par(new=T)

# plot pred accuracy
plot(log(penaltyCp[-1],4),smooth.spline(pred_accuracy[-1], df = 7)$y,axes=F,type="l",ylab="",
     xlab="",col="blue",ylim=c(4,9),lwd=3, lty = 2)
axis(4,ylim=c(4,10),col="black")
mtext("Prediction Accuracy",side=4,line=2.5,adj=1)

# plot x axis
axis(1,xlim=c(-11,2),col="black")
box(bty="o")
mtext(expression(paste("log10(",lambda[C],")")),side=1,line=2.5)

# legend
legend("topright",legend=c("Estimation","Prediction"),lwd=3,col=c("orange","blue"),bty="n", lty = c(1,2))
dev.off()

