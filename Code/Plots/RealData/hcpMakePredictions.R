############################### section 1: compare prediction results: pca and lds #
data = read.csv("../../../Data/Base/hcpmotordata_smoothed.csv",header=T)

require(R.matlab)
pca_est = readMat("../../../Data/Results/RealData/hcp_partial_init_val.mat")
pca_A = pca_est$a
pca_C = pca_est$c
pca_x = pca_est$pcax
pca_pred = NULL
for(i in 1:200){
	pca_x = pca_A %*% pca_x;
	pca_pred = cbind(pca_pred, pca_C %*% pca_x)
}

lds_est = readMat("../../../Data/Results/RealData/hcp_motor_partial_kfs_149.mat")
lds_A = lds_est$kfsA
lds_C = lds_est$kfsC
lds_x = lds_est$Sx[,1000]
lds_pred = NULL
for(i in 1:200){
	lds_x = lds_A %*% lds_x
	lds_pred = cbind(lds_pred, lds_C %*% lds_x)
}

pca_accy = colSums((scale(t(data[1001:1200,])) - scale(pca_pred))^2)
lds_accy = colSums((scale(t(data[1001:1200,])) - scale(lds_pred))^2)
pca_accy = pca_accy/dim(data)[2] * 100
lds_accy = lds_accy/dim(data)[2] * 100

pdf("../../../Figures/pdf/hcp_pred_accy.pdf")
par(mar=c(5, 5, 4, 5) + 0.1)
plot(lds_accy, type="l", lwd=3, main="",lty = 1, col ="green", axes=FALSE,xlab="",ylab="",ylim=c(0,0.4))
axis(2,ylim = c(0,0.4),col="green",las=1,at=0.05 *(0:8))
mtext("PLDS Prediction MSE",side=2,line = 3,col="green",adj=0.75,col.axis="green")
axis(1,xlim=c(1,200))
mtext("Time",side = 1, line = 2.5)
box()
par(new=TRUE)
plot(pca_accy, type="l", lwd=3,lty = 2, col ="red",axes=FALSE,xlab="",ylab="",ylim=c(0,200))
axis(4,ylim = c(0,200),col="red",las=1,col.axis="red")
mtext("PCA Prediction MSE", side=4, line = 3, col="red",adj=0.75)
legend("topright",legend=c("PCA","PLDS"),lty=2:1,col=c("red","green"),lwd=3)
dev.off()


# plot a sample of predicted time series, compare with truth & pca predictions
variances = diag(lds_est$rr)
sorted = sort(variances,index.return=T)$ix
idx = sorted[1:20]
	#sd = sqrt(sum(variances[idx])/(length(idx))^2)
sd = sqrt(variances[idx])
truth = scale(t(data[1001:1200,]))
lds_predn = scale(lds_pred)
pca_predn = scale(pca_pred)
scales = apply(data[1001:1200,],1,sd)

d = dim(lds_A)[1]
var_1 = list()
var_1[[1]] = diag(rep(1,d))
for(i in 2:200){
	var_1[[i]] = lds_A %*% var_1[[i-1]] %*% t(lds_A)
}

var_2 = var_1
for(i in 1:200){
	print(i)
	var_2[[i]] = (lds_C %*% var_2[[i]] %*% t(lds_C))[idx,idx]/(scales[i])^2
}

var_3 = var_2
for(i in 2:200){
	var_3[[i]] = var_3[[i]] + var_3[[i-1]]
}
for(i in 1:200){
	var_3[[i]] = var_3[[i]] + diag(sd^2)/(scales[i])^2
}

var_4 = rep(0, 200)
for(i in 1:200){
	var_4[i] = mean(as.vector(var_3[[i]]))
}

middle = colMeans(lds_predn[idx,])

upper = middle + 0.84 * sqrt(var_4) # 60% confidence interval
lower = middle - 0.84 * sqrt(var_4)

pca_middle = colMeans(pca_predn[idx,])
truth_middle = colMeans(truth[idx,])


pdf("../../../Figures/pdf/hcpSampleTS.pdf")
par(mar=c(5, 5, 4, 5) + 0.1)
combined = cbind(upper,middle,lower,truth_middle)
matplot(combined[1:30,],col= c("green","green","green","blue"),type="l",lty=c(2,1,2,1),main="",xlab="",ylab="",lwd=2,axes=FALSE,ylim=c(-3.36,-3.23))
axis(2,ylim=c(-3.36,-3.23))
mtext("PLDS Prediction/True Value",side = 2, line = 3,adj = 0.75)
axis(1,xlim=c(1,30))
mtext("Time",side = 1, line = 2.5)
box()
par(new=TRUE)
plot(pca_middle[1:30], type="l", lty = 1, col = "red",main="",xlab="",ylab="",lwd=2,axes =FALSE,ylim = c(-3.36,0.05))
axis(side = 4, ylim = c(-3.36,0.05))
mtext("PCA Prediction",side = 4, line = 3,adj = 0.75)
legend("topright",legend=c("PLDS Prediction","True Value","PCA Prediction"),col=c("green","blue","red"),lwd = 2)
dev.off()

