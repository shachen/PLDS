# load data
require(R.matlab)
source("../../../Code/Functions/amari.R")
data = readMat("../../../Data/Results/RealData/kki_motor_kfs_11.mat")
pdata = readMat("../../../Data/Results/RealData/kki_motor_pkfs_11.mat")

# extract A and C matrices
kfsA = data$kfsA;
kfsC = data$kfsC;
a11 = kfsA[[1]]
a12 = kfsA[[2]]
a21 = kfsA[[3]]
a22 = kfsA[[4]]

pkfsA = pdata$pkfsA;
pkfsC = pdata$pkfsC;
a11 = pkfsA[[1]]
a12 = pkfsA[[2]]
a21 = pkfsA[[3]]
a22 = pkfsA[[4]]



# calculate amari errors
m = 11;
amari(a11,a12)
amari(a11,a21)
amari(a11,a22)
amari(a12,a22)
amari(a12,a21)
amari(a22,a21)

# calculate canonical correlations
source("../../../Code/Functions/canon_cor.R")
c1112 = cor(a11,a12)
c1121 = cor(a11,a21)
c1122 = cor(a11,a22)
c1221 = cor(a12,a21)
c1222 = cor(a12,a22)
c2122 = cor(a21,a22)
I = diag(rep(1,m))
sum(diag(pMatrix.min(c1112,I)$A))
sum(diag(pMatrix.min(c1121,I)$A))
sum(diag(pMatrix.min(c1122,I)$A))
sum(diag(pMatrix.min(c1221,I)$A))
sum(diag(pMatrix.min(c1222,I)$A))
sum(diag(pMatrix.min(c2122,I)$A))

# plot mari errors & canonical correlations
require(fields)
amari_mat = matrix(c(0,0.88,1.05,1.02,0.88,0,1.08,1.09,1.05,1.08,0,0.98,1.02,1.09,0.98,0),ncol=4,nrow=4)
cc = -log(c(c1112,c1121,c1122,c1221,c1222,c2122)/m)
cc_mat = matrix(c(0,0.076,0.105,0.095,0.076,0,0.095,0.095,0.105,0.095,0,0.085,0.095,0.095,0.085),ncol = 4, nrow=4)

for(i in 1:4){
	for(j in 1:4){
		if(i < j){
			cc_mat[i,j]=0.105
			amari_mat[i,j] = 1.1
		}
	}
}

cc_mat = t(cc_mat[4:1,])
amari_mat = t(amari_mat[4:1,])

pdf("../../../Figures/pdf/A-matrices-similarity.pdf")
par(mar=c(5,4.5,4,7))
image(cc_mat,xaxt="n",yaxt="n",col=heat.colors(64),main="Similarities Among Estimated A Matrices")
axis(1,at=c(0,1/3,2/3,1),labels=c(expression(A[11]),expression(A[12]),expression(A[21]),expression(A[22])))
axis(2,at=c(0,1/3,2/3,1),labels=c(expression(A[22]),expression(A[21]),expression(A[12]),expression(A[11])))
image.plot(cc_mat,legend.only=T,col=heat.colors(64))
dev.off()

par(mar=c(5,4.5,4,7))
image(amari_mat,xaxt="n",yaxt="n",col=heat.colors(64),main="Similarities Among Estimated A Matrices")
axis(1,at=c(0,1/3,2/3,1),labels=c(expression(A[11]),expression(A[12]),expression(A[21]),expression(A[22])))
axis(2,at=c(0,1/3,2/3,1),labels=c(expression(A[22]),expression(A[21]),expression(A[12]),expression(A[11])))
image.plot(amari_mat,legend.only=T,col=heat.colors(64))

# connectivity graph
require(igraph)
m = 11
motor_pkfs = readMat("../../../Data/Results/RealData/kki_motor_pkfs_11.mat");
pkfsA = motor_pkfs[[1]]
thresh = 0.8
index = which(abs(pkfsA[[1]])>=quantile(abs(pkfsA[[1]]),thresh))
index = setdiff(index, (0:10)*12+1)
intensity = pkfsA[[1]][index]
cols = rep("blue",length(intensity))
cols[which(intensity<=0)] = "red"
graphmat = (abs(pkfsA[[1]])>=quantile(abs(pkfsA[[1]]),thresh))
diag(graphmat) = 0
graph = graph.adjacency(graphmat)
E(graph)$color = cols
pdf("../../../Figures/pdf/ConnectivityGraph_11.pdf")
plot(graph,main="Connectivity Graph",edge.width=4*abs(intensity),edge.arrow.size=0.7,vertex.color="green", vertex.frame.color="darkgreen")
dev.off()


# 4 connectivity graph in one plot
pdf("../../../Figures/pdf/ConnectivityGraph_11_all.pdf")
par(mfrow=c(2,2))
for(i in 1:4){
	thresh = 0.8
	index = which(abs(pkfsA[[i]])>=quantile(abs(pkfsA[[i]]),thresh))
	index = setdiff(index, (0:(m-1))*(m+1)+1)
	intensity = pkfsA[[i]][index]
	cols = rep("blue",length(intensity))
	cols[which(intensity<=0)] = "red"
	graphmat = (abs(pkfsA[[i]])>=quantile(abs(pkfsA[[i]]),thresh))
	diag(graphmat) = 0
	graph = graph.adjacency(graphmat)
	E(graph)$color = cols
	plot(graph,main=paste("Subject ",floor((i-0.01)/2)+1,", ","Scan", i%%2 + 1),edge.width=4*abs(intensity),edge.arrow.size=0.7,vertex.color="green", vertex.frame.color="darkgreen")
}
dev.off()

# interactive graph plots

i = 1
thresh = 0.8
index = which(abs(pkfsA[[i]])>=quantile(abs(pkfsA[[i]]),thresh))
index = setdiff(index, (0:(m-1))*(m+1)+1)
intensity = pkfsA[[i]][index]
cols = rep("blue",length(intensity))
cols[which(intensity<=0)] = "red"
graphmat = (abs(pkfsA[[i]])>=quantile(abs(pkfsA[[i]]),thresh))
diag(graphmat) = 0
graph = graph.adjacency(graphmat)
E(graph)$color = cols
tkplot(graph,main=paste("Subject ",floor((i-0.01)/2)+1,", ","Scan", i%%2 + 1),edge.width=4*abs(intensity),edge.arrow.size=0.7,vertex.color="green", vertex.frame.color="darkgreen")

