################
# task: write the columns of C matrix into nifti files
require(R.matlab)
require(igraph)
source('../../Functions/canon_cor.R')
motor_kfs = readMat("../../../Data/Results/RealData/kki_motor_kfs_11.mat");
motor_pkfs = readMat("../../../Data/Results/RealData/kki_motor_pkfs_11.mat");
kfsA = list();
kfsC = list();
pkfsA = list();
pkfsC = list();
ndata = length(motor_kfs$kfsA);
for(i in 1:ndata){
	kfsA[[i]] = motor_kfs$kfsA[[i]]
	kfsC[[i]] = motor_kfs$kfsC[[i]]
	pkfsA[[i]] = motor_pkfs$pkfsA[[i]];
	pkfsC[[i]] = motor_pkfs$pkfsC[[i]];
}
require(AnalyzeFMRI)
mask3D <- f.read.nifti.volume("../../../Data/Raw/wRL_pre.nii")[,,,1]
mask = which(mask3D!=0)
d = dim(kfsA[[1]])[1] 
for(i in 1:d){
	img <- rep(0, prod(dim(mask3D)))
	imgp <- rep(0, prod(dim(mask3D)))
	img[mask] <-kfsC[[1]][,i]
	imgp[mask] <-pkfsC[[1]][,i]
	dim(img) <- dim(mask3D)
	dim(imgp) <- dim(mask3D)
	f.write.nifti(img, paste("../../../Data/Results/RealData/nifti/kfs-col-",i,sep=""), size="float",nii=TRUE)
	f.write.nifti(imgp, paste("../../../Data/Results/RealData/nifti/pkfs-col-",i,sep=""), size="float",nii=TRUE)
}

