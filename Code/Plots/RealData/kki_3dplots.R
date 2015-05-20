library(fmri)
library(AnalyzeFMRI)
library(rgl)
library(misc3d)
mask <- f.read.nifti.volume("../../../Data/Raw/MNI152_T1_2mm_brain_mask.nii")[,,,1]
motormask = f.read.nifti.volume("../../../Data/Raw/wRL_pre.nii")[,,,1]
idx = which(motormask!=0)

N = 11 # for 11 latent states
tq = 0.9
colors = rainbow(N)
contour3d(mask, level= 1, smooth = 20, fill = TRUE, mask = array(TRUE, dim(mask)), alpha = .1, add = TRUE)
contour3d(motormask, level= 1, smooth = 20, fill = TRUE, mask = array(TRUE, dim(mask)), alpha = .25, add = TRUE)
for(i in 1:N)
{
	if(!(i %in% c(7,8,10,11))){
		img <- f.read.nifti.volume(paste("../../../Data/Results/RealData/nifti/kfs-col-",i,".nii",sep=''))[,,,1]
		thresh=quantile(abs(img[idx]),tq)
		timg <- (abs(img) > thresh) + 0

		contour3d(timg, level= 1, smooth = 20, fill = TRUE, mask = array(TRUE, dim(timg)), alpha = 1, add = TRUE, color = colors[i])
	}
}
# do a few snapshots
rgl.viewpoint(theta=180,phi=90,zoom=0.5)
rgl.snapshot("../../../Figures/png/view1.png")
rgl.viewpoint(theta=0,phi=0,zoom=0.6)
rgl.snapshot("../../../Figures/png/view2.png")
rgl.viewpoint(theta=0,phi=-90,zoom=0.6)
rgl.snapshot("../../../Figures/png/view3.png")
print("You Can Also Manually Rotate The 3D Plot And Get Some Other Snapshots!")
rgl.quit()




