% compare LDS and PLDS with zero penalty
% they should yield the same results

% add path of core functions
originpath = path;
newpath = genpath('../../Functions');
path(originpath,newpath);

load('../../../Data/Results/Simulation/p-300-d-10-T-100-data-regA-sparseC.mat');

[aa,cc,qq,rr,pipi,vv,Sx] = kfs_learn(y,a,c,q,r,Pi,v,tol,miter);

lambdaA = 0;
lambdaC = 0;
[aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);

difference = abs(aa - aap);
sum(difference(:))

difference = abs(cc - ccp);
sum(difference(:))

lambdaA = 10e-9;
lambdaC = 0.001;
[aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);

path(originpath);
