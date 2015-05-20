% simulation 1

% add path of core functions
originpath = path;
newpath = genpath('../../Functions');
path(originpath,newpath);

% load data
load('../../../Data/Results/Simulation/p-300-d-10-T-100-data-regA-sparseC.mat');

penaltyC = linspace(-11,2,14);
penaltyC = [0 exp(log(4)*penaltyC)];
penaltyA = linspace(-30,-17,14);
penaltyA = [0 exp(log(4)*penaltyA)];
[~,npenal] = size(penaltyC);

cormat={};

for i = 1:npenal
    lambdaA = penaltyA(i);
    lambdaC = penaltyC(i);
    
    [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
    cormat{i,1}=corr(A,aap);
    cormat{i,2}=corr(C,ccp);
end

save(['../../../Data/Results/Simulation/p-',num2str(p),'-d-',num2str(d),'-T-',num2str(T),'-sim1-result.mat'],'p','d','T','penaltyA','penaltyC','cormat');

% next step: use the simulation1.R for further analysis and plots
