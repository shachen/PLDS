function kkiMotorAnalyze()
%   analyze real data with KFS and penalized KFS model
%   KIRBY 21 Motor Network Analysis

% add path of core functions
originpath = path;
newpath = genpath('../../Functions');
path(originpath,newpath);

%   import data
yy=[];
%   use fsl smoothed data
datapath = '../../../Data/Base/';
yy(1,:,:)=csvread([datapath,'KKI2009-01.csv'],1,0);
yy(2,:,:)=csvread([datapath,'KKI2009-25.csv'],1,0);
yy(3,:,:)=csvread([datapath,'KKI2009-04.csv'],1,0);
yy(4,:,:)=csvread([datapath,'KKI2009-11.csv'],1,0);

%   relavant dimensions
[ndata,T,p]=size(yy);
m=11;

%   results containers
kfsA={};
kfsC={};
pkfsA={};
pkfsC={};

load('../../../Data/Results/RealData/kki_motor_init.mat');

%   analyze
for i = 1:4
y = squeeze(yy(i,:,:))';
%[U,S,V] = svd(y,'econ');
%a = eye(m);
%c = U(:,1:m) * sqrt(S(1:m,1:m));
q=eye(m);
r=spdiags(ones(p,1),0,p,p);
Pi=zeros(m,1);
v=eye(m)*10e-3;
tol = 10e-3;
miter = 20;

lambdaA = 0.0000000001;
lambdaC = 0.00000001;

if i == 1
    a = a11;
    c = c11;
elseif i == 2
    a = a12;
    c = c12;
elseif i == 3
    a = a21;
    c = c21;
else
    a = a22;
    c = c22;
end
       

[aa,cc,qq,rr,pipi,vv,Sx]=kfs_learn(y,a,c,q,r,Pi,v,tol,miter);
[aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);

kfsA{i}=aa;
kfsC{i}=cc;
pkfsA{i}=aap;
pkfsC{i}=ccp;
end

save('../../../Data/Results/RealData/kki_motor_pkfs_11.mat','pkfsA','pkfsC')
save('../../../Data/Results/RealData/kki_motor_kfs_11.mat','kfsA','kfsC')
path(originpath);
end
