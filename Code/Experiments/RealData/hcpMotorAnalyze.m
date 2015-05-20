function hcpMotorAnalyze()
%%  apply PLDS algorithm on HCP motor data and save results

    %   analyze first 1000 time points of real data with KFS and penalized KFS model
    %   HCP Data Motor Region

    %   import data
    %   use fsl smoothed data
    datapath = '../../../Data/Base/';
    yy(:,:)=csvread([datapath,'hcpmotordata_smoothed.csv'],1,0);
    
    %   relavant dimensions
    [T,p]=size(yy);
    m = 149;

    %   analyze
    y = squeeze(yy(:,:))';
    %[U,S,V] = svd(y,'econ');
    %a = eye(m);
    %c = U(:,1:m) * sqrt(S(1:m,1:m));
    q=eye(m);
    r=spdiags(ones(p,1),0,p,p);
    Pi=zeros(m,1);
    v=eye(m)*10e-3;
    tol = 10e-3;
    miter = 20;

    load('../../../Data/Results/RealData/hcp_init_val.mat');


    lambdaA = 0.000001;
    lambdaC = 0.00001;

    [aa,cc,qq,rr,pipi,vv,Sx]=kfs_learn(y,a,c,q,r,Pi,v,tol,miter);
    [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
    kfsA=aa;
    kfsC=cc;
    pkfsA=aap;
    pkfsC=ccp;

    save('../../../Data/Results/RealData/hcp_motor_partial_kfs_149.mat','kfsA','kfsC','Sx','qq','rr')
    save('../../../Data/Results/RealData/hcp_motor_partial_pkfs_149.mat','pkfsA','pkfsC','Sxp','qqp','rrp')
end
