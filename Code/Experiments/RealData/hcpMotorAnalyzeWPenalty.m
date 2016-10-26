function hcpMotorAnalyzeWPenalty(lambdaA, lambdaC)
%%  apply PLDS algorithm on HCP motor data and save results

    %   analyze first 1000 time points of real data with KFS and penalized KFS model
    %   HCP Data Motor Region
    %   lambdaA, lambdaC: penalty terms

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


%     lambdaA = 0.000001;
%     lambdaC = 0.00001;

    [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
    pkfsA=aap;
    pkfsC=ccp;

    save(strcat('../../../Data/Results/RealData/hcp_motor_partial_pkfs_149_', num2str(lambdaA), '_', num2str(lambdaC), '.mat'),'pkfsA','pkfsC','Sxp','qqp','rrp')
end

% run the following to get analysis result
% values = [0.00000000001, 0.0000000001, 0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001];
% for i = 1:7
%      hcpMotorAnalyzeWPenalty(values(i), values(i))
% end