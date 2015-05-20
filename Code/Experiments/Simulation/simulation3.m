% prediction accuracy study

% set seed for reproducibility
rng(100);

% calculate future data
load('../../../Data/Results/Simulation/p-300-d-10-T-100-data-regA-sparseC.mat');
T_future = 10;
future_x = zeros(m,T_future);
future_y = zeros(n,T_future);

tmp_x = x(:,T);
tmp_y = y(:,T);

for t = 1:T_future
    ww = mvnrnd(zeros(m,1),Q);
    tmp_x = A * tmp_x + ww';
    future_x(:,t) = tmp_x;
    
    vtmp = normrnd(0,1,n,1);
    vv = sqrtRdiag .* vtmp;   
    future_y(:,t) = C * tmp_x + vv;
end

% make estimations & predictions
penaltyC = linspace(-11,2,14);
penaltyC = [0 exp(log(4)*penaltyC)];
penaltyC = penaltyC([1:2 5:13 15]);
penaltyA = linspace(-30,-17,14);
penaltyA = [0 exp(log(4)*penaltyA)];
penaltyA = penaltyA([1:2 5:13 15]);
[~,npenal] = size(penaltyC);

accuracy = zeros(1,npenal);
for i = 1:npenal
    lambdaA = penaltyA(i);
    lambdaC = penaltyC(i);
    [aap,ccp,qqp,rrp,pipip,vvp,Sxp]=kfs_learn_p(y,a,c,q,r,Pi,v,tol,miter,lambdaA,lambdaC);
    x_next = x(:,T);
    for j = 1:T_future
        x_next = aap * x_next;
        y_next = ccp * x_next;
        accuracy(1,i) = accuracy(1,i) + abs(corr(y_next,future_y(:,j)));
    end
end

plot(smooth(accuracy))
% save results
save(['../../../Data/Results/Simulation/p-',num2str(p),'-d-',num2str(d),'-T-',num2str(T),'-sim3-result.mat'],'p','d','T','penaltyA','penaltyC','accuracy','T_future');
