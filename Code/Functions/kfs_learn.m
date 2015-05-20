function [A,C,Q,R,Pi,V,Sx]=kfs_learn(y,a,c,q,r,pi,v,tol,miter)
%   reference: Zoubin etal, Parameter Estimation for LDS
%   EM algo to learn parameters in linear dynamical system(LDS)
%   author: Shaojie CHEN
%   refer to KFS code for model details
%
%   use as: [A,C,Q,R,Pi,V]=kfs_learn(y,a,c,q,r,pi,v,tol,miter)
%   input:
%   y: observed data, dimension p by T
%   a,c,q,r,pi,v: inital parameter values
%   miter: max number of iterations
%   tol: tolerance for difference between two successive iterations
%
%   output:
%   A,C,Q,R,Pi,V: final parameter estimations

if nargin < 9
    tol=1e-6;
end
if nargin < 8
    miter=100;
end

[~,d]=size(q);
[p,T]=size(y);
A=a;C=c;Q=q;R=r;Pi=pi;V=v;
iter=1;dif=1;
Pt=zeros(d,d,T);
Ptt_1=zeros(d,d,T-1);

while (iter<miter)&&(dif>tol)
    disp(iter);
    %   E step:
    [~,~,~,~,Sx,Sv,Scov]=KFS(A,C,Q,R,Pi,V,y);
    for t=1:T
        Pt(:,:,t)=Sv(:,:,t)+Sx(:,t)*Sx(:,t)';
    end
    for t=1:(T-1)
        Ptt_1(:,:,t)=Scov(:,:,t)+Sx(:,t+1)*Sx(:,t)';
    end

    %   M step:
    Cnew=(y*Sx')/(sum(Pt,3));
    vRnew = (sum(y.*y,2)-sum((Cnew*Sx).*y,2))/T;
    %	if we want R to be diagonal
    Rnew = spdiags(vRnew,0,p,p);
    %	if we want R to be mutiple of Identity
    %	Rnew = mean(vRnew) * spdiags(ones(p,1),0,p,p); 
    Anew=(sum(Ptt_1,3))/(sum(Pt(:,:,1:(T-1)),3));
    Pinew=Sx(:,1);
    
    %	constraints:
    Qnew=eye(d);
    Vnew = V;
    
    %   updates
    dif=norm(Cnew-C)+norm(diag(Rnew-R))+norm(Anew-A)+norm(Qnew-Q)+norm(Pinew-Pi)+norm(Vnew-V);
    iter=iter+1;    
    C=Cnew;R=Rnew;A=Anew;Q=Qnew;Pi=Pinew;
    V=Vnew;
end
end
