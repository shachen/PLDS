function [A,C,Q,R,Pi,V,Sx]=kfs_learn_p(y,a,c,q,r,pi,v,tol,miter,lambdaA,lambdaC)
%   this is the parameter learning function with penalty terms
%   this version might not work for very high dimensional data
%   reference: Zoubin etal, Parameter Estimation for LDS
%   EM algo to learn parameters in linear dynamical system(LDS)
%   author: Shaojie CHEN
%   refer to KFS code for model details

%   input:
%   y: observed data, dimension p by T
%   a,c,q,r,pi,v: inital parameter values
%   miter: max number of iterations
%   tol: tolerance for difference between two successive iterations
%   lambda1, lambda2 are tuning parameters for A and C respectively

%   output:
%   A,C,Q,R,Pi,V: final parameter estimations

if nargin < 9
    tol=1e-2;
end
if nargin < 8
    miter=20;
end

[~,d]=size(q);
[p,T]=size(y);
A=a;C=c;Q=q;R=r;Pi=pi;V=v;
iter=1;dif=1;
Pt=zeros(d,d,T);
Ptt_1=zeros(d,d,T-1);
fistiter = 25;

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
    %   estimate C
    l2regA = sum(Pt,3);
    l2regb=Sx*(sqrt(inv(R))*y)';
    Cvec = l2reg(l2regA,l2regb,lambdaC*sqrt(diag(R)));
    Cnew = sqrt(R)*(reshape(Cvec,d,p))';

    %   estimate R
    vRnew = (sum(y.*y,2)-sum((Cnew*Sx).*y,2))/T;
    Rnew = spdiags(vRnew,0,p,p);
    %   if we want R to be sigma^2 * I
    %   Rnew = mean(vRnew) * spdiags(ones(p,1),0,p,p);
    
    %   estimate A
    %   Anew=(sum(Ptt_1,3))/(sum(Pt(:,:,1:(T-1)),3));
    fistA = sum(Pt(:,:,1:(T-1)),3);
    fistb = reshape(sum(Ptt_1,3)',d*d,1);
    %   fistL: refer to Example 2.2 in FISTA paper
    eigs = eig(sum(Pt(:,:,1:(T-1)),3));
    fistL = 2*eigs(d,1);
    Avec=fista(fistA,d,fistb,lambdaA,fistL,fistiter,10e-5); 
    Anew=(reshape(Avec,d,d))';
                
    Pinew=Sx(:,1);
    %   V=Vnew never changes
    Vnew = V;
    
    %   constraints:
    Qnew=eye(d);
    
    %   updates
    dif=normest(Cnew-C)+norm(diag(Rnew-R))+norm(Anew-A)+norm(Qnew-Q)+norm(Pinew-Pi)+norm(Vnew-V);
    iter=iter+1;
    %   fprintf('Iteration %d: dif=%f.\n',iter,dif);
    %   fprintf('I%d:A=%f,C=%f,Q=%f,R=%f,Pi=%f,V=%f\n',iter,A,C,Q,R,Pi,V);
    C=Cnew;R=Rnew;A=Anew;Q=Qnew;Pi=Pinew;V=Vnew;
end

end
