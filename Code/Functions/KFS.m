function [Fv1,Fv2,Fx1,Fx2,Sx,Sv,Scov]=KFS(A,C,Q,R,x0,v0,y)
%   kalman filter and smoother
%   author: Shaojie Chen
%   Version 1.2
%
%   The model:
%   x(t+1) = A x(t) + w(t), x(t) is d by 1
%   y(t) = C x(t) + v(t), y(t) is p by 1
%   w(t) ~ iid N(0,Q)
%   v(t) ~ iid N(0,R)
%   t = 1,2,...,T
%   x(0) ~ N(x0,v0)
%
%   Input:
%   A,C,Q,R,x0,v0 are as above
%   y is observation matrix, dimension p by T
%
%   Output:
%   Fx1,Fx2,Sx are d by T matrices
%   Fv1,Fv2,Sv,Scov are d by d by T arrays
%   Fv1, (:,:,t)  is Var(x(t)|y(1:t))
%   Fv2, (:,:,t)  is Var(x(t+1)|y(1:t))
%   Fx1, column t is E(x(t)|y(1:t))
%   Fx2, column t is E(x(t+1)|y(1:t))
%   Sx,  column t is E(x(t)|y(1:T))
%   Sv,  (:,:,t)  is Var(x(t)|y(1:T))
%   Scov,(:,:,t)  is Cov(x(t),x(t+1)|y(1:T))

%   initilizations
[~,T]=size(y);
[d,~]=size(A);
Fv1=zeros(d,d,T);
Fv2=zeros(d,d,T);
Fx1=zeros(d,T);
Fx2=zeros(d,T);
Sx=zeros(d,T);
Sv=zeros(d,d,T);
Scov=zeros(d,d,T);

%   CRC is used a lot, so save it separately to avoid multiple calculations
CRC = C' * inv(R) * C;
invR = inv(R);
%   Filter
for t=1:T
    if t==1
        %Fv1, Fx1 Calculated with Woodbury Matrix Identity
        Fv1(:,:,t)=v0-(v0*CRC*v0 - v0*CRC*inv(inv(v0)+CRC)*CRC*v0);
        Fv2(:,:,t)=A*Fv1(:,:,t)*A'+Q;
        Fx1(:,t)=x0+v0*C'*inv(R)*(y(:,t)-C*x0)-v0*CRC*inv(inv(v0)+CRC)*C'*inv(R)*(y(:,t)-C*x0);
        Fx2(:,t)=A*Fx1(:,t);
    else
    % here introduce temp1, temp2 to reduce duplicate computation
        v2=Fv2(:,:,t-1);
        temp1 = v2*CRC*inv(inv(v2)+CRC);     
        Fv1(:,:,t)=v2-v2*CRC*v2+temp1*CRC*v2;
        Fv2(:,:,t)=A*Fv1(:,:,t)*A'+Q;
        x2=Fx2(:,t-1);
        temp2 = C'*invR * (y(:,t)-C*x2);
        Fx1(:,t)=x2 + v2*temp2-temp1*temp2;
    end
end

%   Smoother
for i=1:T
    t=T-i+1;
    if t==T
        Sv(:,:,t)=Fv1(:,:,t);
        Jt=Fv1(:,:,t)*A'/(Fv2(:,:,t));
        Scov(:,:,t)=Fv2(:,:,t)*Jt';
        Sx(:,t)=Fv1(:,t)+Jt*(Fv2(:,t)-A*Fv1(:,t));
    else
        Jt=Fv1(:,:,t)*A'/(Fv2(:,:,t));
        Sv(:,:,t)=Fv1(:,:,t)+Jt*(Sv(:,:,t+1)-Fv2(:,:,t))*Jt';
        Scov(:,:,t)=Sv(:,:,t+1)*Jt';
        Sx(:,t)=Fx1(:,t)+Jt*(Sx(:,t+1)-A*Fx1(:,t));
    end
end

end
