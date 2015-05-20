% this code generates data for simulation
% this code assumes that R is diagonal
% option = 1: sparse C, regular A; option = 2: smooth C, sparse A
% option = 3: smooth & sparse C, regular A;

function simdata_gen(p,d,T,maxiter,option)
% set seed for reproducibility
rng(100);

% dimensions
n=str2num(p);m=str2num(d);T=str2num(T);option=str2num(option);

% C
if option == 1 || option == 3
    C=randn(n,m);
    for i= 1:n
       for j = 1:m
           if abs(C(i,j)) < 1
               C(i,j) = 0;
           end
       end
    end    
    if option == 3
        for i = 1:m
            idx = find(C(:,i));
            vals = C(idx,i);
            vals = sort(vals);
            C(idx,i) = vals;
        end
    end
end
if option == 2
    C = zeros(n,m);
    k = n / m;
    for i = 1:m
        if i ~= 1
            C(1:(k*(i-1)),i) = -sort(abs(randn(k*(i-1),1)),'descend');
        end
        C((1+k*(i-1)):n,i) = sort(abs(randn(n-k*(i-1),1)),'ascend');
    end
    clear i j;
end


% R
R=2*spdiags(ones(n,1),0,n,n);

% Q
Q=eye(m);


% A
cond_num = 50;
P = randn(m,m);
[U,S,V]=svd(P);
S(S~=0)=linspace(cond_num,1,m);
P=U*S*V;
Aeigen=repmat([0.9,0.8],1,m);
Aeigen=Aeigen(1:m);
A=(P\diag(Aeigen))*P;
clear U S V Aeigen cond_num;
if option == 2
    A = A + diag(ones(m,1)*4) + diag(ones(m-1,1)*2,1) + diag(ones(m-1,1)*2,-1);
    for i = 1:m
        for j = 1:m
            if abs(A(i,j)) < 0.1
                A(i,j) = 0;
            end
        end
    end
    eigvals = eig(A);
    A = A/max(abs(eigvals))*0.95;
end

% x0
init_x=zeros(m,1);

% x y
x=zeros(m,T);
y=zeros(n,T);
sqrtRdiag=sqrt(diag(R));
for t=1:T
    if t==1
        w=mvnrnd(zeros(m,1),Q);
        x(:,t)=A*init_x+w';
    else
        w=mvnrnd(zeros(m,1),Q);
        x(:,t)=A*x(:,t-1)+w';
    end
    vtmp=normrnd(0,1,n,1);
    v=sqrtRdiag .* vtmp;
    y(:,t)=C*x(:,t)+v;
end
clear w v vtmp;

% starting values
a=eye(m);
c=randi(10,[n,m]);
q=eye(m);
r=spdiags(ones(n,1),0,n,n);
Pi=zeros(m,1);
v=eye(m);
tol = 10e-3;
miter = str2num(maxiter);

if option == 1
    save(strcat('../../../Data/Results/Simulation/p-',int2str(n),'-d-',int2str(m),'-T-',int2str(T),'-data-regA-sparseC.mat'));
end
if option == 2
    save(strcat('../../../Data/Results/Simulation/p-',int2str(n),'-d-',int2str(m),'-T-',int2str(T),'-data-sparseA-smoothC.mat'));
end
if option == 3
    save(strcat('../../../Data/Results/Simulation/p-',int2str(n),'-d-',int2str(m),'-T-',int2str(T),'-data-regA-smoothsparseC.mat'));
end

end
