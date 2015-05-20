function [x]=fista(AtA,p,b,lambda,L,maxiter,tol)
%   version 2.1
%   the Fast Iterative Shrinkage-Thresholding Algorithm (fista)
%   solve the following problem:
%   argmin object(x) = x' * bigAtA * x -2 * b' * x + lambda * |x| = f + g
%   bigAtA is a block diagnal matrix, with p identical diagnol elements: AtA
%   L: Liptchiz constant for  f
%   maxiter,tol: max number of iterations; convergence tolerance

%   scale lambda
%   scale by the geometric mean of eigen values of AtA
eigvals = eig(AtA);
%lambda = lambda * exp(mean(log(abs(eigvals))));

[~,d] = size(AtA);
x=reshape(AtA\reshape(b,d,p),d*p,1);
y=x;
dif = 1;
iter = 0;
tau = 1 / L;
t = 1;
while (dif >= tol) && (iter <= maxiter)
    %iter
    %grad = 2 * A' * (A * y - b);
    %grad = 2 * mat_vec_prod(smallAtA,p,y) - 2 * bigAtb;
    grad = 2 * mat_vec_prod(AtA,p,y) - 2 * b;
    newx = (y - tau*grad - tau*lambda) .* ((y - tau*grad - tau*lambda)>0)...
        + (y-tau*grad + tau*lambda) .* ((y-tau*grad + tau*lambda) < 0);
    newt = (1 + sqrt(1+4*t^2))/2;
    y = newx + (t - 1)/newt * (newx - x);
    iter = iter + 1;
    dif = norm(newx-x);
    x = newx;
    t = newt;
end    
end
