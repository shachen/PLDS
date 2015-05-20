function [x] = l2reg (A,b,lambda)
%   version 2.0
%   l2reg solves the following L2 regularization problem:
%   argmin f(x) = x' * BigA * x - 2 * b' * x + lambda * x' * x
%   where: 
%   BigA: block diagnol matrix, with all p diagonal elements equaling A
%   A: d by d positive definite matrix
%   b: d by p matrix, converted from a pd by 1 vector
%   lambda: a vector of length p, penalty parameters for different parts of
%   vector x
%   x: the optimizer

%   solve the problem
%   differentiate:
%   2 * BigA * x - 2 * b + 2 * lambda * I * x = 0
%   x = inv(BigA + lambda * I ) * b

%   testing
%   this code is compared to the version 1.0
%   and got the same results

%   scale lambda
%   scale by the geometric mean of eigen values of A
eigvals = eig(A);
%lambda = lambda * exp(mean(log(abs(eigvals))));
lambda = lambda * min(abs(eigvals(1,1)));

[d,p] = size(b);

x = zeros(d,p);
for i = 1:p
    x(:,i) = (A + lambda(i,1) * eye(d))\b(:,i);
end
x = reshape(x,p*d,1);
end