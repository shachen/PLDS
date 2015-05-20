function [prodvec]=mat_vec_prod(A,p,y)
% mat_vec_prod calculates the product of a big block diagonal matrix 
% and a vector. the matrix is of format bigA=kron(eye(p),A)
% prodvec = bigA * y
[r,c]=size(A);
newy=reshape(y,c,p);
prodmat = A * newy;
prodvec = reshape(prodmat,r*p,1);
end