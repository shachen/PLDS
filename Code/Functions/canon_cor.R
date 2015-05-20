# This function calculate the canonical correlation in polynomial time O(n^4)
# The Hungurian algorithm solves the LSAP problem in polynomial time
library(combinat)
require(clue)
pMatrix.min <- function(A, B) {
# finds the permutation P of A such that ||PA - B|| is minimum in Frobenius norm
# In Frobenious norm ||PA-I||^2 = Trace[(PA-I)(PA-I)^T]=Trace(AA^T)-2Trace(PA)+Trace(I) 
# Uses the linear-sum assignment problem (LSAP) solver in the "clue" package
# Returns P%*%A and the permutation vector `pvec' such that 
# A[pvec, ] is the permutation of A closest to B
	n <- nrow(A)
	D <- matrix(NA, n, n)
	for (i in 1:n) {
	for (j in 1:n) {
	D[j, i] <- (sum((B[j, ] - A[i, ])^2))
	} }
vec <- c(solve_LSAP(D))
list(A=A[vec,], pvec=vec)
}

##########
#An example:
#A <- matrix(sample(1:25, size=25, rep=FALSE),  5,  5)
#B <- diag(1, nrow(A)) # this choice of B maximizes the trace of permuted A
#X <- pMatrix.min(A,B)
#sum(diag(A))  # original square matrix
#sum(diag(X$A))  # permuted A such that its trace is maximum among all permutations
#perm = permn(1:nrow(A))
#maxdiag = 0;
#for(i in 1:length(perm)){
#	tmpx = A[perm[[i]],];
#	tmpdiag = sum(diag(tmpx));
#	if (tmpdiag > maxdiag){
#		maxdiag = tmpdiag;
#	}
#}
##########
