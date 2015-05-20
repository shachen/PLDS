# A function for computing the Amari error of two matrices (true and estimated)

  amari=function(W.h,W,st=FALSE){

# standardize the two matrices to have column norm 1

	  if(st==TRUE){
		  for(i in 1:m){
		  W.h[,i]=W.h[,i]/sqrt(sum(W.h[,i]^2))
		  W[,i]=W[,i]/sqrt(sum(W[,i]^2))
  	  }}

	  P=abs(solve(W)%*%W.h)
	  max.row=apply(P,1,max)
	  max.col=apply(P,2,max)
	  sum.row=apply(P,1,sum)
	  sum.col=apply(P,2,sum)
		  
	  return((sum(sum.row/max.row-1)+sum(sum.col/max.col-1))/(2*m))
  }
