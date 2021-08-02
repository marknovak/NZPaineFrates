# Classic Jaccard index
# input is frequency
J.class<-function(i,j){
  rem <- i==0 & j==0 | is.na(i) & is.na(j)
  i <- i[!rem]
  j <- j[!rem]
  if(length(i) == 0 | length(j) == 0){
    return(NA)
    }
  else{
	a <- sum(i>0 & j>0)
	if(a==0){return(0)}
	else
		b <- length(setdiff(which(i>0), which(j>0)))
		c <- length(setdiff(which(j>0), which(i>0)))
	return(a/(a+b+c))
		}
}

# Jaccard index based on abundance - Chao et al. (2005, Ecology Letters)
# input is frequency
J.abd<-function(i,j){
  rem <- i==0 & j==0 | is.na(i) & is.na(j)
  i <- i[!rem]
  j <- j[!rem]
  if(length(i) == 0 | length(j) == 0){
    return(NA)
  }
  else{
    i <- prop.table(i)
    j <- prop.table(j)
  	O <- sum(i>0 & j>0)
  	if(O==0){return(0)}
  	else
  		U <- sum(i[(i>0 & j>0)])
  		V <- sum(j[(i>0 & j>0)])
  	return(U*V/(U+V-U*V))
  }
}

# Jaccard estimator based on abundance - Chao et al. (2005, Ecology Letters)
# input is frequency
II<-function(x){ifelse(x==1,1,0)} # Indicator function

J.est<-function(i,j){
  rem <- i==0 & j==0 | is.na(i) & is.na(j)
  i <- i[!rem]
  j <- j[!rem]
  if(length(i) == 0 | length(j) == 0){
    return(NA)
  }
  else{
	O<-sum(i>0 & j>0)
	if(O==0){return(0)}
	else
		D12 <- which(i>0 & j>0)
		n <- sum(i); m<-sum(j)
		f_1p <- sum(i==1 & j>=1)
		f_2p <- max(c(sum(i==2 & j>=1), 1))
		f_p1 <- sum(i>=1 & j==1)
		f_p2 <- max(c(sum(i>=1 & j==2), 1))
		U.est <- min(sum(i[D12]/n) + ((m-1)/m) * f_p1/(2*f_p2) * sum((i[D12]/n)*II(j[D12])) ,1)
		V.est <- min(sum(j[D12]/m) + ((n-1)/n) * f_1p/(2*f_2p) * sum((j[D12]/m)*II(i[D12])) ,1)
	return(U.est*V.est/(U.est+V.est-U.est*V.est))
  }
}
