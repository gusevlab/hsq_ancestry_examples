# Library required for deltamethod computation of SE
library('msm')

# Utility function for calculating log of determinant
`logdet` <-
function(p) {
	det_l <- determinant(p,log=T)

	if ( det_l$sign[1] == 1 ) return(det_l$modulus[1])
	else return(det_l$modulus[1] * -1);
}

# Execute Haseman-Elston Regression

`HEreg` <-
function( A , y ) {
  r <- length(A)
  N <- length(y)
  yprod = y %*% t( y )
  yvals = yprod[ lower.tri(yprod,diag=F) ]
  
  
  # Matrix for all entries
  K_mat = matrix(NA,nrow=length(yvals),ncol=r)
  for ( i in 1:r ) {
    K_mat[,i] = (A[[i]])[ lower.tri(A[[i]],diag=F) ]
  }
  
  reg = lm( yvals ~ K_mat )
  return( list( "h2" = reg$coef[2] ) )
}

# Execute Average-Information ML (no fixed effects) until convergence
# A = list of GRM
# y = vector of phenotype entries
# Var = initial variance components (percent)

`aiML` <-
function( A, y , Var , verbose=TRUE ){

	r <- length(A) + 1
	N <- length(y)

	# Add matrix of residuals to A
	A[[r]] <- diag(N)

	AI <- matrix(0, ncol=r, nrow=r)
	S <- matrix(0, ncol=r, nrow=r)
	s <- matrix(0, ncol=1, nrow=r)

	l_dif <- 10
	it <- 0

	Var <- c(var(y)) * Var

	# Perform a single iteration of EM-based REML to initiate parameters
	if(verbose) cat("Performing a single iteration of EM-REML\n")
	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv

	if(verbose)
	cat("Prior values from EM-REML:")
	for ( i in 1:r ) {
		Var[i] <- (Var[i]^2 * t(y) %*% P %*% A[[i]] %*% P %*% y + sum(diag(Var[i]*diag(N) - Var[i]^2 * P %*% A[[i]])) )/N
		if(verbose)
		cat(" ",Var[i],sep='')
	}
	if(verbose)
	cat('\n')

	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv
	logL <- -0.5 * ( logdet(V) + t(y) %*% P %*% y )

	if(verbose)
	cat ("EM:\t",logL,'\n',sep='')

	# Iterate AI REML until convergence
	# while ( abs(l_dif) >= 10^-4 & it < 100 ){

	# ** GCTA style:
	while ( it < 100 & ( abs(l_dif) >= 10^-4 | (abs(l_dif)<10^-2 & l_dif < 0)) ){

		it <- it + 1

		# Average information matrix
		for ( i in 1:r ) {
			for( ii in 1:r ) {
				if ( i == r & ii == r ) AI[r,r] <- t(y) %*% P %*% P %*% P %*% y
				else if ( i == r ) AI[r,ii] <- t(y) %*% P %*% P %*% A[[ii]] %*% P %*% y
				else if ( ii == r ) AI[i,r] <- t(y) %*% P %*% A[[i]] %*% P %*% P %*% y
				else AI[i,ii] <- t(y) %*% P %*% A[[i]] %*% P %*% A[[ii]] %*% P %*% y
			}
		}
		AI <- 0.5*AI

		# Vector of first derivatives of log likelihood  function
		for ( i in 1:r ) {
			if ( i == r ) s[r,1] <- sum(diag(( P ))) - ( t(y) %*% P %*% P %*% y )
			else s[i,1] <- sum(diag(( P %*% A[[i]] ))) - ( t(y) %*% P %*% A[[i]] %*% P %*% y )
		}
		s <- -0.5*s

		# New variance components from AI and likelihood
		# Var <- Var + solve(AI) %*% s

		# ** GCTA style:
		if ( l_dif > 1 ) Var <- Var + 0.316*(solve(AI) %*% s)
		else Var <- Var + solve(AI) %*% s

		# Re-calculate V and P matrix
		V <- 0
		for ( i in 1:r ) V <- V + A[[i]] * Var[i]
		Vinv <- solve(V)
		P <- Vinv 

		# Likelihood of the MLM (ignoring constants)
		new_logL <- -0.5 * ( logdet(V) + t(y) %*% P %*% y )
		l_dif <- new_logL - logL
		logL <- new_logL

		if(verbose) {
		cat(it,'\t',logL,sep='')
		for( i in 1:r ) cat( '\t',Var[i],sep='' )
		cat('\n')
		}
	}
	if(verbose)
	cat('\n')

	# Calculate matrix for standard errors (same as AI matrix but w/out y)
	for( i in 1:r ) {
		for ( ii in 1:r ) {
			S[i,ii] <- sum(diag(P %*% A[[i]] %*% P %*% A[[ii]] ))
		}
	}
	S <- 0.5*S
	Sinv <- solve(S)

	if(verbose){
	for( i in 1:r ) cat( "V(G",i,")\t",Var[i],'\t',sqrt(Sinv[i,i]),'\n',sep="")
	}
	
	# Construct string equation for delta method "~x1+x2 ... +xn"
	SE.eq <- "~"
	for(i in 1:r) {
		if ( i == 1 ) SE.eq <- paste(SE.eq,"x",i,sep='')
		else SE.eq <- paste(SE.eq,"+x",i,sep='')
	}
	SE.p <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)

	if( verbose )
	cat( "Vp\t",sum(Var),'\t',SE.p,'\n',sep="")

	SE.i <- rep(0,r)
	
	for( i in 1:r ) {
		# Construct string equation for delta method "~xi/(x2 ... +xn)"
		SE.eq <- paste("~x",i,"/(",sep='')
		for( ii in setdiff(1:r,i) ) {
			if ( ii == setdiff(1:r,i)[1] ) SE.eq <- paste (SE.eq,"x",i,sep='')
			SE.eq <- paste (SE.eq,"+x",ii,sep='')
		}
		SE.eq <- paste(SE.eq,")",sep='')

		SE.i[i] <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)
		
		if(verbose)
		cat( "V(G",i,")/Vp\t",Var[i]/sum(Var),'\t',SE.i[i],'\n',sep='')
	}
	
	return( list( "h2" = Var[1]/sum(Var) , "se" = SE.i[1] ))
}

# Execute Average-Information REML until convergence
# A = list of GRM
# y = vector of phenotype entries
# X = matrix of fixed effects
# Var = initial variance components (percent)

`aiREML` <-
function( A, y, X , Var , verbose = FALSE ){

	r <- length(A) + 1
	N <- length(y)

	# Add matrix of residuals to A
	A[[r]] <- diag(N)

	AI <- matrix(0, ncol=r, nrow=r)
	S <- matrix(0, ncol=r, nrow=r)
	s <- matrix(0, ncol=1, nrow=r)

	l_dif <- 10
	it <- 0

	Var <- var(y)[1,] * Var

	# Perform a single iteration of EM-based REML to initiate parameters
	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv - Vinv %*% X %*% solve( t(X) %*% Vinv %*% X ) %*% t(X) %*% Vinv

	#if(verbose)
	#cat("Prior values from EM-REML:")
	#for ( i in 1:r ) {
	#	Var[i] <- (Var[i]^2 * t(y) %*% P %*% A[[i]] %*% P %*% y + sum(diag(Var[i]*diag(N) - Var[i]^2 * P %*% A[[i]])) )/N
	#	if(verbose)
	#	cat(" ",Var[i],sep='')
	#}
	#if(verbose)
	#cat('\n')

	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv - Vinv %*% X %*% solve( t(X) %*% Vinv %*% X ) %*% t(X) %*% Vinv
	logL <- -0.5 * ( logdet(V) + logdet(t(X) %*% Vinv %*% X) + t(y) %*% P %*% y )

	# Calculate using fixed effects
	# b <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
	# logL <- -0.5 * ( t(y - X * b) %*% Vinv %*% (y - X * b) + logdet(V) + logdet(t(X) %*% Vinv %*% X) )

	#if(verbose)
	#cat ("EM:\t",logL,'\n',sep='')

	# Iterate AI REML until convergence
	# while ( abs(l_dif) >= 10^-4 & it < 100 ){

	# ** GCTA style:
	#while ( it < 100 & ( abs(l_dif) >= 10^-4 | (abs(l_dif)<10^-2 & l_dif < 0)) ){
        while ( it < 50 & ( abs(l_dif) >= 10^-4 | (abs(l_dif)<10^-2 & l_dif < 0)) ){

		it <- it + 1

		# Average information matrix
		for ( i in 1:r ) {
			for( ii in 1:r ) {
				if ( i == r & ii == r ) AI[r,r] <- t(y) %*% P %*% P %*% P %*% y
				else if ( i == r ) AI[r,ii] <- t(y) %*% P %*% P %*% A[[ii]] %*% P %*% y
				else if ( ii == r ) AI[i,r] <- t(y) %*% P %*% A[[i]] %*% P %*% P %*% y
				else AI[i,ii] <- t(y) %*% P %*% A[[i]] %*% P %*% A[[ii]] %*% P %*% y
			}
		}
		AI <- 0.5*AI

		# Vector of first derivatives of log likelihood  function
		for ( i in 1:r ) {
			if ( i == r ) s[r,1] <- sum(diag(( P ))) - ( t(y) %*% P %*% P %*% y )
			else s[i,1] <- sum(diag(( P %*% A[[i]] ))) - ( t(y) %*% P %*% A[[i]] %*% P %*% y )
		}
		s <- -0.5*s

		# New variance components from AI and likelihood
		# Var <- Var + solve(AI) %*% s

		# ** GCTA style:
		if ( l_dif > 1 ) Var <- Var + 0.316*(solve(AI) %*% s)
		else Var <- Var + solve(AI) %*% s

		# Re-calculate V and P matrix
		V <- 0
		for ( i in 1:r ) V <- V + A[[i]] * Var[i]
		Vinv <- solve(V)
		P <- Vinv - Vinv %*% X %*% solve( t(X) %*% Vinv %*% X ) %*% t(X) %*% Vinv

		# Likelihood of the MLM (ignoring constants)
		new_logL <- -0.5 * ( logdet(V) + logdet(t(X) %*% Vinv %*% X) + t(y) %*% P %*% y )
		l_dif <- new_logL - logL
		logL <- new_logL

		if(verbose) {
		cat(it,'\t',logL,sep='')
		for( i in 1:r ) cat( '\t',Var[i],sep='' )
		cat('\n')
		}
	}
	#if(verbose)
	#cat('\n')

        if (1 == 2) {
                                        # Calculate matrix for standard errors (same as AI matrix but w/out y)
          for( i in 1:r ) {
            for ( ii in 1:r ) {
              S[i,ii] <- sum(diag(P %*% A[[i]] %*% P %*% A[[ii]] ))
            }
          }
          S <- 0.5*S
          Sinv <- solve(S)
        
        
          #if(verbose){
          #  for( i in 1:r ) cat( "V(G",i,")\t",Var[i],'\t',sqrt(Sinv[i,i]),'\n',sep="")
          #}
	
        
                                        # Construct string equation for delta method "~x1+x2 ... +xn"
          SE.eq <- "~"
          for(i in 1:r) {
            if ( i == 1 ) SE.eq <- paste(SE.eq,"x",i,sep='')
            else SE.eq <- paste(SE.eq,"+x",i,sep='')
          }
          SE.p <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)
          if(verbose) {
            cat( "Vp\t",sum(Var),'\t',SE.p,'\n',sep="")
          }
          
          SE.i <- rep(0,r)
          
          for( i in 1:r ) {
                                        # Construct string equation for delta method "~xi/(x2 ... +xn)"
            SE.eq <- paste("~x",i,"/(",sep='')
            for( ii in setdiff(1:r,i) ) {
              if ( ii == setdiff(1:r,i)[1] ) SE.eq <- paste (SE.eq,"x",i,sep='')
              SE.eq <- paste (SE.eq,"+x",ii,sep='')
            }
            SE.eq <- paste(SE.eq,")",sep='')
            
            SE.i[i] <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)
            if(verbose){
              cat( "V(G",i,")/Vp\t",Var[i]/sum(Var),'\t',SE.i[i],'\n',sep='')
            }
          }
        }
        
	#return( list( "h2" = Var[1]/sum(Var) , "se" = SE.i[1], "vars" = Var))
        return( list( "h2" = Var[1]/sum(Var) , "vars" = Var,"logL"=logL))
}
