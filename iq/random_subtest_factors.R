library("MASS")
library("RColorBrewer")

clr = rev(brewer.pal(6,"Blues"))

N = 500
M_val = c(10,50)
seeds = 50
rsq_val = c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.15)
rsq_est = matrix(NA,nrow=length(rsq_val),ncol=length(M_val))
factor_cor = rsq_est

for ( m in 1:length(M_val) ) {
  M = M_val[m]
  for ( r in 1:length(rsq_val) ) {
    rsq = rsq_val[r]
    sigma = matrix(sqrt(rsq),nrow=M,ncol=M)
    sigma_2 = matrix(sqrt(rsq),nrow=M*2,ncol=M*2)
    diag(sigma) = 1
    diag(sigma_2) = 1
    
    est = rep(NA,seeds)
    est_factor = rep(NA,seeds)
    for ( s in 1:seeds ) {
      mat = mvrnorm( N , mu = rep(0,M) , Sigma = sigma )
      pca = prcomp(t(scale(mat)),center=F,scale=F)
      est[s] = pca$sdev[1]^2 / sum(pca$sdev^2)
      # compute two factors
      mat = mvrnorm( N , mu = rep(0,M*2) , Sigma = sigma_2 )
      pca_a = prcomp(t(scale(mat[,1:M])),center=F,scale=F)
      pca_b = prcomp(t(scale(mat[,(M+1):(M+M)])),center=F,scale=F)
      est_factor[s] = cor( pca_a$rotation[,1] , pca_b$rotation[,1] )
    }
    rsq_est[r,m] = mean(est)
    factor_cor[r,m] = mean(est_factor)
    cat( rsq , rsq_est[r,m] , factor_cor[r,m] , '\n' )
  }
}

par(mfrow=c(1,2))
plot( rsq_val , rsq_est[,1] , type="o" , pch=19 , cex=0.5, lwd=2, bty="n" , col=clr[1] , las=1 , ylim=c(0,0.5),xlab="R2 between pairs of subtests", ylab="Variance explained by PC1")
lines( rsq_val , rsq_est[,2] , col=clr[2] ,lwd=2 )
points( rsq_val , rsq_est[,2] , col=clr[2] , pch=19 , cex=0.5 )
legend("topleft",legend=paste(M_val,"subtests"),col=clr[1:2],bty="n",pch=19)

plot( rsq_val , factor_cor[,1] , type="o" , cex=0.5, pch=19 , lwd=2, bty="n" , col=clr[1] , las=1 , ylim=c(0,1),xlab="R2 between pairs of subtests", ylab="Correlation of PC1 from two tests")
lines( rsq_val , factor_cor[,2] , col=clr[2] , lwd=2)
points( rsq_val , factor_cor[,2] , col=clr[2] , pch=19 , cex=0.5 )
