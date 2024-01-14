# Visualizing the relationship between Fst and heritable trait differences
# Warning lots of runs by default, will be slow

library(RColorBrewer)

par(mfrow=c(1,4))

for ( M in c(100,1e3) ) {

seeds = 1e3
N = 100
p0 = 0.5

gen_power = (seq(0,6,by=0.25))
gens = seq(1000,8000,by=1000)
  
hsq = 0.2
Ne = 10e3
maf_anc = rep(p0,M)
drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
all_drift = drift_var

pop_lbl = c(rep(1,N),rep(0,N))
mean_fst_hudson = matrix(NA,nrow=1,ncol=length(all_drift))
mean_fst_nei73 = matrix(NA,nrow=1,ncol=length(all_drift))
mean_rsq = matrix(NA,nrow=1,ncol=length(all_drift))
se_rsq = matrix(NA,nrow=1,ncol=length(all_drift))
mean_gap = matrix(NA,nrow=1,ncol=length(all_drift))

for( i in 1:length(all_drift) ) {
  
  drift = all_drift[i]
  fst_hudson = rep(NA,seeds)
  fst_nei73 = rep(NA,seeds)
  rsq = rep(NA,seeds)
  gap = rep(NA,seeds)
  
  for ( s in 1:seeds ) {
    maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift) )
    maf_pop2 = maf_anc + rnorm(M,0,sqrt(drift) )
    
    maf_pop2[ maf_pop2 > 1 ] = 1
    maf_pop2[ maf_pop2 < 0 ] = 0
    
    maf_pop1[ maf_pop1 > 1 ] = 1
    maf_pop1[ maf_pop1 < 0 ] = 0
    
    # hudson fst
    cur_fst = ((maf_pop1 - maf_pop2)^2) / ( maf_pop1*(1-maf_pop2) + maf_pop2*(1-maf_pop1))
    fst_hudson[s] = mean(cur_fst,na.rm=T)
    
    # nei fst (Edge)
    h_s = apply( rbind( maf_pop1^2 + (1-maf_pop1)^2 , maf_pop2^2 + (1-maf_pop2)^2 ) , 2 , mean)
    h_t = apply( rbind( apply( rbind( maf_pop1 , maf_pop2 ) , 2 , mean ) , apply( rbind( 1 - maf_pop1 , 1 - maf_pop2 ) , 2 , mean ) )^2 , 2 , sum )
    cur_fst = (h_s - h_t) / (1 - h_t )
    fst_nei73[s] = mean(cur_fst,na.rm=T)
    
    # generate heritable phenotype
    b = rnorm(M,0,1)
    x = matrix(NA,nrow=N*2,ncol=M)
    for ( snp in 1:M ) {
      x[1:N,snp] = rbinom(N,2,maf_pop1[snp])
      x[(N+1):(N+N),snp] = rbinom(N,2,maf_pop2[snp])
    }
    x = scale(x)
    x[ , apply(is.na(x),2,sum) != 0 ] = 0
    y1 = sqrt(hsq) * scale(x %*% b) + rnorm(N,0,sqrt(1-hsq))
    rsq[s] = summary(lm( y1 ~ pop_lbl ))$adj.r.sq
    y1 = scale(y1) * (3) + 70
    gap[s] = mean(y1[pop_lbl==0]) - mean(y1[pop_lbl==1])
  }
  mean_fst_hudson[1,i] = mean(fst_hudson,na.rm=T)
  mean_fst_nei73[1,i] = mean(fst_nei73,na.rm=T)
  mean_rsq[1,i] = mean(rsq,na.rm=T)
  se_rsq[1,i] = sd(rsq,na.rm=T)/sqrt(sum(!is.na(rsq)))
  mean_gap[1,i] = mean(abs(gap),na.rm=T)
  cat( gens[i] , mean_fst_hudson[1,i] , mean_fst_nei73[1,i] ,mean_rsq[1,i] , '\n' )
}


clr = brewer.pal(3,"Set2")
plot(gens,mean_fst_hudson[1,],xlim=c(0,max(gens)),ylim=c(0,max(mean_fst_hudson)),type="l",lwd=2,col=clr[1],las=1,bty="n",xlab="Generations",ylab="Fst or Rsq",main=paste(M," SNPs",sep=''))
points(gens,mean_fst_hudson[1,],pch=19,col=clr[1],cex=0.75)
lines(gens,mean_fst_nei73[1,],lwd=2,col=clr[2])
points(gens,mean_fst_nei73[1,],pch=19,col=clr[2],cex=0.75)
lines(gens,mean_rsq[1,],lwd=2,col=clr[3])
points(gens,mean_rsq[1,],pch=19,col=clr[3],cex=0.75)
legend("topleft",legend=c("Fst Hudson","Fst Nei","Rsq"),col=clr,bty="n",lwd=3)

label = which.min( abs(mean_fst_nei73 - 0.15) )
text( gens[label] , mean_rsq[label] , sprintf("%1.3f",mean_rsq[label]) , col=clr[3] , pos=3 )

max_val = 0.04 # max(mean_rsq)
plot( mean_fst_nei73*hsq , mean_rsq , xlab="Fst Nei * Heritability" , ylab="Rsq Observed", xlim=c(0,max_val) , ylim=c(0,max_val),las=1,bty="n",pch=19,cex=0.75)
abline(0,1,lty=3)
arrows( mean_fst_nei73*hsq , mean_rsq - 1.96*se_rsq , mean_fst_nei73*hsq , mean_rsq + 1.96*se_rsq , len=0 , col=clr[3] , lwd=5 )
points( mean_fst_nei73*hsq , mean_rsq , pch=19 , cex=0.75)

}