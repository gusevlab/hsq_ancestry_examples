# PCA phase transition

library("RColorBrewer")

M = 5e3
p0 = 0.5
Ne = 10e3
n_steps = c(seq(5,50,by=5),seq(60,100,by=10),200)
seeds = 10

all_gens = c(63,89,124)
fst_est = rep(NA,length(all_gens))
fst_seeds = 10
for ( g in 1:length(all_gens) ) {
  gens = all_gens[g]
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  maf_anc = rep(p0,M)
  all_fst = rep(NA,fst_seeds)
  for ( s in 1:fst_seeds ) {
    maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift_var) )
    maf_pop2 = maf_anc + rnorm(M,0,sqrt(drift_var) )
    maf_pop2[ maf_pop2 > 1 ] = 1
    maf_pop2[ maf_pop2 < 0 ] = 0
    maf_pop1[ maf_pop1 > 1 ] = 1
    maf_pop1[ maf_pop1 < 0 ] = 0
    # hudson fst
    fst = ((maf_pop1 - maf_pop2)^2) / ( maf_pop1*(1-maf_pop2) + maf_pop2*(1-maf_pop1))
    # neu fst
    h_s = apply( rbind( maf_pop1^2 + (1-maf_pop1)^2 , maf_pop2^2 + (1-maf_pop2)^2 ) , 2 , mean)
    h_t = apply( rbind( apply( rbind( maf_pop1 , maf_pop2 ) , 2 , mean ) , apply( rbind( 1 - maf_pop1 , 1 - maf_pop2 ) , 2 , mean ) )^2 , 2 , sum )
    fst = (h_s - h_t) / (1 - h_t )
    
    # ---
    fst = mean(fst,na.rm=T)
    all_fst[s] = fst
  }
  fst_est[g] = mean(all_fst)
}
bpp_thresh = as.integer((1 / (fst_est^2))/M/2)
cat( bpp_thresh , '\n' , sep='\t')

all_rsq = matrix(NA,nrow=length(all_gens),ncol=length(n_steps))
all_rsq_se = matrix(NA,nrow=length(all_gens),ncol=length(n_steps))

for ( g in 1:length(all_gens) ) {
  gens = all_gens[g]
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  
  for ( N_i in 1:length(n_steps) ) {
    N = n_steps[N_i]
    
    if ( N <= 50 ) {
      seeds = 2000/N
    } else {
      seeds = 10
    }
    
    pop_lbl = c(rep(1,N),rep(0,N))
    maf_anc = rep(p0,M)
    rsq = rep(NA,seeds)
    
    for ( s in 1:seeds ) {
    maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift_var) )
    maf_pop2 = maf_anc + rnorm(M,0,sqrt(drift_var) )
    
    maf_pop2[ maf_pop2 > 1 ] = 1
    maf_pop2[ maf_pop2 < 0 ] = 0
    
    maf_pop1[ maf_pop1 > 1 ] = 1
    maf_pop1[ maf_pop1 < 0 ] = 0
    
    x = matrix(NA,nrow=N*2,ncol=M)
    
    for ( snp in 1:M ) {
      x[1:N,snp] = rbinom(N,2,maf_pop1[snp])
      x[(N+1):(N+N),snp] = rbinom(N,2,maf_pop2[snp])
    }
    x = scale( x[ , apply(x,2,sd) != 0 ] )
    
    pca = prcomp(t(x),center=F,scale=F)
    rsq[ s ] = summary( lm(pop_lbl ~ pca$rotation[,1] ) )$adj.r.sq
    # cat( gens , N , rsq[ s ] , '\n' )
    }
    all_rsq[g,N_i] = mean(rsq,na.rm=T)
    all_rsq_se[g,N_i] = sd(rsq,na.rm=T)/sqrt(sum(!is.na(rsq)))
    cat( gens , N , all_rsq[g,N_i] , all_rsq_se[g,N_i] , '\n' )
  }
}

clr = tail(brewer.pal(6,"Blues"),n=length(all_gens))

plot(0,0,type="n",xlim=c(0,max(n_steps)) , ylim=c(0,1),xlab="# Samples",ylab="R2 with population label",bty="n",las=1)
for ( g in 1:length(all_gens) ) {
  polygon( c(n_steps , rev(n_steps)) , c((all_rsq[g,] - 1.96*all_rsq_se[g,]),rev(all_rsq[g,] + 1.96*all_rsq_se[g,])) , col=paste(clr[g],"50",sep='') , border=NA )
  lines( n_steps , (all_rsq[g,]) , col=clr[g], lwd=2 )
  pt = which.min( abs(n_steps - bpp_thresh[g]) )
  points( n_steps[pt] , (all_rsq[g,pt]) , pch=19 , col=clr[g] , cex=0.75)
}
legend( "bottomright" , legend=paste("Fst = ",sprintf("%1.4f",fst_est),sep='') , col=clr, lwd=2 , bty='n' )
