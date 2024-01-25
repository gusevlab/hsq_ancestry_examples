# Visualizing the relationship between Fst and heritable trait differences
# Warning lots of runs by default, will be slow

library(RColorBrewer)

num_markers = c(100,1e3)
par(mfrow=c(1,2*length(num_markers)))

seeds = 1e3
N = 500
hsq = 0.1
Ne = 10e3
gens = seq(1e3,8e3,by=1e3)

# uniform frequency spectrum for the ancestral population
maf_anc = runif(M,0.1,0.9)

for ( M in num_markers ) {
  
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  all_drift = drift_var
  
  pop_lbl = c(rep(1,N),rep(0,N))
  mean_fst_hudson = matrix(NA,nrow=1,ncol=length(all_drift))
  mean_fst_nei73 = matrix(NA,nrow=1,ncol=length(all_drift))
  mean_rsq = matrix(NA,nrow=1,ncol=length(all_drift))
  se_rsq = matrix(NA,nrow=1,ncol=length(all_drift))
  
  for( i in 1:length(all_drift) ) {
    
    drift = all_drift[i]
    fst_hudson = rep(NA,seeds)
    fst_nei73 = rep(NA,seeds)
    rsq = rep(NA,seeds)
    
    for ( s in 1:seeds ) {
      maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift) )
      maf_pop2 = maf_anc + rnorm(M,0,sqrt(drift) )
      
      maf_pop2[ maf_pop2 > 1 ] = 1
      maf_pop2[ maf_pop2 < 0 ] = 0
      
      maf_pop1[ maf_pop1 > 1 ] = 1
      maf_pop1[ maf_pop1 < 0 ] = 0
      
      # retain only sites that have not fixed to the same allele
      keep = maf_pop1*maf_pop2 != 0 & (1-maf_pop1)*(1-maf_pop2) != 0
      
      # Calculate hudson fst
      # --- standard calculation
      # cur_fst = ((maf_pop1 - maf_pop2)^2) / ( maf_pop1*(1-maf_pop2) + maf_pop2*(1-maf_pop1))
      # fst_hudson[s] = mean(cur_fst,na.rm=T)
      # --- average of ratios
      fst_hudson[s] = mean( ((maf_pop1[keep] - maf_pop2[keep])^2) ,na.rm=T) / mean(( maf_pop1[keep]*(1-maf_pop2[keep]) + maf_pop2[keep]*(1-maf_pop1[keep])),na.rm=T)
      
      # Calculate nei fst
      h_s = apply( rbind( maf_pop1[keep]^2 + (1-maf_pop1[keep])^2 , maf_pop2[keep]^2 + (1-maf_pop2[keep])^2 ) , 2 , mean)
      h_t = apply( rbind( apply( rbind( maf_pop1[keep] , maf_pop2[keep] ) , 2 , mean ) , apply( rbind( 1 - maf_pop1[keep] , 1 - maf_pop2[keep] ) , 2 , mean ) )^2 , 2 , sum )
      # --- standard calculation:
      # cur_fst = (h_s - h_t) / (1 - h_t )
      # fst_nei73[s] = mean(cur_fst,na.rm=T)
      # --- average of ratios
      fst_nei73[s] = mean((h_s - h_t),na.rm=T)/mean((1 - h_t ),na.rm=T)
      
      # fst_l from edge+rosenberg (equivalent to Hudson's Fst for diploids)
      # ---
      # deltasq = mean( ((maf_pop1[keep] - maf_pop2[keep])^2) )
      # pbar = mean(maf_pop1[keep])
      # qbar = mean(maf_pop2[keep])
      # ps2 = var(maf_pop1[keep])
      # qs2 = var(maf_pop2[keep])
      # numer = 2*deltasq
      # denom = 2*( pbar*(1-pbar) - ps2 + qbar*(1-qbar) - qs2) + 2*deltasq
      # fst_hudson[s] = numer / denom
      # ---
      
      # generate a heritable phenotype
      
      b = rnorm(M,0,1)
      x = matrix(NA,nrow=N*2,ncol=M)
      for ( snp in 1:M ) {
        # draw SNPs from the binomial with population-specific frequencies
        x[1:N,snp] = rbinom(N,2,maf_pop1[snp])
        x[(N+1):(N+N),snp] = rbinom(N,2,maf_pop2[snp])
      }
      # generate phenotype
      y1 = sqrt(hsq) * scale(x %*% b) + rnorm(N,0,sqrt(1-hsq))
      # estimate variance explained by pop label
      rsq[s] = summary(lm( y1 ~ pop_lbl ))$adj.r.sq
    }
    mean_fst_hudson[1,i] = mean(fst_hudson)
    mean_fst_nei73[1,i] = mean(fst_nei73)
    mean_rsq[1,i] = mean(rsq)
    se_rsq[1,i] = sd(rsq)/sqrt(sum(!is.na(rsq)))
    cat( gens[i] , mean_fst_hudson[1,i] , mean_fst_nei73[1,i] ,mean_rsq[1,i] , '\n' )
  }
  
  
  clr = brewer.pal(3,"Set2")
  plot(gens,mean_fst_hudson[1,],xlim=c(0,max(gens)),ylim=c(0,max(mean_fst_hudson)),type="l",lwd=2,col=clr[1],las=1,bty="n",xlab="Generations",ylab="Fst or Rsq",main=paste(M," SNPs",sep=''))
  points(gens,mean_fst_hudson[1,],pch=19,col=clr[1],cex=0.75)
  lines(gens,mean_fst_nei73[1,],lwd=2,col=clr[2])
  points(gens,mean_fst_nei73[1,],pch=19,col=clr[2],cex=0.75)
  lines(gens,mean_rsq[1,],lwd=2,col=clr[3])
  points(gens,mean_rsq[1,],pch=19,col=clr[3],cex=0.75)
  legend("topleft",legend=c("Fst Hudson","Fst Nei","R2 pop"),col=clr,bty="n",lwd=3)
  
  label = which.min( abs(mean_fst_hudson - 0.15) )
  text( gens[label] , mean_rsq[label] , sprintf("%1.3f",mean_rsq[label]) , col=clr[3] , pos=3 )
  
  # this is the expectation from Edge + Rosenberg
  expected_rsq = (mean_fst_hudson)*hsq
  
  max_val = max(c(mean_rsq,expected_rsq))
  plot( expected_rsq , mean_rsq , type="n",xlab="Fst Hudson * Heritability" , ylab="Observed R2", xlim=c(0,max_val) , ylim=c(0,max_val),las=1,bty="n",pch=19,cex=0.75)
  abline(0,1,lty=3)
  polygon( c(expected_rsq ,rev(expected_rsq)), c(mean_rsq - 1.96*se_rsq , rev(mean_rsq + 1.96*se_rsq)) , border=NA , col=paste(clr[3],"50",sep=''))
  points( expected_rsq , mean_rsq , pch=19 , cex=0.75,col="black")
  
}