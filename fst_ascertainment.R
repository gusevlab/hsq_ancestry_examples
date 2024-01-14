# Simulating Fst estimates under different SNP ascertainments

library(RColorBrewer)

par(mfrow=c(1,4))

# generate two source populations
seeds = 200
M = 100

for ( p0 in c(0.1,0.5) ) {
  gen_power = (seq(0,6,by=0.25))
  gens = sort( c(seq(10,10e3,by=1e3) , (10^gen_power)[-1]) )
  gens = gens[gens>10]
  
  Ne = 10e3
  #p0 = 0.5
  maf_anc = rep(p0,M)
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  all_drift = drift_var
  
  clr = brewer.pal(3,"Blues")[2:3]
  all_maf = c(0,0.01)
  
  mean_fst_hudson = matrix(NA,nrow=2,ncol=length(all_drift))
  mean_fst_nei73 = matrix(NA,nrow=2,ncol=length(all_drift))
  
  for ( m in 1:length(all_maf) ) {
    
    min_maf = all_maf[m]
    for( i in 1:length(all_drift) ) {
      
      drift = all_drift[i]
      fst_hudson = rep(NA,seeds)
      fst_nei73 = rep(NA,seeds)
      
      for ( s in 1:seeds ) {
        maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift) )
        maf_pop1[ maf_pop1 > 1 ] = 1
        maf_pop1[ maf_pop1 < 0 ] = 0
        
        keep = maf_pop1 >= min_maf & maf_pop1 <= (1-min_maf)
        #keep = maf_pop1 <= min_maf | maf_pop1 >= (1-min_maf)
        
        # hudson fst
        cur_fst = ((maf_pop1 - maf_anc)^2) / ( maf_pop1*(1-maf_anc) + maf_anc*(1-maf_pop1))
        fst_hudson[s] = mean(cur_fst[keep],na.rm=T)
        
        # nei fst (Edge)
        h_s = apply( rbind( maf_pop1^2 + (1-maf_pop1)^2 , maf_anc^2 + (1-maf_anc)^2 ) , 2 , mean)
        h_t = apply( rbind( apply( rbind( maf_pop1 , maf_anc ) , 2 , mean ) , apply( rbind( 1 - maf_pop1 , 1 - maf_anc ) , 2 , mean ) )^2 , 2 , sum )
        cur_fst = (h_s - h_t) / (1 - h_t )
        fst_nei73[s] = mean(cur_fst[keep],na.rm=T)
        
      }
      mean_fst_hudson[m,i] = mean(fst_hudson,na.rm=T)
      mean_fst_nei73[m,i] = mean(fst_nei73,na.rm=T)
      cat( drift , '\n' )
    }
  }
  
  max = 10e3
  for( is_log in c(F,T) ) {
    if ( is_log ) {
      plot(0,0,type="l",log="x",lwd=3,las=2,ylab="Fst",xlab="Generations (10^x)",bty="n",ylim=c(0,0.5),xlim=range(gens),xaxt="n",main=paste("p_anc = ",p0,sep='') )
      axis(side=1,at=10^(1:6),labels=1:6)
      abline(v=max,lty=2)
    } else {
      plot(0,0,type="l",lwd=3,las=2,ylab="Fst",xlab="Generations (thousands)",bty="n",ylim=c(0,0.5),xlim=c(1,max),xaxt="n",main=paste("p_anc = ",p0,sep=''))
      axis(side=1,at=seq(0,max,by=max/10),labels=seq(0,max/1e3,by=max/10e3))
      abline(v=max,lty=2)
      legend("topleft",legend=c("All","Common","Fst Hudson","Fst Nei 1973"),lwd=3,col=c(clr,"black","black"),lty=c(1,1,1,2),bty="n")
      
    }
    for ( m in 1:2 ) {
      lines(gens,mean_fst_hudson[m,],lty=1,lwd=2 , col=clr[m] )
      lines(gens,mean_fst_nei73[m,] , lty=2 , lwd=2 , col=clr[m] )
    }  
  }
}