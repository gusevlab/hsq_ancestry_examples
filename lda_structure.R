library("topicmodels")
library("RColorBrewer")

M = 200
N = 500
clr = "#3182bd"
seeds = 20

#title = c("a. STRUCTURE (uniform)","b. STRUCTURE (top 50%)","c. STRUCTURE (top 25%)","d. STRUCTURE (uniform)","e. STRUCTURE (top 50%)","f. STRUCTURE (top 25%)","g. PCA (uniform)","h. PCA (top 50%)","i. PCA (top 25%)")
title = c("a. Alpha=0.5, Uniform","b. Alpha=0.5, Top 50%","c. Alpha=0.5, Top 25%","d. Alpha=0.1, Uniform","e. Alpha=0.1, Top 50%","f. Alpha=0.1, Top 25%")
ctr = 0
par(mfrow=c(2,3))
for ( MODE in c("STRUCTURE_0.5","STRUCTURE_0.1") ) {
  for ( min_pct in c(0,0.5,0.75) ) {
    mean_buckets = seq(0,0.9,by=0.1)
    mean_values = matrix(nrow=seeds,ncol=length(mean_buckets))
    mean_rsq = rep(NA,seeds)
    
    for ( s in 1:seeds ) {
      # generate two populations and a continuum
      p0 = runif(M,0.1,0.9)
      p1 = p0 + rnorm(M,0,0.3)
      p1[ p1 < 0 ] = 0
      p1[ p1 > 1 ] = 1
      #fst = ((p0 - p1)^2) / ( p0*(1-p1) + p1*(1-p0))
      
      # continuum
      pct = c(rep(seq(0,0.9,by=0.1)+0.05,4),runif(N-40,min_pct,1))
      # sample each individual
      mat = matrix(NA,nrow=N,ncol=M)
      for ( i in 1:N ) {
        # get the weighted average frequency
        # p_cur = p0 * pct[i] + p1 * (1-pct[i])
        # sample some sites from each population
        p_cur = p0
        keep = sample( 1:M , M*(1-pct[i]) )
        p_cur[ keep ] = p1[ keep ]
        mat[i,] = rbinom(M,2,p_cur)
      }
      mat = mat[,apply(mat,2,sd)!=0]
      
      if ( MODE =="PCA") {
        s_mat = scale(mat)
        pca = prcomp(t(s_mat),center=F,scale=F)
        pc1 = pca$rotation[,1]
        if ( cor(pct,pc1) < 0 ) {
          pc1 = -1*pc1
        }
        if ( s == 1 ) {
          ctr = ctr + 1
          plot( pct , pc1 , type="n",las=1,ylim=c(-0.25,0.1),xlim=c(0,1),xlab="Genetic distance",ylab=MODE,main=title[ctr])
        }
        
      } else if ( MODE == "STRUCTURE_0.5" | MODE == "STRUCTURE_0.1" ) {
        if ( MODE == "STRUCTURE_0.5" ) lda_run = LDA(mat,k=2,method="Gibbs" , control = list(alpha = 0.5) )
        else if ( MODE == "STRUCTURE_0.1" ) lda_run = LDA(mat,k=2,method="Gibbs" , control = list(alpha = 0.1) )
        
        pc1 = posterior(lda_run)$topics[,1]
        if ( cor(pct,pc1) < 0 ) {
          pc1 = 1 - pc1
        }
        if ( s == 1 ) {
          ctr = ctr + 1
          plot( pct , pc1 , type="n",las=1,ylim=c(0,1),xlim=c(0,1),xlab="Genetic distance",ylab="STRUCTURE Pop",main=title[ctr])
        }
      }
      
      points( pct , pc1 , pch=1, col=paste(clr[1],"10",sep=''),cex=0.5)
            
      for ( rng in 1:length(mean_buckets) ){
        keep = pct >= mean_buckets[rng] & pct < mean_buckets[rng]+0.1
        mean_values[s,rng] = mean(pc1[keep])
      }  
      mean_rsq[s] = cor(pct,pc1)^2
    }
    legend("topleft",legend=paste("R2=",sprintf("%1.2f",mean(mean_rsq)),sep=''),bty="n")

    if( MODE == "STRUCTURE_0.5" | MODE == "STRUCTURE_0.1" ) {
      abline(0,1,lty=3)
    }
    for ( rng in 1:length(mean_buckets) ){
      keep = pct >= mean_buckets[rng] & pct < mean_buckets[rng]+0.1
      points( mean_buckets[rng] + 0.05 , mean(mean_values[,rng])  , pch=19 )
    }
    
  }
}
