# Modeling the long term response to selection under the Breeder's Equation and different trait architecture
# See Coop for single generation simulations: https://github.com/cooplab/popgen-notes/blob/master/Rcode/Quant_gen/QT3.R

library("shape")

# recombination shuffles the loci
recomb = TRUE
# samples
N = 5e3

allele.freq = 0.5
sel = 0.25
gens = 12
stop_selection = 10
seeds = 5
par(mfrow=c(2,4))

for ( M in c(10,1000) ) {
for ( hsq in c(0.05,0.6) ) {


all_parental_mean = matrix(NA,nrow=seeds,ncol=gens)
all_offspring_mean = matrix(NA,nrow=seeds,ncol=gens)
all_hsq = matrix(NA,nrow=seeds,ncol=gens)
all_expected_mean = matrix(NA,nrow=seeds,ncol=gens)

for ( s in 1:seeds ) {
# generate parents
mum.hap.1<-matrix(rbinom(N*M,1,allele.freq),nrow=N)
mum.hap.2<-matrix(rbinom(N*M,1,allele.freq),nrow=N)
dad.hap.1<-matrix(rbinom(N*M,1,allele.freq),nrow=N)
dad.hap.2<-matrix(rbinom(N*M,1,allele.freq),nrow=N)
betas = rnorm(M,0,1)

for ( i in 1:gens ) {
  
  # generate genetic value
  mum.geno<-mum.hap.1+mum.hap.2
  mum.gv = mum.geno %*% betas
  dad.geno<-dad.hap.1+dad.hap.2
  dad.gv = dad.geno %*% betas
  
  if ( i == 1 ) {
    genetic_var = var(c(mum.gv,dad.gv))
    genetic_mean = mean(c(mum.gv,dad.gv))
  }
  mum.gv = (mum.gv - genetic_mean)/sqrt(genetic_var)
  dad.gv = (dad.gv - genetic_mean)/sqrt(genetic_var)
  
  mum.y = sqrt(hsq) * mum.gv + rnorm(N,0,sqrt(1-hsq))
  dad.y = sqrt(hsq) * dad.gv + rnorm(N,0,sqrt(1-hsq))
  
  cur_hsq = cor(c(mum.y,dad.y),c(mum.gv,dad.gv))^2
  #cat( mean(c(mum.y,dad.y)) , cur_hsq , '\n' )
  
  all_hsq[s,i] = cur_hsq
  all_offspring_mean[s,i] = mean(c(mum.y,dad.y))
  
  # select on phenotype
  if ( i >= stop_selection ) {
    top.sel.per.mums<- sample( length(mum.y) )
    top.sel.per.dads<- sample( length(dad.y) )
  } else { 
    top.sel.per.mums<- sample(which(mum.y>quantile(mum.y,p=1-sel)),N,replace=T)
    top.sel.per.dads<- sample(which(dad.y>quantile(dad.y,p=1-sel)),N,replace=T)
  }
  all_parental_mean[s,i] = mean(c(mum.y[top.sel.per.mums],dad.y[top.sel.per.dads]))
  
  dad.hap.1 = dad.hap.2[top.sel.per.dads,]
  dad.hap.2 = mum.hap.1[top.sel.per.mums,]
  mum.hap.1 = dad.hap.1[top.sel.per.dads,]
  mum.hap.2 = mum.hap.2[top.sel.per.mums,]
  
  if ( recomb ) {
    for ( m in 1:M ) {
      dad.hap.1[,m] = dad.hap.2[sample(N),m]
      dad.hap.2[,m] = mum.hap.1[sample(N),m]
      mum.hap.1[,m] = dad.hap.1[sample(N),m]
      mum.hap.2[,m] = mum.hap.2[sample(N),m]      
    }
  }
  
  # expected response
  S = mean((mum.y)[which(mum.y>quantile(mum.y,p=1-sel))]) - mean((mum.y))
  all_expected_mean[s,i] = mean(mum.y) + cur_hsq * S
  #cat( "Expected:" , mean(mum.y) + cur_hsq * S , '\n' )
}
}

expected_mean = apply(all_expected_mean,2,mean)
offspring_mean = apply(all_offspring_mean,2,mean)
all_hsq = apply(all_hsq,2,mean)
parental_mean = apply(all_parental_mean,2,mean)

plot( 1:gens , offspring_mean , type="n" , las=1, xlim=c(0,gens),ylim=c(0,10) , pch=19, bty="n" , xlab="Generations" , ylab="Phenotype Mean",main=paste("Initial h2 = ",hsq,", # SNPs = ",M,sep=''))
abline(v=stop_selection,lty=3)
lines( 2:(stop_selection) , expected_mean[1:(stop_selection-1)] , col="gray" , lwd=6 )
points( 2:gens , offspring_mean[2:gens] , pch=19 )
Arrows(x0=1:(stop_selection-1),y0=parental_mean[1:(stop_selection-1)],x1=2:stop_selection,y1=offspring_mean[2:stop_selection],lwd=2,col="#9ecae1",arr.length=0.15,arr.adj=1,arr.type = "triangle")
points( 1:(stop_selection-1) , parental_mean[1:(stop_selection-1)] , pch=19, col="#3182bd" )
legend("topleft",legend=c("S*h^2","Selected Parents","Offspring"),lwd=c(6,NA,NA),pch=c(NA,19,19),col=c("gray","#3182bd","black"),bty="n")

#lines( 1:gens , parental_mean , las=1, type="l" , lwd=3 , lty=3, col="royalblue", bty="n" , xlab="Generations" , ylab="Mean")
plot( 1:gens , all_hsq , type="l", las=1, lwd=3 , bty="n" , xlab="Generations" , ylab="Heritability (h2)",ylim=c(0,0.6),xlim=c(0,gens))
abline(v=stop_selection,lty=3)
points( 1:gens , all_hsq , pch=19 )

}
}