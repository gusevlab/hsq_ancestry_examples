library("pheatmap")
library("RColorBrewer")
library(gridExtra)

clr = brewer.pal(3,"Set1")

M = 10000
p0 = rep(0.5,M)
p_a = p0 + rnorm(M,0,0.1)
p_a[ p_a < 0 ] = 0
p_a[ p_a > 1 ] = 1
p_b = p0 + rnorm(M,0,0.1)
p_b[ p_b < 0 ] = 0
p_b[ p_b > 1 ] = 1

fst = ((p_a - p_b)^2) / ( p_a*(1-p_b) + p_b*(1-p_a))
mean(((p_a - p_b)^2))/mean(( p_a*(1-p_b) + p_b*(1-p_a)))

all_pa = p_a
all_pb = p_b

par(mfrow=c(1,2))
N_1 = 10
for ( N_2 in c(20,1010) ) {
m_val = c(10,100,1000,10e3)
plot( 0 , 0 , type="n" , main=paste("N1 = 10 , N2 = ",N_2-N_1,sep=''),xlim=c(0.9,110e3) , ylim=c(-0.5,0.5) , log="x" , xlab="# SNPs" , xaxt="n" , ylab="PC1",las=1,bty="n")
axis( 1 , at=m_val,lab=m_val)

for ( M in m_val ) {
  p_a = all_pa[1:M]
  p_b = all_pb[1:M]
  mat = matrix(NA,nrow=N_2,ncol=M)
  for ( i in 1:N_1 ) {
    mat[i,] = rbinom(M,2,p_a)
  }
  
  for ( i in (N_1+1):(N_2) ) {
    mat[i,] = rbinom(M,2,p_b)
  }
  
  mat = mat[ , apply(mat,2,sd) != 0 ]
  mat = scale(mat)
  
  pop_lbl = c(rep(1,N_1),rep(0,N_2-N_1))
  pca = prcomp(t(mat),center=T,scale=T)
  pc1 = pca$rotation[,1]
  if ( cor(pc1,pop_lbl) < 0 ) pc1 = pc1*-1
  points( jitter(rep(M,length(pop_lbl)))[pop_lbl==0] , pc1[pop_lbl==0] , pch=19 , col = paste(clr[1],"50",sep='') , cex=0.75 )
  points( jitter(rep(M,length(pop_lbl)))[pop_lbl==1] , pc1[pop_lbl==1] , pch=19 , col = paste(clr[2],"50",sep='') , cex=0.75 )
}
}

# plot heatmaps
N_2 = 20
plot_list=list()
ctr = 1
for ( M in m_val ) {
  p_a = all_pa[1:M]
  p_b = all_pb[1:M]
  mat = matrix(NA,nrow=N_2,ncol=M)
  for ( i in 1:N_1 ) {
    mat[i,] = rbinom(M,2,p_a)
  }
  
  for ( i in (N_1+1):(N_2) ) {
    mat[i,] = rbinom(M,2,p_b)
  }
  
  mat = mat[ , apply(mat,2,sd) != 0 ]
  mat = scale(mat)
  
  K = (mat) %*% t(mat) / (M)
  diag(K) = NA
  plot_list[[ paste("p",ctr,sep='') ]] = pheatmap(K,cluster_rows=F,cluster_cols=F)[[4]]
  ctr = ctr + 1
}
grid.arrange(grobs=plot_list, ncol=4)