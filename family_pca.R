set.seed(42)

# Num individuals
par( mfrow=c(1,4) )
# Num markers
M = 500
# Minor allele frequency
maf = 0.5

LABEL = TRUE

title_prefix = c("A","B","C")
all_N = c(1,5,50)

clr = c("#762a83","#c2a5cf","#5aae61")

for ( i in 1:3 ) {
  # generate some random individuals
  N = all_N[i]
  mate1_hap_a = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  mate1_hap_b = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  
  mate2_hap_a = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  mate2_hap_b = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  
  # generative children
  n_children = 5
  all_children_hap_a = matrix( NA , nrow=n_children , ncol=M )
  all_children_hap_b = matrix( NA , nrow=n_children , ncol=M )
  for ( n in 1:n_children ) {
    # generate offspring
    child_hap_a = mate1_hap_a
    sampled_vars_mate1 = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
    child_hap_a[ !sampled_vars_mate1 ] = mate1_hap_b[ !sampled_vars_mate1 ]
    
    child_hap_b = mate2_hap_a
    sampled_vars_mate2 = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
    child_hap_b[ !sampled_vars_mate2 ] = mate2_hap_b[ !sampled_vars_mate2 ]
    
    all_children_hap_a[n,] = child_hap_a[1,]
    all_children_hap_b[n,] = child_hap_b[1,]
  }
  
  all_labels = c( rep("children",n_children) , rep("parents",1) , rep("unrelated",N-1) , rep("parents",1) , rep("unrelated",N-1) )
  
  all_genos = rbind( all_children_hap_a + all_children_hap_b , mate1_hap_a + mate1_hap_b , mate2_hap_a + mate2_hap_b )
  
  # drop monomorphic
  keep = apply(all_genos,2,sd) != 0
  all_genos = all_genos[ , keep ]
  
  pc <- prcomp( t( scale(all_genos) ) , center = FALSE, scale. = FALSE )
  varexp = pc$sdev^2 / sum(pc$sdev^2)
  
  # visualize
  plot( pc$rotation[,"PC1"] , pc$rotation[,"PC2"] , type="n" , xlab="PC1" , ylab="PC2" , main=paste(title_prefix[i],": ",N*2," unrelated samples\nand 5 children",sep='') , las=1 )
  points(  pc$rotation[all_labels == "unrelated","PC1"] , pc$rotation[all_labels == "unrelated","PC2"] , pch=19 , col=paste(clr[3],"75",sep='') )
  points(  pc$rotation[all_labels == "unrelated","PC1"] , pc$rotation[all_labels == "unrelated","PC2"] , pch=1 , col=clr[3] )
  points(  pc$rotation[all_labels == "parents","PC1"] , pc$rotation[all_labels == "parents","PC2"] , pch=19 , col=clr[1] )
  points(  pc$rotation[all_labels == "parents","PC1"] , pc$rotation[all_labels == "parents","PC2"] , pch=1 , col=clr[1] )
  points(  pc$rotation[all_labels == "children","PC1"] , pc$rotation[all_labels == "children","PC2"] , pch=19 , col=paste(clr[2],"75",sep='') )
  points(  pc$rotation[all_labels == "children","PC1"] , pc$rotation[all_labels == "children","PC2"] , pch=1 , col=clr[2] )
  
}

pc <- prcomp( t( scale(all_genos[ ((n_children+1):nrow(all_genos)) ,]) ) , center = FALSE, scale. = FALSE )
varexp = pc$sdev^2 / sum(pc$sdev^2)
all_labels = c("parents",rep("unrelated",N-1),"parents",rep("unrelated",N-1))
plot( pc$rotation[,"PC1"] , pc$rotation[,"PC2"] , type="n" , xlab="PC1" , ylab="PC2" , main="D: Only unrelated" , las=1 )
points(  pc$rotation[all_labels == "unrelated","PC1"] , pc$rotation[all_labels == "unrelated","PC2"] , pch=19 , col=paste(clr[3],"75",sep='') )
points(  pc$rotation[all_labels == "unrelated","PC1"] , pc$rotation[all_labels == "unrelated","PC2"] , pch=1 , col=clr[3] )
points(  pc$rotation[all_labels == "parents","PC1"] , pc$rotation[all_labels == "parents","PC2"] , pch=19 , col=clr[1] )
points(  pc$rotation[all_labels == "parents","PC1"] , pc$rotation[all_labels == "parents","PC2"] , pch=1 , col=clr[1] )
points(  pc$rotation[all_labels == "children","PC1"] , pc$rotation[all_labels == "children","PC2"] , pch=19 , col=paste(clr[2],"75",sep='') )
points(  pc$rotation[all_labels == "children","PC1"] , pc$rotation[all_labels == "children","PC2"] , pch=1 , col=clr[2] )

legend("bottomright",legend=c("Parents","Children","Unrelated") , cex=1 , pch=19 , col=clr)
