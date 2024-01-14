# Sampling and visualizing the bound on Fst estimates

all_p_b = seq(0.01,0.99,by=0.005)
all_p_a = all_p_b

mat = matrix(NA,nrow=length(all_p_a),ncol=length(all_p_b))
mat_hudson = mat
mat_nei86 = mat
mat_nei73 = mat
mat_common = mat

for ( i in 1:length(all_p_b) ) {
  p_b = all_p_b[i]
  
  mat_common[i,] = apply( rbind(apply(rbind(p_b,all_p_a),2,mean) , apply(rbind(1-p_b,1-all_p_a),2,mean)) ,2,max)
  
  # hudson fst
  fst = ((all_p_a - p_b)^2) / ( all_p_a*(1-p_b) + p_b*(1-all_p_a))
  mat_hudson[i,] = fst
  
  # nei fst
  p_avg = (all_p_a + p_b) / 2
  fst = ((all_p_a - p_b)^2) / ( 2*p_avg*(1-p_avg) )
  mat_nei86[i,] = fst
  
  # nei (Edge)
  h_s = apply( rbind( all_p_a^2 + (1-all_p_a)^2 , p_b^2 + (1-p_b)^2 ) , 2 , mean)
  h_t = apply( rbind( apply( rbind( all_p_a , p_b ) , 2 , mean ) , apply( rbind( 1 - all_p_a , 1 - p_b ) , 2 , mean ) )^2 , 2 , sum )
  fst = (h_s - h_t) / (1 - h_t )
  mat_nei73[i,] = fst
}

max_fst_nei73 = rep(NA,length(all_p_b))
max_fst_nei86 = rep(NA,length(all_p_b))
max_fst_hudson = rep(NA,length(all_p_b))

maf_var = sort(unique(c(mat_common)))
for ( i in 1:length(maf_var) ) {
  p_b = maf_var[i]
  max_fst_hudson[i] = max(mat_hudson[mat_common == p_b],na.rm=T)
  max_fst_nei73[i] = max(mat_nei73[mat_common == p_b],na.rm=T)
  max_fst_nei86[i] = max(mat_nei86[mat_common == p_b],na.rm=T)
}

# employ simple smoothing

prev = max_fst_hudson[1]
for ( i in 2:length(maf_var) ) {
  if ( max_fst_hudson[i] < prev - 0.01 ) {
    max_fst_hudson[i] = NA
  } else {
    prev = max_fst_hudson[i]
  }
}

prev = max_fst_nei73[1]
for ( i in 2:length(maf_var) ) {
  if ( max_fst_nei73[i] < prev - 0.01 ) {
    max_fst_nei73[i] = NA
  } else {
    prev = max_fst_nei73[i]
  }
}

prev = max_fst_nei86[1]
for ( i in 2:length(maf_var) ) {
  if ( max_fst_nei86[i] < prev*0.9 | max_fst_nei86[i] < prev - 0.05 ) {
    max_fst_nei86[i] = NA
  } else {
    prev = max_fst_nei86[i]
  }
}

plot( 0 , 0 , xlim=c(0.5,1) , ylim=c(0,1) , type="l" , lwd=3 , las=1 , bty="n" , ylab="Maximum Fst", xlab="Most common allele frequency (either allele)")

lines(maf_var[!is.na(max_fst_hudson)], max_fst_hudson[!is.na(max_fst_hudson)] , lwd=2 )
lines(maf_var[!is.na(max_fst_nei73)], max_fst_nei73[!is.na(max_fst_nei73)] , lwd=2 , lty=2 )
#lines(maf_var[!is.na(max_fst_nei86)], max_fst_nei86[!is.na(max_fst_nei86)] , lwd=2 , lty=3 )

legend("topright",legend=c("Hudson","Nei 1973"),lty=c(1,2),bty="n",title="Estimator",lwd=3)
