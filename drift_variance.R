# Visualizing variance due to drift and sampling

total_gens = as.integer(2*65e3/30)
gens = (1:total_gens)
xval = gens*30/1e3
Ne = 10e3

plot( 0, 0 , type="n" , lwd=3 , xlim=c(0,2*65) , ylim=c(0.02,0.08),las=1,bty="n",xlab="Years (thousands)",ylab="Allele Frequency",main="Variance due to neutral drift and sampling")

for( p0 in c(0.05) ) {
  
  samp_0 = 10
  samp_t = 100
  drift_var = p0 * (1-p0) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  polygon( c(xval,rev(xval)) , c(p0 - (drift_var),rev(p0 + (drift_var))) , col=NA , border="#08519c" , lty=2 )
  
  samp_0 = 100
  samp_t = 100
  drift_var = p0 * (1-p0) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  polygon( c(xval,rev(xval)) , c(p0 - (drift_var),rev(p0 + (drift_var))) , col=NA , border="#08519c" , lty=3 )
  
  
  samp_0 = 10
  samp_t = 10
  drift_var = p0 * (1-p0) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  polygon( c(xval,rev(xval)) , c(p0 - (drift_var),rev(p0 + (drift_var))) , col=NA , border="#08519c" )
  
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  polygon( c(xval,rev(xval)) , c(p0 - (drift_var),rev(p0 + (drift_var))) , col="#08519c50" , border=NA )
  
  lines( c(0,tail(xval,n=1)) , c(p0,p0) ,col="#08519c" , lwd=2 , lty=1 )
}

legend("topleft",legend=c("drift alone",
                          expression(paste('drift, n'[0], ' = ', 10,', n'[t], ' = ', 10)),
                          expression(paste('drift, n'[0], ' = ', 10,', n'[t], ' = ', 100)),
                          expression(paste('drift, n'[0], ' = ', 100,', n'[t], ' = ', 100))),
       bty="n",
       lty=c(1,1,2,3),lwd=c(3,1,1,1),col=c("#08519c50","#08519c","#08519c","#08519c"))