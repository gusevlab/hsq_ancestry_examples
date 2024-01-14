# Visualizing the trajectory of an allele under selection

par(mfrow=c(1,2))

# Negative selection
s = -0.0007
gens = 65e3/30
xval = (1:gens)*30/1e3

plot( 0, 0 , type="n" , lwd=3 , xlim=c(0,80) , ylim=c(0,0.5),las=1,bty="n",main="Weak negative selection (s = -0.0007)",xlab="Years (thousands)",ylab="Allele Frequency")

for ( p0 in c(0.05,0.25,0.5) ) {
  p1 = rep(NA,gens)
  p1[1] = p0
  for ( g in 2:gens ) {
    p1[g] = p1[g-1] + s*p1[g-1]*(1-p1[g-1])/(1+2*s*p1[g-1])
  }
  
  lines( xval , p1 ,lwd=3 )
  points( head(xval,n=1), head(p1,n=1) , pch=19 , col="white")
  points( head(xval,n=1), head(p1,n=1) )
  points( tail(xval,n=1), tail(p1,n=1) , pch=19 )
  text(  tail(xval,n=1), tail(p1,n=1) , sprintf("%.2f",tail(p1,n=1)) , pos=4,cex=0.85)
}

# Positive selection
# Example 5.1 with positive selection

s = 0.0007
gens = 65e3/30
xval = (1:gens)*30/1e3

plot( 0, 0 , type="n" , lwd=3 , xlim=c(0,80) , ylim=c(0,0.001),las=1,bty="n",main="Weak positive selection (s = +0.0007)",xlab="Years (thousands)",ylab="Allele Frequency",yaxt="n")
axis(2,at=c(0,0.001),las=1)
for ( p0 in c(1/10e3) ) {
  p1 = rep(NA,gens)
  p1[1] = p0
  for ( g in 2:gens ) {
    p1[g] = p1[g-1] + s*p1[g-1]*(1-p1[g-1])/(1+2*s*p1[g-1])
  }
  
  lines( xval , p1 ,lwd=3 )
  points( head(xval,n=1), head(p1,n=1) , pch=19 , col="white")
  points( head(xval,n=1), head(p1,n=1) )
  points( tail(xval,n=1), tail(p1,n=1) , pch=19 )
  text(  tail(xval,n=1), tail(p1,n=1) , sprintf("%.5f",tail(p1,n=1)) , pos=4,cex=0.85)
}
