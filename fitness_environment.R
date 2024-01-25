g = rep(NA,3)
e = rep(NA,3)
f = rep(NA,3)

g[1] = 0
e[1] = 0
f[1] = 0

g[2] = 5
e[2] = 0
f[2] = 2.5

g[3] = 5
e[3] = -5
f[3] = 0

title = c("Baseline","Shift in optimum","Shift in environment")
par(mfrow=c(1,3))

for ( i in 1:3 ) {
xlm = c(-10,10)
ylm = c(-0.75,0.75)
plot( 0 , 0 , type="n" , xlim=xlm , ylim=ylm , bty="n" , las=1 , yaxt="n" , ylab="" , xlab="" , main=title[i])
x = seq(g[i]-6, g[i]+6, by=0.1) 
x = x[x>-10 & x<10]
y = dnorm(x, g[i], 1)
polygon( x , y , col="#4575b475", border="#4575b475" )

x = seq(e[i]-6, e[i]+6, by=0.1) 
x = x[x>-10 & x<10]
y = dnorm(x, e[i], 1)
polygon( x , -1*y , col="#1a985075" , border="#1a985075")
points(f[i],0,pch=19,cex=1.5)
abline(v=0,lty=3)
lines(xlm,c(0,0))
if ( i==1) {
  legend(legend=c("Genetics","Environment","Optimum"),"topleft",col=c("#4575b475","#1a985075","black"),bty="n",pch=c(NA,NA,19),lwd=c(3,3,NA))
}
}