# Visualizing the expected age of an allele.

library("RColorBrewer")
clr = brewer.pal(6,"Blues")[c(3,4,5)]

N = 10e3
p = seq(0.001,1,by=0.001)
t = -4*N*p*log(p)/(1-p)
plot_t = 30*t/1e3
plot(p,plot_t,type="l",lwd=4,col="black",las=1,bty="n",ylab="Expected allele age (thousands of years)",xlab="Derived allele frequency")

breaks = c(65e3,200e3,650e3)
for ( i in 1:length(breaks) ) {
  highlight = which.min( abs(30*t - breaks[i]) )
  points( p[highlight],plot_t[highlight],pch=19,col=clr[i],cex=1.25)
  points( p[highlight],plot_t[highlight],col="black",cex=1.25)
  text( p[highlight],plot_t[highlight] , sprintf("%1.2f",p[highlight]) , pos=4,cex=0.75,col=clr[i])
}

legend("bottomright",legend=rev(c("Out of Africa","Modern humans","Hominin divergence")),pch=19,bty="n",col=rev(clr) )