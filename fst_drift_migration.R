# Visualizing equivalent Fst under different drift/migration parameters

library(RColorBrewer)

par(mfrow=c(1,2))

t = seq(1,10e3,by=1)
all.N = c(100,1e3,10e3)
clr = tail(brewer.pal(6,"Blues"),n=length(all.N))
plot(0,0,main="Fst vs. generations of drift",type="n",xlim=range(t),ylim=c(0,1),log="x",las=1,bty="n",xlab="Generations",ylab="Fst_i")

for ( i in 1:length(all.N) ) {
  fst = 1 - (1 - 1/(2*all.N[i]))^t
  lines(t,fst,col=clr[i],lwd=3)
  x = which.min( abs(fst-0.05) )
  points(t[x],fst[x],pch=19,col=clr[i],cex=1.25)
  points(t[x],fst[x],col="white",cex=1.25)
}

legend("topleft",legend=c(all.N),col=clr,bty="n",lwd=3,title="Ne")

m = 1/(10^seq(1,6,length.out=100))
all.N = c(100,1e3,10e3)
clr = tail(brewer.pal(6,"Blues"),n=length(all.N))
plot(0,0,main="Fst vs. migration rate",type="n",xlim=range(m),ylim=c(0,1),las=1,bty="n",xlab="Migration",ylab="Fst",log="x")

for ( i in 1:length(all.N) ) {
  fst = 1/(4*all.N[i]*m + 1)
  lines(m,fst,col=clr[i],lwd=3)
  x = which.min( abs(fst-0.2) )
  points(m[x],fst[x],pch=19,col=clr[i],cex=1.25)
  points(m[x],fst[x],col="white",cex=1.25)
}

legend("bottomleft",legend=c(all.N),col=clr,bty="n",lwd=3,title="Ne")
