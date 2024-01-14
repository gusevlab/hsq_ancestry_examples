# Visualizing Nei's and Hudson's Fst

par(mfrow=c(1,2))
for ( p_b in c(0.1,0.5)) {
  p_a = seq(0.01,0.99,by=0.01)
  fst = list()
  
  # Hudson
  fst[[1]] = 1 - (p_a*(1-p_a)+p_b*(1-p_b))/(p_a*(1-p_b)+p_b*(1-p_a))

  # Nei:
  p_t = (p_a + p_b) / 2
  p_t * (1-p_t)
  h_s_bar = (p_b * (1-p_b) + p_a*(1-p_a))/2
  fst[[2]] = 1 - (h_s_bar) / ( p_t * (1-p_t))

  plot(0,0,type="n",ylim=c(0,1),xlim=c(0,1),main=paste("p1 = ",sprintf("%1.2f",p_b),sep=''),ylab="Fst",xlab="p2",las=1,bty="n")
  clr = brewer.pal(4,"Set2")[c(1,3)]
  
  for ( i in 1:2 ) {
    lines(p_a,fst[[i]],col=clr[i],lwd=2)
  }
  points(p_b,0,pch=19,col="black")
  points(p_b,0,col="white")
  legend("topleft",legend=c("Hudson","Nei 1973"),col=clr,lty=1,bty="n",lwd=2)
}
