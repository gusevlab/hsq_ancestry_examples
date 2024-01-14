# Visualizing the probability of fixation relative to neutral drift

N = 1e6
Ne = 20e3

s = -0.0007
Ne = seq(100,20e3,by=10)
big_s = 4 * Ne * s
relative_p_fix = big_s / (1 - exp(-big_s) )

plot(Ne/1e3,relative_p_fix,type="l",lwd=3,las=1,bty="n",log="x",xlab="Effective population size (Ne)",main="Probability of fixation relative to neutral drift",ylab="Probability",ylim=c(0,1),xaxt="n")
axis(1,at=c(100,500,1e3,5e3,10e3,20e3)/1e3,labels = c("100","500","1,000","5,000","10,000","20,000"))

s = -1e-4
Ne = seq(100,20e3,by=10)
big_s = 4 * Ne * s
relative_p_fix = big_s / (1 - exp(-big_s) )
lines(Ne/1e3,relative_p_fix,lwd=2,col="#4575b4")

s = -1e-5
Ne = seq(100,20e3,by=10)
big_s = 4 * Ne * s
relative_p_fix = big_s / (1 - exp(-big_s) )
lines(Ne/1e3,relative_p_fix,lwd=2,col="#74add1")

legend("bottomleft",legend=c("s = -0.0007","s = -1e-4","s = -1e-5"),col=c("black","#4575b4","#74add1"),lwd=c(3,2,2),bty='n')