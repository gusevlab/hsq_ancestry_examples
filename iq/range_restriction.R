# Range restriction
# See examples in Sackett + Yang : https://doi.org/10.1037/0021-9010.85.1.112
N = 5e3
x = rnorm(N,0,1)
y = x*0.5 + rnorm(N,0,0.5)

par(mfrow=c(2,5))
lims = c(-4,4)
plot(x,y,pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="True relationship",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = x > 1
plot(x[keep],y[keep],pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on x",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(y[keep]~x[keep])
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = x > 2 | x < -2
plot(x[keep],y[keep],pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on tails of x",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(y[keep]~x[keep])
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = x < 1 & x > -1
plot(x[keep],y[keep],pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on middle of x",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(y[keep]~x[keep])
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = y > 1
plot(x[keep],y[keep],pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on y",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(y[keep]~x[keep])
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

# --- scaled

x = scale(x)
y = scale(y)
plot(x,y,pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="True relationship (scaled)",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = x > 1
plot(scale(x[keep]),scale(y[keep]),pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on x (scaled)",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(scale(y[keep])~scale(x[keep]) )
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = x > 2 | x < -2
plot(scale(x[keep]),scale(y[keep]),pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on tails of x (scaled)",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(scale(y[keep])~scale(x[keep]) )
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = x < 1 & x > -1
plot(scale(x[keep]),scale(y[keep]),pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on middle of x (scaled)",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(scale(y[keep])~scale(x[keep]) )
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")

keep = y > 1
plot(scale(x[keep]),scale(y[keep]),pch=19,cex=0.5,col="#00000025",xlim=lims,ylim=lims,las=1,main="Restriction on y (scaled)",xlab="x",ylab="y")
reg = lm(y~x)
abline(reg,col="royalblue",lty=1,lwd=2)
reg = lm(scale(y[keep])~scale(x[keep]) )
abline(reg,col="royalblue",lty=2,lwd=2)
legend("bottomright",legend=paste("b=",sprintf("%2.1f",reg$coef[2]),sep=''),bty="n",text.col="royalblue")
