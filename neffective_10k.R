# Simulating populations with Ne=10,000 under different growth models
# We will uze `optimize` to search for populations that have equivalent Ne's

# plot side by side
par(mfrow=c(1,2))

# Function for computing the harmonic mean
hmean = function( x ) {
  1/mean(1/x)
}

# Function for generating an exponential population
generate_exp = function( growth , N0 ) {
  gens = 1:2000
  N0 * ( 1 + growth )^(gens-1)
}
test_exp = function( x , N0 = 1000 ) {
  abs( 10e3 - hmean(generate_exp(x,N0)) )
}

gens = 1:2000
r = optimize( test_exp ,c(0,1) , N0 = 1000 )$min
N_all =  generate_exp( r , 1000 )

# Confirm the harmonic mean matches
cat("Exponential population Ne:", hmean(N_all) , '\n' )
plot(gens, N_all ,type="l",las=1,yaxt="n",bty="n",lwd=3,ylab="N",xlab="Generations")
axis(side=2,at=c(10e3,5e6,10e6,20e6),labels=c("10k","5M","10M","20M") ,las=1)
#abline(h=10e3,lty=1,col="royalblue")
lines(c(1,2e3),c(10e3,10e3),col="royalblue")
# Function for generating a bottleneck+exponential population
generate_bneck = function( growth , Ne_bneck = 2500 ) {
  N_all = rep(10e3,2000)
  N_all[ 1001:2000 ] = Ne_bneck * ( 1 + growth )^((1:1000)-1)  
  N_all
}
test_bneck = function( x ) {
  abs( 10e3 - hmean(generate_bneck(x)) )
}

gens = 1:2000
r = optimize( test_bneck ,c(0,1) )$min
N_all = generate_bneck( r )

cat("Bottleneck population Ne:", hmean(N_all) , '\n' )

plot(gens, N_all ,type="l",las=1,yaxt="n",bty="n",lwd=3,ylim=c(0,max(N_all)),ylab="N",xlab="Generations")
axis(side=2,at=c(2.5e3,10e3,50e3,120e3),labels=c("2.5k","10k","50k","120k") ,las=1)
lines(c(1,2e3),c(10e3,10e3),col="royalblue")
