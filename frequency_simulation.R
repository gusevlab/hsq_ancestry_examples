# Visualizing the drift of alleles with/without selection.

library(RColorBrewer)

gens = as.integer(65e3 / 30)

par(mfrow=c(2,2))
params_Ne = c(500,10e3,500,10e3)
params_s = c(0,0,0.005,0.005)
params_main = c("a. High Drift + No Selection","b. Low Drift + No Selection","c. High Drift + Selection","d. Low Drift + Selection")

cols <- brewer.pal(3, "YlGnBu")
cols = c("#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58")
pal <- colorRampPalette(cols)(gens)

for ( types in 1:4 ) {
  Ne = params_Ne[types]
  s = params_s[types]
  plot( 0 , 0 , type="n" , bty="n", xlim=c(0,gens) , ylim=c(0,1) , xlab="Generations", ylab="Frequency", las=1 , main=params_main[types])
  
  # sample a mutation
  for ( i in 1:500 ) {
    g_start = as.integer(runif(1,1,gens))
    freq = rep(NA,gens)
    for ( g in g_start:gens ) {
      if ( g == g_start ) {
        freq[g] = 1/Ne
      } else {
        cur_f = freq[g-1] + s*freq[g-1]*(1-freq[g-1])/(1+2*s*freq[g-1])
        freq[g] = mean(rbinom(Ne,2,cur_f)/2)
      }
      if ( freq[g] == 0 ) {
        break()
      } else if ( freq[g] == 1 ) {
        points( g , freq[g], col=pal[g_start] , pch=19 , cex=0.75 )
        break()
      }
    }
    lines( 1:gens , freq , col=pal[g_start])
    #points( g_start , freq[g_start], col=pal[g_start] , pch=19 , cex=0.5 )
  }
}