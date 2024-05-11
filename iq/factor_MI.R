# simple factor analysis example
# ---

# lavaan for factor analysis
library('lavaan')
library('lavaanPlot')

# --- parameters
# number of people
N = 5e3
# two groups
group = rbinom(N,1,0.5)
# number of factors
M = 100
# number of tests
M_tests = 8
# fraction of factors going into each test
p = 0.5
# measurement error
# error SD for each test will be sampled from between these values
min_err = 0.5
max_err = 1.5
# --- parameters
seeds = 50

all_MODES = 1:5
all_ndiffs = 30

results_gap = matrix( NA , nrow=length(all_MODES) , ncol=length(all_ndiffs) )
results_weak = matrix( NA , nrow=length(all_MODES) , ncol=length(all_ndiffs) )
results_strong = matrix( NA , nrow=length(all_MODES) , ncol=length(all_ndiffs) )

for ( m in 1:5 ) {
  MODE = all_MODES[m]
  for ( n in 1:length(all_ndiffs) ) {
    n_diff = all_ndiffs[n]
    pv_weak = rep(NA,seeds)
    pv_strong = rep(NA,seeds)
    gap = rep(NA,seeds)
    for ( s in 1:seeds ) {
      # 100 different random factors
      latent = matrix(runif(N*M),nrow=N,ncol=M)
      # group differences in some factors
      if ( n_diff != 0 ) {
        if ( MODE == 1 ) {
          # for each individual, shrink *the same* abilities
          latent[group==0,1:n_diff] = 0.75 * latent[group==0,1:n_diff]
        } else if ( MODE == 3 ) {
          for ( ind in which(group == 0) ) {
            # for each individual, shrink a fraction of random abilities
            keep = sample(M,n_diff)
            latent[ind,keep] = 0.75 * latent[ind,keep]
          }
        } else if ( MODE == 2 ) {
          for ( ind in which(group == 0) ) {
            # for each individual, lower the same abilities
            keep = 1:n_diff
            newval = latent[ind,keep] - rnorm(length(keep),0.15,0.1)
            newval[ newval < 0 ] = 0
            newval[ newval > 1 ] = 1
            latent[ind,keep] = newval
          }
        } else if ( MODE == 4 ) {
          for ( ind in which(group == 0) ) {
            # for each individual, lower a random fraction
            keep = sample(M,n_diff)
            newval = latent[ind,keep] - rnorm(length(keep),0.15,0.1)
            newval[ newval < 0 ] = 0
            newval[ newval > 1 ] = 1
            latent[ind,keep] = newval
          }
        }
      }
      
      # 8 tests that sample from each of them
      tests = matrix(NA,nrow=N,ncol=M_tests)
      colnames(tests) = paste("T",1:ncol(tests),sep='')
      shared = rnorm(N,0,0.1)
      test_map = matrix(0,nrow=M,ncol=M_tests)
      for ( i in 1:M_tests ) {
        # sample factors going into this test
        t = sample(1:M,M*p)
        # record the mapping
        test_map[t,i] = 1
        # sample some measurement error for this test
        noise_var = runif(1,min_err,max_err)
        # lower the performance on all tests
        if ( MODE == 5 ) {
          tests[,i] = scale(apply(latent[,t],1,mean)) + rnorm(N,0,noise_var) + rnorm(N,group*(n_diff/30),1)
        } else {
          tests[,i] = scale(apply(latent[,t],1,mean)) + rnorm(N,0,noise_var)
        }
      }
  
      tests = scale(tests)
      
      # variance explained
      #pca = prcomp(tests,center=T,scale=T)
      #varexp = pca$sdev^2/sum(pca$sdev^2)
      #names(varexp) = 1:M_tests
      #barplot(varexp,border=NA,col="black",las=1,xlab="Factor",ylab="Variance Explained")
      #cat( "Variance explained:" , varexp , '\n' )
      
      try({
      # CFA
      model <- 'G =~ NA*T1+T2+T3+T4+T5+T6+T7+T8
                G ~~ 1*G'
      fit <- cfa(model, data=tests, meanstructure=T, missing="fiml")
      # look at group diffs
      g = scale(predict( fit , tests ))
      gap[s] = mean(g[group==0]) - mean(g[group==1])
      
      # MI testing
      tests = cbind(tests,group)
      fit1 <- cfa(model, data=tests, meanstructure=T, missing="fiml",group="group")
      # weak invariance
      fit2 <- cfa(model, data=tests, meanstructure=T, missing="fiml",group="group",group.equal = "loadings")
      # strong invariance
      fit3 <- cfa(model, data=tests, meanstructure=T, missing="fiml",group="group",group.equal = c("intercepts", "loadings"))
      # model comparison tests
      
      tst = lavTestLRT(fit1, fit2)
      pv_weak[s] = tst$"Pr(>Chisq)"[2]
      tst = lavTestLRT(fit1, fit3)
      pv_strong[s] = tst$"Pr(>Chisq)"[2]
      },silent=T)
      #cat( s , mean(pv_weak<0.05,na.rm=T) , mean(pv_strong<0.05,na.rm=T) , '\n' )
    }
    results_gap[m,n] = mean(gap,na.rm=T)
    results_weak[m,n] = mean(pv_weak<0.05,na.rm=T)
    results_strong[m,n] = mean(pv_strong<0.05,na.rm=T)
    cat( "***" , n_diff , mean(gap,na.rm=T) , mean(pv_weak<0.05,na.rm=T) , mean(pv_strong<0.05,na.rm=T) ,'\n' )
  }
}


library("RColorBrewer")
mode_names = c("Subtest Shrink","Subtest Difference","Random Shrink","Random Difference","Mean Shift")
clr = brewer.pal(4,"Blues")[c(2:4)]
par(mar = c(5, 10, 2, 2))
barplot( t(cbind(results_weak,results_strong))*100 , border=clr[3] , xlab="Power (%)" , col=clr[1:2] , las=1 , names=mode_names , horiz=T , beside=T)
legend("topright",legend=c("Loadings (Weak)","Loadings + Intercepts (Strong)"),fill=clr,bty="n")

# par(mfrow=c(1,length(all_MODES)))
# clr = brewer.pal(length(all_MODES),"Set1")
# for ( MODE in all_MODES ) {
#   plot( abs(results_gap[MODE,]) , results_strong[MODE,] , xlab="Group difference" , ylab="Power", type="n" , ylim=c(0,1) , las=1 , xlim=c(0,2) , main=mode_names[MODE])
#   lines( abs(results_gap[MODE,]) , results_strong[MODE,] , lty=1 , col=clr[MODE])
#   points( abs(results_gap[MODE,]) , results_strong[MODE,] , pch=19,lty=1 , col=clr[MODE])
#   lines( abs(results_gap[MODE,]) , results_weak[MODE,] , lty=2 , col=clr[MODE])
#   points( abs(results_gap[MODE,]) , results_weak[MODE,] , pch=19, lty=2 , col=clr[MODE])
# }
