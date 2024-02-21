library("pROC")

library("RColorBrewer")

M = 500
p0 = rep(0.5,M)
#p0 = runif(M,0.01,0.99)
Ne = 10e3

all_gens = c(2.5e3,3e3,3.7e3,300e3)
fst_est = rep(NA,length(all_gens))
fst_seeds = 10
for ( g in 1:length(all_gens) ) {
  gens = all_gens[g]
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  maf_anc = p0
  all_fst = rep(NA,fst_seeds)
  for ( s in 1:fst_seeds ) {
    maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift_var) )
    maf_pop2 = maf_anc + rnorm(M,0,sqrt(drift_var) )
    maf_pop2[ maf_pop2 > 1 ] = 1
    maf_pop2[ maf_pop2 < 0 ] = 0
    maf_pop1[ maf_pop1 > 1 ] = 1
    maf_pop1[ maf_pop1 < 0 ] = 0
    # hudson fst
    fst = ((maf_pop1 - maf_pop2)^2) / ( maf_pop1*(1-maf_pop2) + maf_pop2*(1-maf_pop1))
    # neu fst
    #h_s = apply( rbind( maf_pop1^2 + (1-maf_pop1)^2 , maf_pop2^2 + (1-maf_pop2)^2 ) , 2 , mean)
    #h_t = apply( rbind( apply( rbind( maf_pop1 , maf_pop2 ) , 2 , mean ) , apply( rbind( 1 - maf_pop1 , 1 - maf_pop2 ) , 2 , mean ) )^2 , 2 , sum )
    #fst = (h_s - h_t) / (1 - h_t )
    
    # ---
    fst = mean(fst,na.rm=T)
    all_fst[s] = fst
  }
  fst_est[g] = mean(all_fst)
}
fst_est = sprintf("%1.2f",fst_est)
cat( fst_est , '\n' , sep='\t')

df = data.frame( "Fst" = vector() , "AUC" = vector() , "Method" = vector() )
for ( g in 1:length(all_gens) ) {
  gens = all_gens[g]
  drift_var = p0 * (1-p0) * ( gens / (2*Ne) )
  maf_anc = p0
  
  pop_lbl = c(rep(1,N),rep(0,N))
  maf_pop1 = maf_anc + rnorm(M,0,sqrt(drift_var) )
  maf_pop2 = maf_anc + rnorm(M,0,sqrt(drift_var) )
  
  maf_pop2[ maf_pop2 > 1 ] = 1
  maf_pop2[ maf_pop2 < 0 ] = 0
  
  maf_pop1[ maf_pop1 > 1 ] = 1
  maf_pop1[ maf_pop1 < 0 ] = 0
  
  x = matrix(NA,nrow=N*2,ncol=M)
  
  all_roc = rep(NA,M)
  rand_roc = rep(NA,M)
  for ( snp in 1:M ) {
    x[1:N,snp] = rbinom(N,2,maf_pop1[snp])
    x[(N+1):(N+N),snp] = rbinom(N,2,maf_pop2[snp])
    
    if ( mean(x[,snp])/2 > 0.5 ) {
      y = x[,snp] == 0
    } else {
      y = x[,snp] == 2
    }
    if ( sd(y) != 0 ) {
      cur_roc = roc(response = y , predictor = pop_lbl ,ci=T,quiet=T)
      all_roc[snp] = cur_roc$auc
      
      cur_roc = roc(response = y , predictor = sample(pop_lbl) ,ci=T,quiet=T)
      rand_roc[snp] = cur_roc$auc
    }
  }
  
  # compute AUC
  cat( fst_est[g] , mean(all_roc,na.rm=T) , mean(rand_roc,na.rm=T) , '\n' )
  df = rbind( df , data.frame( "Fst" = fst_est[g] , "AUC" = all_roc , "Method" = "Ancestry" ) )
  df = rbind( df , data.frame( "Fst" = fst_est[g] , "AUC" = rand_roc , "Method" = "Random" ) )
  
  
  # now try with families
  mate = sample(1:N)
  x_mum = x[1:N,]
  x_dad = x[sample(1:N),]
  x_kid = x_mum
  x_kid2 = x_mum
  all_roc_parents = rep(NA,M)
  all_roc_sibling = rep(NA,M)
  for ( snp in 1:M ) {
    x_kid[,snp] = rbinom(N,1, (x_mum[,snp])/2 ) + rbinom(N,2, (x_dad[,snp])/2 )
    x_kid2[,snp] = rbinom(N,1, (x_mum[,snp])/2 ) + rbinom(N,2, (x_dad[,snp])/2 )

    if ( mean(x_kid[,snp])/2 > 0.5 ) {
      y = x_kid[,snp] == 0
    } else {
      y = x_kid[,snp] == 2
    }

    if ( sd(y) != 0 ) {
      # train a simple predictor to include interactions
      #reg = lm(y ~ x_mum[,snp]+x_dad[,snp]+x_mum[,snp]*x_dad[,snp])
      #pred = fitted(reg)
      pred = x_mum[,snp]+x_dad[,snp]
      cur_roc = roc(response = y , predictor = pred ,ci=F,quiet=T)
      all_roc_parents[snp] = cur_roc$auc
      
      #cur_roc = roc(response = y , predictor = x_mum[,snp]+x_dad[,snp] ,ci=F,quiet=T)
      #all_roc_parents[snp] = cur_roc$auc
      cur_roc = roc(response = y , predictor = x_kid2[,snp] ,ci=F,quiet=T)
      all_roc_sibling[snp] = cur_roc$auc
    }
  }
  df = rbind( df , data.frame( "Fst" = fst_est[g] , "AUC" = all_roc_parents , "Method" = "Family: Parents" ) )
  df = rbind( df , data.frame( "Fst" = fst_est[g] , "AUC" = all_roc_sibling , "Method" = "Family: Sibling" ) )
  
}


# --- plot the model results
library("ggplot2")
# remove outliers
df$Method = factor(df$Method, levels=c("Family: Parents","Family: Sibling","Ancestry","Random"))
df$Fst = as.factor(df$Fst)
df = df[ !is.na(df$AUC) & df$AUC > 0.4 , ]
ggplot(df, aes(x=Fst, y=AUC , fill=Method )) + 
  theme_bw() + stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
                                 geom="pointrange", aes(color=Method),
                                 shape = 18, size = 0.75,
                                 position = position_dodge(width = 0.5)) + 
  ylim(0.5,1) + scale_colour_manual(values = c("#d53e4f","#fc8d59","#3288bd","black")) + theme(legend.position="top")
