set.seed(42)

source("reml_func.R")

# --- parameters
# Num individuals in the population
N = 5e3
# Num individuals to use for estimation
N_est = 1e3
# Num markers
M = 100
# Num causal
M_causal = 50
# minor allele frequency
maf = 0.5
# variance in child trait explained by genetic variation (i.e. direct heritability)
var_direct = 0.2
# variance in child trait explained by indirect effects
var_indirect = 0.1
# sibling indirect effects (typically set to zero)
var_sib = 0
# number of iterations to run
seeds = 5
# Put "HE" or "REML" here for algorithm to use
METHOD = "REML"
# Spousal correlation
AM_cor = 0
# Assortative Mating on the underlying genetic value
AM_GENO = FALSE
# --- done parameters

# matrices to store the estimates and standard errors
est_h2 = matrix(NA,nrow=seeds,ncol=4)
est_se = matrix(NA,nrow=seeds,ncol=4)
colnames(est_h2) = c("True","POP_h2","RDR_fam","RDR_sib")
colnames(est_se) = c("True","POP_h2","RDR_fam","RDR_sib")

for ( s in 1:seeds ) {
  # ---
  # sample effect sizes for direct effect (follow Young et al. and just sample from std norm)
  betas = rnorm(M,0,1)
  # set only some variants to be causal
  betas[ sample(M,M-M_causal) ] = 0
  
  # sample 0/1 haplotypes
  mate1_hap_a = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  mate1_hap_b = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  
  # sample a random mate
  mate2_hap_a = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  mate2_hap_b = matrix( rbinom(N*M,1,maf) , nrow=N , ncol=M )
  
  ## --- Assortative Mating until equilibrium
  if ( AM_cor != 0 ) {
  for ( gens in 1:10 ) {
    
    # compute genetic value for direct phenotype
    gv_mate1_direct = scale( mate1_hap_a + mate1_hap_b ) %*% betas
    gv_mate2_direct = scale( mate2_hap_a + mate2_hap_b ) %*% betas

    pheno_mate1 = sqrt(var_direct) * scale(gv_mate1_direct) + rnorm(N,0,sqrt(1-var_direct))
    pheno_mate2 = sqrt(var_direct) * scale(gv_mate2_direct) + rnorm(N,0,sqrt(1-var_direct))
    
    # rank the phenotypes and find similar ranks
    if ( AM_cor != 0 ) {
      if ( AM_GENO ) {
        ord1 = order(gv_mate1_direct)
        ord2 = order( AM_cor * scale(gv_mate2_direct) + rnorm(N,0,sqrt(1-AM_cor^2)))
      } else {
        ord1 = order(pheno_mate1)
        ord2 = order( AM_cor * scale(pheno_mate2) + rnorm(N,0,sqrt(1-AM_cor^2)))
      }
    } else {
      ord1 = sample(1:N)
      ord2 = sample(1:N)
    }
    
    mate1_hap_a = mate1_hap_a[ord1,]
    mate1_hap_b = mate1_hap_b[ord1,]
    mate2_hap_a = mate2_hap_a[ord2,]
    mate2_hap_b = mate2_hap_b[ord2,]
    pheno_mate1 = pheno_mate1[ord1]
    pheno_mate2 = pheno_mate2[ord2]
    gv_mate1_direct = gv_mate1_direct[ord1]
    gv_mate2_direct = gv_mate2_direct[ord2]
    
    # --- mate and generate kids
    
    child1_hap_a = mate1_hap_a
    sampled_vars = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
    child1_hap_a[ !sampled_vars ] = mate1_hap_b[ !sampled_vars ]
    
    child1_hap_b = mate2_hap_a
    sampled_vars = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
    child1_hap_b[ !sampled_vars ] = mate2_hap_b[ !sampled_vars ]
    
    # --- find new mate and generate kids
    if ( AM_cor != 0 ) {
      if ( AM_GENO ) {
        ord1 = order(gv_mate1_direct)
        ord2 = order( AM_cor * scale(gv_mate2_direct) + rnorm(N,0,sqrt(1-AM_cor^2)))
      } else {
        ord1 = order(pheno_mate1)
        ord2 = order( AM_cor * scale(pheno_mate2) + rnorm(N,0,sqrt(1-AM_cor^2)))
      }
    } else {
      ord1 = sample(1:N)
      ord2 = sample(1:N)
    }
    
    mate1_hap_a = mate1_hap_a[ord1,]
    mate1_hap_b = mate1_hap_b[ord1,]
    mate2_hap_a = mate2_hap_a[ord2,]
    mate2_hap_b = mate2_hap_b[ord2,]
    pheno_mate1 = pheno_mate1[ord1]
    pheno_mate2 = pheno_mate2[ord2]
    gv_mate1_direct = gv_mate1_direct[ord1]
    gv_mate2_direct = gv_mate1_direct[ord2]
    
    child2_hap_a = mate1_hap_a
    sampled_vars = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
    child2_hap_a[ !sampled_vars ] = mate1_hap_b[ !sampled_vars ]
    
    child2_hap_b = mate2_hap_a
    sampled_vars = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
    child2_hap_b[ !sampled_vars ] = mate2_hap_b[ !sampled_vars ]
    
    # update to the next generation
    mate1_hap_a = child1_hap_a
    mate1_hap_b = child1_hap_b
    mate2_hap_a = child2_hap_a
    mate2_hap_b = child2_hap_b
  }
  }
  ### ----

  # --- generate one child per family
  child_hap_a = mate1_hap_a
  sampled_vars_mate1 = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
  child_hap_a[ !sampled_vars_mate1 ] = mate1_hap_b[ !sampled_vars_mate1 ]
  
  child_hap_b = mate2_hap_a
  sampled_vars_mate2 = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
  child_hap_b[ !sampled_vars_mate2 ] = mate2_hap_b[ !sampled_vars_mate2 ]
  
  # --- generate a sibling
  sib_hap_a = mate1_hap_a
  sampled_vars_mate1 = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
  sib_hap_a[ !sampled_vars_mate1 ] = mate1_hap_b[ !sampled_vars_mate1 ]
  
  sib_hap_b = mate2_hap_a
  sampled_vars_mate2 = matrix( as.logical(rbinom(N*M,1,0.5)) , nrow=N , ncol=M )
  sib_hap_b[ !sampled_vars_mate2 ] = mate2_hap_b[ !sampled_vars_mate2 ]
  
  # compute phenotype for sib
  gv_sib_direct = scale( sib_hap_a + sib_hap_b ) %*% betas
  gv_sib_indirect = scale( mate1_hap_a + mate1_hap_b + mate2_hap_a + mate2_hap_b ) %*% betas
  pheno_sib = sqrt(var_direct) * scale(gv_sib_direct) + sqrt(var_indirect) * scale(gv_sib_indirect)
  # add noise to make phenotype have variance 1
  pheno_sib = pheno_sib + rnorm(N,0,sqrt(1-var(pheno_sib)))
  
  # compute phenotype for child
  gv_child_direct = scale( child_hap_a + child_hap_b ) %*% betas
  gv_child_indirect = scale( mate1_hap_a + mate1_hap_b + mate2_hap_a + mate2_hap_b ) %*% betas
  pheno_child = sqrt(var_direct) * scale(gv_child_direct) + sqrt(var_indirect) * scale(gv_child_indirect) + sqrt(var_sib) * scale(pheno_sib)
  # add noise to make phenotype have variance 1
  pheno_child = pheno_child + rnorm(N,0,sqrt(1-var(pheno_child)))
  
  # --- Estimate:
  
  # truth
  reg = summary(lm( pheno_child ~ scale(gv_child_direct) + scale(gv_child_indirect) ))
  est_h2[s,"True"] = reg$coef[2,1]^2
  
  # --- estimate heritability
  keep = sample( N , N_est )
  
  # kinship for the kids
  x_child = scale( child_hap_a + child_hap_b )[keep,]
  k_child = x_child %*% t(x_child) / ncol(x_child)
  K = list()
  K[[1]] = k_child

  # --- estimate total/popoulation heritability (just use HE for speed)
  #if ( METHOD == "REML" ) {
  #  reml = aiML( A = K , y = scale(pheno_child[keep]) , Var = c(0.5,0.5) , verbose=FALSE )
  #  est_h2[s,"POP_h2"] = reml$h2
  #  est_se[s,"POP_h2"] = reml$se
  #} else {
    reml = HEreg( A = K , y = scale(pheno_child[keep]) )
    est_h2[s,"POP_h2"] = reml$h2
  #}
  
  # kinship for the parents
  x_parents = sqrt(2) * scale( mate1_hap_a + mate1_hap_b + mate2_hap_a + mate2_hap_b )[keep,]
  k_par = x_parents %*% t(x_parents) / (2*ncol(x_parents))
  k_opar = ( x_child %*% t(x_parents) +  x_parents %*% t(x_child) ) / (2*ncol(x_parents))
  K = list()
  K[[1]] = k_child
  K[[2]] = k_par
  K[[3]] = k_opar
  # --- estimate classic RDR
  if ( METHOD == "REML" ) {
    reml = aiML( A = K , y = scale(pheno_child[keep]) , Var = c(0.25,0.25,0.25,0.25) , verbose=FALSE )
    est_h2[s,"RDR_fam"] = reml$h2
    est_se[s,"RDR_fam"] = reml$se
  } else {
    reml = HEreg( A = K , y = scale(pheno_child[keep]) )
    est_h2[s,"RDR_fam"] = reml$h2
  }
  
  # RDR from siblings
  # kinship for the kids (already generated)
  x_child = scale( child_hap_a + child_hap_b )[keep,]
  k_child = x_child %*% t(x_child) / ncol(x_child)
  # kinship for the parents using the sibling mean
  x_parents = sqrt(2) * scale( child_hap_a + child_hap_b + sib_hap_a + sib_hap_b )[keep,]
  k_par = x_parents %*% t(x_parents) / (2*ncol(x_parents))
  k_opar = ( x_child %*% t(x_parents) +  x_parents %*% t(x_child) ) / (2*ncol(x_parents))
  # estimate
  K = list()
  K[[1]] = k_child
  K[[2]] = k_par
  K[[3]] = k_opar
  if ( METHOD == "REML") {
    reml = aiML( A = K , y = scale(pheno_child[keep]) , Var = c(0.25,0.25,0.25,0.25) , verbose=FALSE )
    est_h2[s,"RDR_sib"] = reml$h2
    est_se[s,"RDR_sib"] = reml$se
  } else {
    reml = HEreg( A = K , y = scale(pheno_child[keep]) )
    est_h2[s,"RDR_sib"] = reml$h2
  }
  cat( AM_cor , s , round(apply(est_h2,2,mean,na.rm=T),3) , " s.e. " , round(apply( est_h2 , 2 , sd , na.rm=T ) / sqrt(apply( !is.na(est_h2) , 2 , sum , na.rm=T )),3) , '\n' , sep='\t' )
  
}

write.table( cbind(apply(est_h2,2,mean,na.rm=T) , 
                   apply(est_se,2,mean,na.rm=T) , 
                   apply( est_h2 , 2 , sd , na.rm=T ) / sqrt(apply( !is.na(est_h2) , 2 , sum , na.rm=T )) ) , quote=F , sep='\t' , row.names=T , col.names = c("Mean_h2","Mean_SE","SE_h2") )
