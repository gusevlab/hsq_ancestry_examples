# simple sampling theory analysis example for factor replication
# ---

# lavaan for factor analysis
library('lavaan')

# --- parameters
# number of people
N = 1e3
# number of processes
M = 100
# number of tests
M_tests = 50
# fraction of processes going into each test
p = 0.25
# measurement error
# error SD for each test will be sampled from between these values
min_err = 0
max_err = 0.5
# --- parameters

# 100 different random processes
latent = matrix(runif(N*M),nrow=N,ncol=M)

tests = matrix(NA,nrow=N,ncol=M_tests)
test_map = matrix(0,nrow=M,ncol=M_tests)
for ( i in 1:M_tests ) {
  # sample processes going into this test
  t = sample(1:M,M*p)
  # sample some measurement error for this test
  noise_var = runif(1,min_err,max_err)
  tests[,i] = scale(apply(latent[,t],1,mean)) + rnorm(N,0,noise_var)
}

# variance explained
pca = prcomp(t(scale(tests)),center=F,scale=F)
varexp = pca$sdev^2/sum(pca$sdev^2)
cat(varexp[1],'\n')
g_1 = pca$rotation[,1]

# make a new test battery
tests = matrix(NA,nrow=N,ncol=M_tests)
for ( i in 1:M_tests ) {
  # sample processes going into this test
  t = sample(1:M,M*p)
  # sample some measurement error for this test
  noise_var = runif(1,min_err,max_err)
  tests[,i] = scale(apply(latent[,t],1,mean)) + rnorm(N,0,noise_var)
}

# variance explained
pca = prcomp(t(scale(tests)),center=F,scale=F)
varexp = pca$sdev^2/sum(pca$sdev^2)
cat(varexp[1],'\n')

g_2 = pca$rotation[,1]
cat( cor(g_1,g_2) , '\n' )
