# simple factor analysis example of sampling theory
# ---

# lavaan for factor analysis
library('lavaan')

# --- parameters
# number of people
N = 1e3
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

# 100 different random factors
latent = matrix(runif(N*M),nrow=N,ncol=M)

pdf("g_latent.pdf",height=4,width=5)
pheatmap(latent,cluster_rows=F,cluster_cols=F,treeheight_row=0,treeheight_col=0,border=NA,col = colorRampPalette(c("white", "black"))(50))
dev.off()

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
  tests[,i] = scale(apply(latent[,t],1,mean)) + rnorm(N,0,noise_var)
}

pdf("g_tests.pdf",height=4,width=5)
pheatmap(test_map,xlab="Tests",ylab="Latent Factors",border="black",cluster_rows=T,cluster_cols=T,treeheight_row=0,treeheight_col=0,col = colorRampPalette(c("white", "black"))(50))
dev.off()

# positive correlation matrix
pdf("g_test_cor.pdf",height=4,width=5)
pheatmap(cor(tests)[,rev(1:M_tests)],display_numbers=T,number_format = "%.2f",cluster_rows=F,cluster_cols=F,treeheight_row=0,treeheight_col=0,col = colorRampPalette(c("white", "firebrick3"))(50))
dev.off()

# variance explained
pca = prcomp(tests,center=T,scale=T)
varexp = pca$sdev^2/sum(pca$sdev^2)
names(varexp) = 1:M_tests
pdf("g_var_exp.pdf",height=4,width=5)
barplot(varexp,border=NA,col="black",las=1,xlab="Factor",ylab="Variance Explained")
dev.off()

# CFA
model <- 'G =~ NA*T1+T2+T3+T4+T5+T6+T7+T8
          G ~~ 1*G'
fit <- cfa(model, data=scale(tests), meanstructure=T, missing="fiml")
print(summary(fit, rsquare=T, fit.measures=F))
cat( "CFA metrics:\n" )
print(fitMeasures(fit, c("cfi","tli","rmsea","srmr") ))
cat( "Variance explained:" , varexp , '\n' )


# loadings
library('lavaanPlot')
pl = lavaanPlot(model=fit,coefs = TRUE,edge_options = list(color = "grey"))
embed_plot_pdf(pl, "g_factor_graph.pdf")