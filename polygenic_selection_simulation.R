# Visualizing the drift of alleles with/without selection.

gens = 2000
Ne = 10e3
s = -0.0007
M = 5e3

df = data.frame( "Type" = vector() , "Generation" = vector() , "Frequency" = vector() )

for ( stype in c("Neutral","Directional","Stabilizing")) {
freq_mat = matrix( NA , nrow=gens , ncol=M )
freq_mat[1,] = runif(M,0.01,0.99)

for ( g in 2:gens ) {
  keep = !is.na(freq_mat[g-1,]) & freq_mat[g-1,] > 0 & freq_mat[g-1,] < 1
  cur_f = freq_mat[g-1,keep]
  if (stype == "Directional" ) {
    cur_f = cur_f + s*cur_f*(1-cur_f)/(1+2*s*cur_f)
  } else if (stype == "Stabilizing" ) {
    cur_f = cur_f - s*cur_f*(1-cur_f)*(cur_f-0.5)
  }
  cur_f[ cur_f < 0 ] = 0
  cur_f[ cur_f > 1 ] = 1
  
  # add drift variance
  cur_f = cur_f + rnorm(length(cur_f),0,sqrt( (cur_f*(1-cur_f))/(2*Ne) ))
  cur_f[ cur_f < 0 ] = 0
  cur_f[ cur_f > 1 ] = 1
  # fill in the rest for fixed alleles
  freq_mat[g:gens, which(keep)[(cur_f == 0)] ] = 0
  freq_mat[g:gens, which(keep)[(cur_f == 1)] ] = 1
  freq_mat[g,keep] = cur_f
}

# for plotting
df = rbind(df , data.frame( "Type" = stype , "Generation" = rep(1,M) , "Frequency" = freq_mat[1,]))
for ( g in seq(100,gens,by=100) ) {
  keep = freq_mat[g,]>0.01 & freq_mat[g,]<0.99
  if(sum(keep)>10) {
  df = rbind(df , data.frame( "Type" = stype , "Generation" = rep(g,sum(keep)) , "Frequency" = freq_mat[g,keep] ))
  }
}
}

library("ggplot2")
library(ggridges)
df$Generation = as.factor(df$Generation)
df$Type = as.factor(df$Type)

nb.cols <- length(unique(df$Generation))
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(nb.cols)

ggplot(df, aes(x = Frequency, y = Generation , fill = Generation)) + 
  geom_density_ridges2(show.legend=F,panel_scaling = FALSE, alpha = 0.8 ) + 
  facet_wrap(~Type) + 
  scale_y_discrete(limits=rev) +
  theme_minimal() + scale_fill_manual(values = mycolors) + xlim(c(0,1))

