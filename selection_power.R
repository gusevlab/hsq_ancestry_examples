# Visualizing the power to distinguish selection from drift

library(gridExtra)
library(pwr)
library(pheatmap)
library(RColorBrewer)

pdf("selection_power.pdf",width=20,height=6)

total_gens = as.integer(65e3/30)
gens = total_gens
Ne = 10e3

p0 = seq(0.05,0.95,by=0.05)
samp_0 = 1e3
samp_t = 1e3
all_p1 = seq(0.05,0.95,by=0.05)
power = matrix(NA,nrow=length(p0),ncol=length(all_p1))
rownames(power) = sprintf("%1.2f",p0)
colnames(power) = sprintf("%1.2f",all_p1)
s = power

for ( i in 1:length(all_p1) ) {
  p1 = all_p1[i]
  drift_var0 = p0 * (1-p0) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  drift_var1 = p1 * (1-p1) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  chisq = (p1 - p0)^2 / ( drift_var0 + drift_var1)
  power[,i] = pwr.chisq.test(w=chisq,N=1,df=1,sig.level=0.05)$power
  
  # from 5.3d
  s[,i] = (1/gens)*log( p1*(1-p0) / (p0*(1-p1)))
}

s_lab = matrix("",nrow=nrow(s),ncol=ncol(s))
s_lab[ power > 0.2 ] = sprintf("%1.2f",power[ power > 0.2 ])

plot1 = pheatmap(power,cluster_cols=F,cluster_rows=F,breaks=seq(0,1,by=0.1) , color = colorRampPalette((brewer.pal(n = 7, name = "Blues")))(10),display_numbers=s_lab,number_color="black")

s_lab = matrix("",nrow=nrow(s),ncol=ncol(s))
s_lab[ abs(s) > 0.0005] = "*"
s_lab[ abs(s) < 0.0005 & abs(s) > 0.0005] = "**"
diag(s_lab) = ""

#p2 = pheatmap( -log10(abs(s)) ,cluster_cols=F,cluster_rows=F,breaks=seq(2,5,length.out=10) , color = colorRampPalette(rev(brewer.pal(n = 9, name = "Blues")))(10),display_numbers=s_lab)

plot2 = pheatmap( -log10(abs(s)) ,cluster_cols=F,cluster_rows=F,breaks=seq(2,5,length.out=10) , color = c("#2d004b","#2d004b","#2d004b","#542788","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffff","#ffffff"))

# ---- add second plot

p0 = seq(0.05,0.95,by=0.05)
samp_0 = 100
samp_t = 100
all_p1 = seq(0.05,0.95,by=0.05)
power = matrix(NA,nrow=length(p0),ncol=length(all_p1))
rownames(power) = sprintf("%1.2f",p0)
colnames(power) = sprintf("%1.2f",all_p1)
s = power

for ( i in 1:length(all_p1) ) {
  p1 = all_p1[i]
  drift_var0 = p0 * (1-p0) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  drift_var1 = p1 * (1-p1) * ( 1/(2*samp_0) + 1/(2*samp_t) + (gens / (2*Ne))*(1-1/(2*samp_t)) )
  chisq = (p1 - p0)^2 / ( drift_var0 + drift_var1)
  power[,i] = pwr.chisq.test(w=chisq,N=1,df=1,sig.level=0.05)$power
  
  # from 5.3d
  s[,i] = (1/gens)*log( p1*(1-p0) / (p0*(1-p1)))
}

s_lab = matrix("",nrow=nrow(s),ncol=ncol(s))
s_lab[ power > 0.2 ] = sprintf("%1.2f",power[ power > 0.2 ])

plot3 = pheatmap(power,cluster_cols=F,cluster_rows=F,breaks=seq(0,1,by=0.1) , color = colorRampPalette((brewer.pal(n = 7, name = "Blues")))(10),display_numbers=s_lab,number_color="black")

# ---

plot_list=list()
plot_list[['p1']]=plot3[[4]]
plot_list[['p2']]=plot1[[4]]
plot_list[['p3']]=plot2[[4]]
grid.arrange(grobs=plot_list, ncol=3)

dev.off()