tbl = read.csv("~/Dropbox (Partners HealthCare)/TWT/LW/FINAL/iq/knyspel_plomin_ace.csv")

tbl$Loadings_BiFactor = c(0.610,
                          0.435,
                          0.671,
                          0.625,
                          0.715,
                          0.703,
                          0.393,
                          0.431,
                          0.679,
                          0.665,
                          0.590,
                          0.585,
                          0.641,
                          0.539)

tbl$Loadings_BiFactor_Twins = c(0.598,
                                0.454,
                                0.689,
                                0.656,
                                0.735,
                                0.695,
                                0.438,
                                0.454,
                                0.702,
                                0.661,
                                0.619,
                                0.590,
                                0.678,
                                0.576)

tbl$BiFactor_A = c(0.361,
                   0.262,
                   0.622,
                   0.581,
                   0.579,
                   0.552,
                   0.438,
                   0.413,
                   0.495,
                   0.361,
                   0.440,
                   0.463,
                   0.673,
                   0.530)


options(ggrepel.max.overlaps = Inf)

p1 <- ggplot(tbl, aes(x=BiFactor_A, y=Network_A, label = Measure)) +
  geom_point() + geom_text_repel() + xlab("Subtest heritability (bifactor)") + ylab("Subtest heritability (network)") +
  theme_bw() + theme(legend.position = "none")

p2 <- ggplot(tbl, aes(Loadings_BiFactor_Twins, BiFactor_A, label = Measure)) +
  geom_smooth(method = "lm") + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = "none") + xlab("g loading (bifactor)") + ylab("Subtest heritability (bifactor)")

p3 <- ggplot(tbl, aes(Loadings_BiFactor_Twins, Network_A, label = Measure)) +
  geom_smooth(method = "lm") + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = "none") + xlab("g loading (bifactor)") + ylab("Subtest heritability (network)")

gridExtra::grid.arrange(p1, p2, p3 , ncol = 3)

par(mfrow=c(1,3))

x = tbl$BiFactor_A
y = tbl$Network_A
plot(x,y,xlim=c(0,1),ylim=c(0,1),las=1,xlab="A (Bifactor model)", ylab="A (network model)",main="Heritability",bty="n")
text(x,y,tbl$Measure,pos=4,cex=0.8)
reg = lm(y ~ x)
abline(reg,lty=3)
reg = cor.test(y,x)
legend("topleft",legend=c(paste("Cor = ",sprintf("%2.2f",reg$est)),paste("p = ",sprintf("%2.2f",reg$p.value),sep='')) ,bty="n")

x = tbl$Loadings_BiFactor_Twins
y = tbl$BiFactor_A
plot(x,y,xlim=c(0,1),ylim=c(0,1),las=1,ylab="A (Bifactor model)", xlab="g loading (Bifactor model)",main="Factor heritability versus g loading",bty="n")
text(x,y,tbl$Measure,pos=4,cex=0.8)
reg = lm(y ~ x)
abline(reg,lty=3)
reg = cor.test(y,x)
legend("topleft",legend=c(paste("Cor = ",sprintf("%2.2f",reg$est)),paste("p = ",sprintf("%2.2f",reg$p.value),sep='')) ,bty="n")

x = tbl$Loadings_BiFactor_Twins
y = tbl$Network_A
plot(x,y,xlim=c(0,1),ylim=c(0,1),las=1,ylab="A (network model)", xlab="g loading (Bifactor model)",main="Network heritability versus g loading",bty="n")
text(x,y,tbl$Measure,pos=4,cex=0.8)
reg = lm(y ~ x)
abline(reg,lty=3)
reg = cor.test(y,x)
legend("topleft",legend=c(paste("Cor = ",sprintf("%2.2f",reg$est)),paste("p = ",sprintf("%2.2f",reg$p.value),sep='')) ,bty="n")

cor.test(tbl$Loadings_BiFactor_Twins,tbl$BiFactor_A)
cor.test(tbl$Loadings_BiFactor_Twins,tbl$Network_A)