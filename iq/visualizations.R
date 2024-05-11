library("pheatmap")

# -- kan et al results

# WISC
subtest = c("Vocabulary","Information","Comprehension","Similarities","Arithmetic","Picture\nCompletion","Picture\nArrangement","Block\nDesign","Coding","Digit\nSpan","Object\nAssembly")
subtest_g = c(0.80,0.78,0.68,0.77,0.66,0.53,0.52,0.61,0.35,0.45,0.50)
# WAIS
subtest_g2 = c(0.88,0.85,0.81,0.82,0.74,0.63,0.63,0.67,0.55,0.58,0.55)
culture = c(35,22,15,9,8,3,2,1,0,0,0)

avg_g = subtest_g + subtest_g2
cor.test(rank(avg_g),rank(culture))
x = rank(avg_g)
y = rank(culture)
reg = lm(y~x)
plot(x,y,pch=19,bty="n",xlab="g loading (rank)",las=1,ylab="Culture loading (rank)",xlim=c(0,12),ylim=c(0,12))
abline(reg,lty=3)
summ = summary(reg)
legend("topleft",legend=paste("R2 =",format(summ$r.sq,digits=2),", p =",sprintf("%2.1e",summ$coef[2,4])),bty="n")
text(x,y,subtest,pos=3,cex=0.75)

# --- fawns-richie ukbb
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0231627#sec046

inp = read.csv("fawns_richie_pone.0231627.s004.csv",as.is=T,head=T)
inp = inp[ , grep("UKB",colnames(inp)) ]
inp_time1 = inp[ , grep("Time1",colnames(inp)) ]
inp_time2 = inp[ , grep("Time2",colnames(inp)) ]
keep = apply(is.na(inp_time2),1,sum) == 0
inp_time1 = inp_time1[keep,]
inp_time2 = inp_time2[keep,]

pca1 = prcomp( t(scale(inp_time1)) , center=F , scale=F )
pca2 = prcomp( t(scale(inp_time2)) , center=F , scale=F )
N = nrow(pca1$rotation)
x = rank(pca1$rotation[,1])/N
y = rank(pca2$rotation[,1])/N
shift = abs(rank(pca1$rotation[,1]) - rank(pca2$rotation[,1])) > 10
clr = rep("#666666",length(x))
clr[shift] = "#e6550d"
plot( x , y , pch=19 , bty="n" , col=clr , las=1 , xlab="Test 1 factor (rank)"  , ylab="Test 2 factor (rank)" )
reg = lm(y~x)
abline(reg,lty=3)
summ = summary(reg)
legend("bottomright",legend=paste("R2 =",format(summ$r.sq,digits=2),", p =",sprintf("%2.1e",summ$coef[2,4])),bty="n")

# --- Murayama et al. results
 
library("ggplot2")
tbl = read.csv("Murayama.csv")
colnames(tbl)[1] = "Feature"

df = data.frame("Feature" = c(tbl$Feature,tbl$Feature) , "Grade" = as.factor(c(rep(5,nrow(tbl)),rep(7,nrow(tbl)))) ,"Estimate" = c(tbl$Gr5_Mean,tbl$Gr7_Mean) , "Mean_Estimate" = c(tbl$Gr5_Mean+tbl$Gr7_Mean,tbl$Gr5_Mean+tbl$Gr7_Mean) ,  "SE" = c(tbl$Gr5_SE,tbl$Gr7_SE) , "p" = c(tbl$Gr5_p,tbl$Gr7_p))

ggplot(data=df, aes(x=reorder(Feature, -Mean_Estimate), y=Estimate, fill=Grade)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.7),color="black") + 
  geom_text(aes(label=p),vjust=0.75,position = position_dodge(width=0.7), size=5,hjust=1.5) + 
  coord_flip() + theme_bw() + scale_fill_brewer(palette="Blues") + xlab("")
  
# --- DeLaFuente et al Results

mat = read.csv("delafuente_joined_gws.zscore.csv",as.is=T,head=T)
mat[ abs(mat) < 1.96 ] = 0
mat[ abs(mat) >= 1.96 & abs(mat) < 5.45 ] = 1
mat[ abs(mat) > 5.45 ] = 2
pheatmap( t(mat) , treeheight_row = 0 , treeheight_col = 0 , border="pink" , color = colorRampPalette(c("white", "red"))(100) )

mat = read.csv("delafuente_joined_gws.zscore.csv",as.is=T,head=T)
snps = read.table("delafuente_all.unique.gws",as.is=T)[,1]
qstat = read.table("delafuente_all.unique.gws.gFactor.tab",as.is=T)
qstat = qstat[ match(snps,qstat[,1]) , ]

qstat_sig = qstat[,ncol(qstat)] < 0.05/nrow(qstat)

# -- johnson et al results [PMID 21299289]
outcomes = c("BMI","Constraints","Anxiety","Depression","Diseases","Physical activity","Smoking","Drinking")
# correlations with age 11 IQ
rho_iq = c(-0.14,-0.18,-0.12,-0.11,-0.06,0.08,-0.11,0.12)
# correlations with education
rho_ea = c(-0.11,-0.12,-0.08,-0.09,-0.12,0.13,-0.14,0.12)

# -- de la fuente results
hsq = c(0.155,0.040,0.074,0.110,0.149,0.114,0.212)
gl = c(0.826,0.651,0.308,0.831,0.976,0.853,0.717)

# population level hsq
# Matrix 0.155 (0.040)
# Memory 0.040 (0.002)
# RT 0.074 (0.003)
# Symbol Digit 0.110 (0.008)
# Trails-B 0.149 (0.009)
# Tower 0.114 (0.038)
# VNR 0.212 (0.008)

# loadings
# Matrix 0.826
# Memory 0.651
# RT 0.308
# Symbol Digit 0.831
# Trails-B 0.976
# Tower 0.853
# VNR 0.717