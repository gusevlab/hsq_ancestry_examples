library('metafor')
library("RColorBrewer")

clr = brewer.pal(n=4,"Set1")
clr[3] = "#cccccc"

files = c("schmidt_hunter_growth_afqt.csv","schmidt_hunter_growth_work_sample.csv","schmidt_hunter_growth_knowledge.csv","schmidt_hunter_growth_supervisor.csv")
names = c("IQ test (AFQT)","Work Sample","Job Knowledge","Supervisor Rating")
xvals = c(mean(c(1,4)),mean(c(5,11)),mean(c(12,16)),mean(c(17,22)),mean(c(23,60)),61)
par(mfrow=c(1,4))
for ( f in 1:length(files) ) {
  
tbl = read.csv(paste("~/Dropbox (Partners HealthCare)/TWT/",files[f],sep=''),as.is=T,head=T)

top = grep("Upper" ,colnames(tbl))
bottom = grep("Lower" ,colnames(tbl))

rows = which(tbl[,1] == "M")

est = matrix(nrow=length(rows),ncol=2)
est_se = matrix(nrow=length(rows),ncol=2)

all_mean = matrix(nrow=length(rows),ncol=1)
all_sd = matrix(nrow=length(rows),ncol=1)

ctr = 1
for ( r in rows ) {
  
  all_mean[ctr,1] = mean(as.numeric(unlist(tbl[r,c(bottom,top)])))
  # under normal distribution, top/bottom half have SDs 0.605 smaller than the full distribution
  all_sd[ctr,1] = mean(as.numeric(unlist(tbl[r+1,c(bottom,top)])))/0.605
  
  df = data.frame( "yi" = as.numeric(unlist(tbl[r,top])) , "vi" = (as.numeric(unlist(tbl[(r+1),top])) / sqrt(as.numeric(unlist(tbl[(r+2),top]))))^2 )
  reg = rma(yi,vi,data=df,method="FE")
  
  df = data.frame( "yi" = as.numeric(unlist(tbl[r,bottom])) , "vi" = (as.numeric(unlist(tbl[(r+1),bottom])) / sqrt(as.numeric(unlist(tbl[(r+2),bottom]))))^2 )
  reg_bottom = rma(yi,vi,data=df,method="FE")
  cat( reg$beta , reg$se , reg_bottom$beta , reg_bottom$se , '\n' )
  est[ctr,1] = reg$beta
  est[ctr,2] = reg_bottom$beta
  
  est_se[ctr,1] = reg$se
  est_se[ctr,2] = reg_bottom$se
  ctr = ctr + 1
}

yrng = range(c(est+est_se,est-est_se,all_mean[,1]+all_sd[,1],all_mean[,1]-all_sd[,1]))
if ( f != 4 ) yrng = c(35,65)
else yrng = c(0,100)

plot(0 , 0 , xaxt="n", xlab="" , ylab="Score", bty="n" , type="n" , xlim=range(xvals) , ylim=yrng , las=2, main=names[f])
axis(1,at=xvals,lab=c("1-4m","5-11m","12-16m","17-22m","23-60m","61+"),las=2)
#lines(1:6,all_mean[,1],col=clr[3],lwd=2)
#points(1:6,all_mean[,1],col=clr[3],pch=19)
#polygon( x = c(1:6,rev(1:6)) , y = c( all_mean[,1] + all_sd[,1], rev(all_mean[,1] - all_sd[,1])) , border=paste(clr[3],"50",sep='') , col=paste(clr[3],"25",sep='') )

lines(xvals,est[,1],col=clr[1],lwd=2)
points(xvals,est[,1],col=clr[1],pch=19)
polygon( x = c(xvals,rev(xvals)) , y = c( est[,1] + 1.96*est_se[,1], rev(est[,1] - 1.96*est_se[,1])) , border=NA , col=paste(clr[1],"25",sep='') )
lines(xvals,est[,2],col=clr[2],lwd=2)
points(xvals,est[,2],col=clr[2],pch=19)
polygon( x = c(xvals,rev(xvals)) , y = c( est[,2] + 1.96*est_se[,2], rev(est[,2] - 1.96*est_se[,2])) , border=NA , col=paste(clr[2],"25",sep='') )

}