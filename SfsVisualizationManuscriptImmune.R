
#and Immune genes
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/NewBoots/Dplex_All_Immune_sfs", rep = 130)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/NewBoots/Dplex_Mod_sfs", rep = 130)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/NewBoots/Dplex_Eff_sfs", rep = 130)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/NewBoots/Dplex_Sig_sfs", rep = 130)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/NewBoots/Dplex_Rec_sfs", rep = 130)


#can we iterate loading files into a list for bootstrapping
####start of plotting code
###
##
#

#immune
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Ef/")
dpef<- lapply(Sys.glob("Dplex_Eff_out_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Si/")
dpsi<- lapply(Sys.glob("Dplex_Sig_out_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Re/")
dpre<- lapply(Sys.glob("Dplex_Rec_out_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Mo/")
dpmo<- lapply(Sys.glob("Dplex_Mod_out_strap*"), parseOutput)

#And Controls
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Ef_Co/")
dpefc<- lapply(Sys.glob("Eff_Con_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Si_Co/")
dpsic<- lapply(Sys.glob("Sig_Con_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Re_Co/")
dprec<- lapply(Sys.glob("Rec_Con_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_Mo_Co/")
dpmoc<- lapply(Sys.glob("Mod_Con_strap*"), parseOutput)

#modify_depth can be used to subset nested lists, but requires the stupidly named "purrr" package
library("purrr")


#alphs<-unlist(modify_depth(edf,2,"alpha"))
#immune
dpmoa<-unlist(modify_depth(dpmo,2,"alpha"))
dpefa<-unlist(modify_depth(dpef,2,"alpha"))
dpsia<-unlist(modify_depth(dpsi,2,"alpha"))
dprea<-unlist(modify_depth(dpre,2,"alpha"))


dpmoca<-unlist(modify_depth(dpmoc,2,"alpha"))
dpefca<-unlist(modify_depth(dpefc,2,"alpha"))
dpsica<-unlist(modify_depth(dpsic,2,"alpha"))
dpreca<-unlist(modify_depth(dprec,2,"alpha"))

#let's bracket these bois
#sem<-sd(alphs)/sqrt(length(alphs))
#cis<-c(mean(alphs)-2*sem,mean(alphs)+2*sem)



dpmoasem<-sd(dpmoa)/sqrt(length(dpmoa))
dpmoacis<-c(mean(dpmoa)-2*dpmoasem,mean(dpmoa)+2*dpmoasem)

dpefasem<-sd(dpefa)/sqrt(length(dpefa))
dpefacis<-c(mean(dpefa)-2*dpefasem,mean(dpefa)+2*dpefasem)

dpsiasem<-sd(dpsia)/sqrt(length(dpsia))
dpsiacis<-c(mean(dpsia)-2*dpsiasem,mean(dpsia)+2*dpsiasem)

dpreasem<-sd(dprea)/sqrt(length(dprea))
dpreacis<-c(mean(dprea)-2*dpreasem,mean(dprea)+2*dpreasem)


dpmocasem<-sd(dpmoca)/sqrt(length(dpmoca))
dpmocacis<-c(mean(dpmoca)-2*dpmocasem,mean(dpmoca)+2*dpmocasem)

dpefcasem<-sd(dpefca)/sqrt(length(dpefca))
dpefcacis<-c(mean(dpefca)-2*dpefcasem,mean(dpefca)+2*dpefcasem)

dpsicasem<-sd(dpsica)/sqrt(length(dpsica))
dpsicacis<-c(mean(dpsica)-2*dpsicasem,mean(dpsica)+2*dpsicasem)

dprecasem<-sd(dpreca)/sqrt(length(dpreca))
dprecacis<-c(mean(dpreca)-2*dprecasem,mean(dpreca)+2*dprecasem)


#let's do an alpha graphic

# bottom, left, top, right margins
#figure 5
tiff("/Users/Andrew/Documents/SpermCompFigs/Figure5_immune.tiff", units="in", width = 8, height=6, res=300)
par(mai=c(0.06,0.94,0.42,0.06))
dumxd<-c(0.5,0.6,0.7,0.8)
dumxcon<-c(0.51,0.61,0.71,0.81)
#plot(-9,-9,xlim=c(0.45,0.85),ylim=c(0,1),ylab="",xlab="",xaxt="n",
#main="Proportion of adaptive substitutions",panel.first = rect(-9,0.2842102,9,0.4189698, col='lightgray', border=NA),las=1,
#cex.main=2,cex.axis=1.5)
plot(-9,-9,xlim=c(0.45,0.85),ylim=c(0,1),ylab="",xlab="",xaxt="n",
main="Estimated proportion of adaptive substitutions",las=1,cex.main=1.7,cex.axis=1.5)
points(dumxd,c(mean(dprea),mean(dpsia),mean(dpmoa),mean(dpefa)), cex=2)
points(dumxd,c(mean(dprea),mean(dpsia),mean(dpmoa),mean(dpefa)), cex=1.8)
points(dumxcon,c(mean(dpreca),mean(dpsica),mean(dpmoca),mean(dpefca)), pch=15,cex=2)
segments(dumxd,c(dpreacis[1],dpsiacis[1],dpmoacis[1],dpefacis[1]),dumxd,c(dpreacis[2],dpsiacis[2],dpmoacis[2],dpefacis[2]),lwd=2)
segments(dumxd-0.005,c(dpreacis[1],dpsiacis[1],dpmoacis[1],dpefacis[1]),dumxd+0.005,c(dpreacis[1],dpsiacis[1],dpmoacis[1],dpefacis[1]),lwd=2)
segments(dumxcon,c(dprecacis[1],dpsicacis[1],dpmocacis[1],dpefcacis[1]),dumxcon,c(dprecacis[2],dpsicacis[2],dpmocacis[2],dpefcacis[2]),lwd=2)
segments(dumxd-0.005,c(dpreacis[2],dpsiacis[2],dpmoacis[2],dpefacis[2]),dumxd+0.005,c(dpreacis[2],dpsiacis[2],dpmoacis[2],dpefacis[2]),lwd=2)
segments(dumxcon,c(dprecacis[1],dpsicacis[1],dpmocacis[1],dpefcacis[1]),dumxcon,c(dprecacis[2],dpsicacis[2],dpmocacis[2],dpefcacis[2]),lwd=2)
segments(dumxcon-0.005,c(dprecacis[1],dpsicacis[1],dpmocacis[1],dpefcacis[1]),dumxcon+0.005,c(dprecacis[1],dpsicacis[1],dpmocacis[1],dpefcacis[1]),lwd=2)
segments(dumxcon-0.005,c(dprecacis[2],dpsicacis[2],dpmocacis[2],dpefcacis[2]),dumxcon+0.005,c(dprecacis[2],dpsicacis[2],dpmocacis[2],dpefcacis[2]),lwd=2)
mtext(expression(alpha),side=2,las=1,line=3,cex=2.5)
text(0.7,0.76,"Modulation",cex=2)
text(0.8,0.22,"Effector",cex=2)
text(0.6,0.5,"Sigaling",cex=2)
text(0.5,0.6,"Recognition",cex=2)
legend(0.43,0.15, pch=c(1,15),c("Immune","Control"),cex=1.8)
#text(0.6,0.35,"Mongue et al. Background")
dev.off()

#how do we sig test this? MWU?
wilcox.test(dpmoa,dpmoca)
#W = 5214, p-value = 0.01396

wilcox.test(dpefa,dpefca)
#W = 8807, p-value = 1.717e-06

wilcox.test(dpsia,dpsica)
#W = 10297, p-value = 3.063e-15

wilcox.test(dprea,dpreca)
#W = 2430, p-value = 7.447e-16

#can we get the bootstrapped DFE?

#####immune
dpmotops<-length(dpmo)
dpmodfemat<-data.frame(matrix(NA, nrow = dpmotops, ncol = 7))
for(i in 1:dpmotops)
{
dpmodfemat[i,]<-unlist(lapply(dpmo[[i]],getDiscretizedDFE))
}

dpmobsest<-colMeans(dpmodfemat)

names(dpmobsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpmosdsdfe<-apply(dpmodfemat, 2, sd)
dpmosemdfe<-dpmosdsdfe/sqrt(dpmotops)

dumx<-seq(0.5,7.5,1)
dpmociHi<-dpmobsest+2*dpmosemdfe
dpmociLo<-dpmobsest-2*dpmosemdfe


dpeftops<-length(dpef)
dpefdfemat<-data.frame(matrix(NA, nrow = dpeftops, ncol = 7))
for(i in 1:dpeftops)
{
dpefdfemat[i,]<-unlist(lapply(dpef[[i]],getDiscretizedDFE))
}

dpefbsest<-colMeans(dpefdfemat)

names(dpefbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpefsdsdfe<-apply(dpefdfemat, 2, sd)
dpefsemdfe<-dpefsdsdfe/sqrt(dpeftops)

dumx<-seq(0.5,7.5,1)
dpefciHi<-dpefbsest+2*dpefsemdfe
dpefciLo<-dpefbsest-2*dpefsemdfe


dpsitops<-length(dpsi)
dpsidfemat<-data.frame(matrix(NA, nrow = dpsitops, ncol = 7))
for(i in 1:dpsitops)
{
dpsidfemat[i,]<-unlist(lapply(dpsi[[i]],getDiscretizedDFE))
}

dpsibsest<-colMeans(dpsidfemat)

names(dpsibsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpsisdsdfe<-apply(dpsidfemat, 2, sd)
dpsisemdfe<-dpsisdsdfe/sqrt(dpsitops)

dumx<-seq(0.5,7.5,1)
dpsiciHi<-dpsibsest+2*dpsisemdfe
dpsiciLo<-dpsibsest-2*dpsisemdfe

dpretops<-length(dpre)
dpredfemat<-data.frame(matrix(NA, nrow = dpretops, ncol = 7))
for(i in 1:dpretops)
{
dpredfemat[i,]<-unlist(lapply(dpre[[i]],getDiscretizedDFE))
}

dprebsest<-colMeans(dpredfemat)

names(dprebsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpresdsdfe<-apply(dpredfemat, 2, sd)
dpresemdfe<-dpresdsdfe/sqrt(dpretops)

dumx<-seq(0.5,7.5,1)
dpreciHi<-dprebsest+2*dpresemdfe
dpreciLo<-dprebsest-2*dpresemdfe

##Immune Controls 

dpmoctops<-length(dpmoc)
dpmocdfemat<-data.frame(matrix(NA, nrow = dpmoctops, ncol = 7))
for(i in 1:dpmoctops)
{
dpmocdfemat[i,]<-unlist(lapply(dpmoc[[i]],getDiscretizedDFE))
}

dpmocbsest<-colMeans(dpmocdfemat)

names(dpmocbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpmocsdsdfe<-apply(dpmocdfemat, 2, sd)
dpmocsemdfe<-dpmocsdsdfe/sqrt(dpmoctops)

dumx<-seq(0.5,7.5,1)
dpmocciHi<-dpmocbsest+2*dpmocsemdfe
dpmocciLo<-dpmocbsest-2*dpmocsemdfe


dpefctops<-length(dpefc)
dpefcdfemat<-data.frame(matrix(NA, nrow = dpefctops, ncol = 7))
for(i in 1:dpefctops)
{
dpefcdfemat[i,]<-unlist(lapply(dpefc[[i]],getDiscretizedDFE))
}

dpefcbsest<-colMeans(dpefcdfemat)

names(dpefcbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpefcsdsdfe<-apply(dpefcdfemat, 2, sd)
dpefcsemdfe<-dpefcsdsdfe/sqrt(dpefctops)

dumx<-seq(0.5,7.5,1)
dpefcciHi<-dpefcbsest+2*dpefcsemdfe
dpefcciLo<-dpefcbsest-2*dpefcsemdfe


dpsictops<-length(dpsic)
dpsicdfemat<-data.frame(matrix(NA, nrow = dpsictops, ncol = 7))
for(i in 1:dpsictops)
{
dpsicdfemat[i,]<-unlist(lapply(dpsic[[i]],getDiscretizedDFE))
}

dpsicbsest<-colMeans(dpsicdfemat)

names(dpsicbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpsicsdsdfe<-apply(dpsicdfemat, 2, sd)
dpsicsemdfe<-dpsicsdsdfe/sqrt(dpsictops)

dumx<-seq(0.5,7.5,1)
dpsicciHi<-dpsicbsest+2*dpsicsemdfe
dpsicciLo<-dpsicbsest-2*dpsicsemdfe

dprectops<-length(dprec)
dprecdfemat<-data.frame(matrix(NA, nrow = dprectops, ncol = 7))
for(i in 1:dprectops)
{
dprecdfemat[i,]<-unlist(lapply(dprec[[i]],getDiscretizedDFE))
}

dprecbsest<-colMeans(dprecdfemat)

names(dprecbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dprecsdsdfe<-apply(dprecdfemat, 2, sd)
dprecsemdfe<-dprecsdsdfe/sqrt(dprectops)

dumx<-seq(0.5,7.5,1)
dprecciHi<-dprecbsest+2*dprecsemdfe
dprecciLo<-dprecbsest-2*dprecsemdfe


#the DFE comparison plots

library("plotrix")
#we have to define our gap bounds
from<-0.175
to<-0.675

effcols<-c("brown4","brown3","brown1","gray66","darkseagreen2","darkseagreen3","darkseagreen4")



#######immune DFE plot

from<-0.27
to<-0.64
par(mfrow=c(4,1))
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.42,0.06))
gap.barplot(dpmobsest,xlim=c(0.5,7.4),main="Immune DFEs",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.45),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpmociHi,dumx+0.5,dpmociLo,lwd=2)
segments(ebx1+0.5,dpmociHi,ebx2+0.5,dpmociHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpmociHi[1]-(to-from),ebx2[1]+0.5,dpmociHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpmociHi[1]-(to-from),dumx[1]+0.5,dpmociLo[1]-(to-from),lwd=2)
text(3,0.22,"Modulation", cex=1.5, font=2)

#from<-0.278
#to<-0.6
gap.barplot(dpefbsest,xlim=c(0.5,7.4),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.45),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpefciHi,dumx+0.5,dpefciLo,lwd=2)
segments(ebx1+0.5,dpefciHi,ebx2+0.5,dpefciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpefciHi[1]-(to-from),ebx2[1]+0.5,dpefciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpefciHi[1]-(to-from),dumx[1]+0.5,dpefciLo[1]-(to-from),lwd=2)
text(3,0.35,"Effectors", cex=1.5, font=2)


#from<-0.2
#to<-0.8
gap.barplot(dpsibsest,xlim=c(0.5,7.4),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.45),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpsiciHi,dumx+0.5,dpsiciLo,lwd=2)
segments(ebx1+0.5,dpsiciHi,ebx2+0.5,dpsiciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpsiciHi[1]-(to-from),ebx2[1]+0.5,dpsiciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpsiciHi[1]-(to-from),dumx[1]+0.5,dpsiciLo[1]-(to-from),lwd=2)
text(3,0.33,"Signaling", cex=1.5, font=2)

#from<-0.15
#to<-0.7
gap.barplot(dprebsest,xlim=c(0.5,7.4),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.45),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpreciHi,dumx+0.5,dpreciLo,lwd=2)
segments(ebx1+0.5,dpreciHi,ebx2+0.5,dpreciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpreciHi[1]-(to-from),ebx2[1]+0.5,dpreciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpreciHi[1]-(to-from),dumx[1]+0.5,dpreciLo[1]-(to-from),lwd=2)
text(3,0.33,"Recognition", cex=1.5, font=2)


#######immune DFE plots zoomed in on pos/neutral variation

par(mfrow=c(4,1))
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.42,0.06))
gap.barplot(dpmobsest,xlim=c(2.67,7.4),main="Immune DFEs",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.1),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpmociHi,dumx+0.5,dpmociLo,lwd=2)
segments(ebx1+0.5,dpmociHi,ebx2+0.5,dpmociHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpmociHi[1]-(to-from),ebx2[1]+0.5,dpmociHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpmociHi[1]-(to-from),dumx[1]+0.5,dpmociLo[1]-(to-from),lwd=2)
text(6.5,0.08,"Modulation", cex=1.5, font=2)


gap.barplot(dpefbsest,xlim=c(2.67,7.4),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.1),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpefciHi,dumx+0.5,dpefciLo,lwd=2)
segments(ebx1+0.5,dpefciHi,ebx2+0.5,dpefciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpefciHi[1]-(to-from),ebx2[1]+0.5,dpefciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpefciHi[1]-(to-from),dumx[1]+0.5,dpefciLo[1]-(to-from),lwd=2)
text(6.5,0.08,"Effectors", cex=1.5, font=2)



gap.barplot(dpsibsest,xlim=c(2.67,7.4),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.1),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpsiciHi,dumx+0.5,dpsiciLo,lwd=2)
segments(ebx1+0.5,dpsiciHi,ebx2+0.5,dpsiciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpsiciHi[1]-(to-from),ebx2[1]+0.5,dpsiciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpsiciHi[1]-(to-from),dumx[1]+0.5,dpsiciLo[1]-(to-from),lwd=2)
text(6.5,0.08,"Signaling", cex=1.5, font=2)

gap.barplot(dprebsest,xlim=c(2.67,7.4),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.1),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpreciHi,dumx+0.5,dpreciLo,lwd=2)
segments(ebx1+0.5,dpreciHi,ebx2+0.5,dpreciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpreciHi[1]-(to-from),ebx2[1]+0.5,dpreciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpreciHi[1]-(to-from),dumx[1]+0.5,dpreciLo[1]-(to-from),lwd=2)
text(6.5,0.08,"Recognition", cex=1.5, font=2)


#####Side-by-side with controls
#######immune DFE plot
#Figure 4
tiff("/Users/Andrew/Documents/SpermCompFigs/Figure4_immune.tiff", units="in", width = 8.5, height=9, res=300)
from<-0.27
to<-0.64
par(mfcol=c(4,2))
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.42,0.06))
gap.barplot(dprebsest,xlim=c(0.5,6.3),main="Immune DFEs",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpreciHi,dumx+0.5,dpreciLo,lwd=2)
segments(dumx+0.55,dpreciHi,dumx+0.45,dpreciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpreciHi[1]-(to-from),dumx[1]+0.45,dpreciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpreciHi[1]-(to-from),dumx[1]+0.5,dpreciLo[1]-(to-from),lwd=2)
text(4.5,0.36,"Recognition", cex=2.5, font=2)

#from<-0.278
#to<-0.6
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.12,0.06))
gap.barplot(dpsibsest,xlim=c(0.5,6.3),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpsiciHi,dumx+0.5,dpsiciLo,lwd=2)
segments(dumx+0.55,dpsiciHi,dumx+0.45,dpsiciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpsiciHi[1]-(to-from),dumx[1]+0.45,dpsiciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpsiciHi[1]-(to-from),dumx[1]+0.5,dpsiciLo[1]-(to-from),lwd=2)
text(4.5,0.36,"Signaling", cex=2.5, font=2)

#from<-0.2
#to<-0.8
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.12,0.06))
gap.barplot(dpmobsest,xlim=c(0.5,6.3),main="",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpmociHi,dumx+0.5,dpmociLo,lwd=2)
segments(dumx+0.55,dpmociHi,dumx+0.45,dpmociHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpmociHi[1]-(to-from),dumx[1]+0.45,dpmociHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpmociHi[1]-(to-from),dumx[1]+0.5,dpmociLo[1]-(to-from),lwd=2)
text(4.5,0.36,"Modulation", cex=2.5, font=2)

#from<-0.15
#to<-0.7
# bottom, left, top, right margins
par(mai=c(0.36,0.64,0.12,0.06))
gap.barplot(dpefbsest,xlim=c(0.5,6.3),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="Proportion",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpefciHi,dumx+0.5,dpefciLo,lwd=2)
segments(dumx+0.55,dpefciHi,dumx+0.45,dpefciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpefciHi[1]-(to-from),dumx[1]+0.45,dpefciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpefciHi[1]-(to-from),dumx[1]+0.5,dpefciLo[1]-(to-from),lwd=2)
text(4.5,0.36,"Effectors", cex=2.5, font=2)
mtext(side=1,text = c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)"), 
at = c(0.44,1.7,3.02,4,5.05,6.08),line=1,cex=0.95)
#######immune controls DFE plot
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.42,0.06))
gap.barplot(dprecbsest,xlim=c(0.5,6.3),main="Control DFEs",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dprecciHi,dumx+0.5,dprecciLo,lwd=2)
segments(dumx+0.55,dprecciHi,dumx+0.45,dprecciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dprecciHi[1]-(to-from),dumx[1]+0.45,dprecciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dprecciHi[1]-(to-from),dumx[1]+0.5,dprecciLo[1]-(to-from),lwd=2)
#text(3,0.33,"Control", cex=1.5, font=2)

#from<-0.278
#to<-0.6
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.12,0.06))
gap.barplot(dpsicbsest,xlim=c(0.5,6.3),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpsicciHi,dumx+0.5,dpsicciLo,lwd=2)
segments(dumx+0.55,dpsicciHi,dumx+0.45,dpsicciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpsicciHi[1]-(to-from),dumx[1]+0.45,dpsicciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpsicciHi[1]-(to-from),dumx[1]+0.5,dpsicciLo[1]-(to-from),lwd=2)
#text(3,0.33,"Control", cex=1.5, font=2)



#from<-0.2
#to<-0.8
# bottom, left, top, right margins
par(mai=c(0.06,0.64,0.12,0.06))
gap.barplot(dpmocbsest,xlim=c(0.5,6.3),main="",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpmocciHi,dumx+0.5,dpmocciLo,lwd=2)
segments(dumx+0.55,dpmocciHi,dumx+0.45,dpmocciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpmocciHi[1]-(to-from),dumx[1]+0.45,dpmocciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpmocciHi[1]-(to-from),dumx[1]+0.5,dpmocciLo[1]-(to-from),lwd=2)
#text(3,0.22,"Control", cex=1.5, font=2)

#from<-0.15
#to<-0.7
# bottom, left, top, right margins
par(mai=c(0.36,0.64,0.12,0.06))
gap.barplot(dpefcbsest,xlim=c(0.5,6.3),main=" ",gap=c(from,to),col=effcols,ytics=seq(0,1,0.1),
xaxlab=rep("",7),xlab="",ylim=c(0,0.55),ylab="",xaxt="n",cex.lab=1.5,las=1,cex.axis=1.4,cex.main=1.8)
segments(dumx+0.5,dpefcciHi,dumx+0.5,dpefcciLo,lwd=2)
segments(dumx+0.55,dpefcciHi,dumx+0.45,dpefcciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(dumx[1]+0.55,dpefcciHi[1]-(to-from),dumx[1]+0.45,dpefcciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpefcciHi[1]-(to-from),dumx[1]+0.5,dpefcciLo[1]-(to-from),lwd=2)
#text(3,0.35,"Control", cex=1.5, font=2)
mtext(side=1,text = c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)"), 
at = c(0.44,1.7,3.02,4,5.05,6.08),line=1,cex=0.95)
dev.off()

