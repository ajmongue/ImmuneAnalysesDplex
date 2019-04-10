#first some file loading and data slicing
library("data.table")
danmaster<-read.csv(file="/Users/Andrew/Downloads/evoconfiles/NewKanzen.csv",header=T,stringsAsFactors=F)
dan_append<-read.csv(file="/Users/Andrew/Documents/danausorthologsformanduca.csv",header=T)
dan_immune<-read.csv(file="/Users/Andrew/Downloads/All_Monarch_Immune_genes_trim_autosomal.csv",header=T,stringsAsFactors=F)
colnames(dan_append)<-c("Gene","ZorA","Chromosome","Manduca.Sperm.Ortho")
dan_ctrl<-read.csv(file="/Users/Andrew/Documents/Immune_Control_IDs.csv",header=T,stringsAsFactors=F)
danmaster<-as.data.frame(danmaster,stringsAsFactors=F)
dantranscend<-merge(danmaster,dan_append,by=c("Gene","Chromosome","ZorA"),all.x=T, all.y=T)
danomega<-merge(dantranscend,dan_immune, by="Gene", all.x=T, all.y=T)
dan_class<-danomega[,19]
dan_class[is.na(dan_class)]<-"Non-immune"
danomega[,19]<-dan_class
danomega<-merge(dan_ctrl,danomega, by="Gene", all.x=T, all.y=T)


#now let's break down the Dn/Ds by immune class

#Recognition
#Dn
wilcox.test(((reco$dN/reco$Non.sites)),((recog$dN/recog$Non.sites)))
W = 1774, p-value = 0.09668
#Ds
wilcox.test(((reco$dS/reco$Syn.sites)),((recog$dS/recog$Syn.sites)))
#W = 2208, p-value = 0.7972
#Dn/Ds
wilcox.test(((reco$dN/reco$Non.sites)/(reco$dS/reco$Syn.sites)),((recog$dN/recog$Non.sites)/(recog$dS/recog$Syn.sites)))
#W = 1156, p-value = 0.01816
#significant!


#Modulation
#Dn
wilcox.test(((moco$dN/moco$Non.sites)),((modul$dN/modul$Non.sites)))
#W = 11062, p-value = 0.3125
#Ds
wilcox.test(((moco$dS/moco$Syn.sites)),((modul$dS/modul$Syn.sites)))
#W = 14112, p-value = 0.2217
#Dn/Ds
wilcox.test(((moco$dN/moco$Non.sites)/(moco$dS/moco$Syn.sites)),((modul$dN/modul$Non.sites)/(modul$dS/modul$Syn.sites)))
#W = 6036.5, p-value = 0.008601
#Hark, a significant result!


#Effector
#Dn
wilcox.test(((efco$dN/efco$Non.sites)),((effec$dN/effec$Non.sites)))
#W = 7918.5, p-value = 0.4833
#Ds
wilcox.test(((efco$dS/efco$Syn.sites)),((effec$dS/effec$Syn.sites)))
#W = 7243, p-value = 0.9559
#Dn/Ds
wilcox.test(((efco$dN/efco$Non.sites)/(efco$dS/efco$Syn.sites)),((effec$dN/effec$Non.sites)/(effec$dS/effec$Syn.sites)))
#W = 2866, p-value = 0.7635



#Signalling 
#Dn
wilcox.test(((sico$dN/sico$Non.sites)),((signa$dN/signa$Non.sites)))
#W = 27713, p-value = 0.3633
#Ds
wilcox.test(((sico$dS/sico$Syn.sites)),((signa$dS/signa$Syn.sites)))
#W = 23423, p-value = 0.01393
#significant

#Dn/Ds
wilcox.test(((sico$dN/sico$Non.sites)/(sico$dS/sico$Syn.sites)),((signa$dN/signa$Non.sites)/(signa$dS/signa$Syn.sites)))
#W = 23427, p-value = 0.352


#now let's plot it
#Figure 3
# bottom, left, top, right margins
tiff("/Users/Andrew/Documents/SpermCompFigs/Figure3_immune.tiff", units="in", width = 8.5, height=6, res=300)
par(mai=c(0.86,0.94,0.42,0.06))
boxplot((dan_bg$dN/dan_bg$Non.sites)/(dan_bg$dS/dan_bg$Syn.sites),outline=F, col=, at = -1, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,5.3), ylim=c(0,0.4),cex.axis=2,cex.lab=1.8,cex.main=1.8,
main="Divergence of immune genes compared to controls",ylab="Dn/Ds",xlab="")
axis(1, at=c(1.25, 2.45, 3.65, 4.85),labels=F)
text(x=c(1.25, 2.45, 3.65, 4.85),y=-0.035,par("usr")[3],labels=c("Recognition","Signaling","Modulation","Effector"),srt=30,pos=2,xpd=T,cex=1.2)
boxplot((recog$dN/recog$Non.sites)/(recog$dS/recog$Syn.sites), add=T, at = 1,outline=F,notch=F, las=1,cex.axis=1.3)
boxplot((reco$dN/reco$Non.sites)/(reco$dS/reco$Syn.sites), add=T, at = 1.5,outline=F,notch=F, las=1,cex.axis=1.3,col="grey")
boxplot((modul$dN/modul$Non.sites)/(modul$dS/modul$Syn.sites), add=T, at = 3.4,outline=F, notch=F, las=1,cex.axis=1.3)
boxplot((moco$dN/moco$Non.sites)/(moco$dS/moco$Syn.sites), add=T, at = 3.9,outline=F, notch=F, las=1,cex.axis=1.3,col="grey")
boxplot((effec$dN/effec$Non.sites)/(effec$dS/effec$Syn.sites), add=T, at = 4.6,outline=F, notch=F, las=1,cex.axis=1.3)
boxplot((efco$dN/efco$Non.sites)/(efco$dS/efco$Syn.sites), add=T, at = 5.1,outline=F, notch=F, las=1,cex.axis=1.3,col="grey")
boxplot((signa$dN/signa$Non.sites)/(signa$dS/signa$Syn.sites), add=T, at = 2.2,outline=F, notch=F, las=1,cex.axis=1.3)
boxplot((sico$dN/sico$Non.sites)/(sico$dS/sico$Syn.sites), add=T, at = 2.7,outline=F, notch=F, las=1,cex.axis=1.3,col="grey")
#text(1,0.12,"Rec",cex=1.4)
text(1,0.305,"*",cex=3.5)
#text(3.4,0.14,"Mod",cex=1.4)
text(3.4,0.29,"*",cex=3.5)
text(3.4,0.32,"*",cex=3.5)
#text(4.6,0.07,"Eff",cex=1.4)
#text(2.2,0.07,"Sig",cex=1.4)
#mtext(side=1, line=4, adj=0.5, "Gene Class", cex=2)
dev.off()

library("vioplot")
vioplot((recog$dN/recog$Non.sites)/(recog$dS/recog$Syn.sites), at = 1,outline=F,notch=T, las=1,cex.axis=1.3,na.rm=T)


