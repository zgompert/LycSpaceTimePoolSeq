## PCA and Fst by chomromsome 
## this code was written for the first two sets of sequence data

library(data.table)

a1f<-list.files(pattern="ad1_filt")
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_filt_lycpool_chrom","",a1f)
chrom<-gsub(".txt","",temp)

pdf("WG_LG_PCAs.pdf",width=7,height=10.5)
par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.3;ca<-1.1;cm<-1.4
for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n<-a1+a2
	p<-a2/(a1+a2) ## non-ref

	p[is.na(p)]<-0.001
	## pca
	pc<-prcomp(t(p),center=TRUE,scale=FALSE)
	o<-summary(pc)
	pct<-round(o$importance[2,1:3] * 100,1)

	plot(pc$x[,1],pc$x[,2],pch=rep(15:20,each=7),col=ids[,2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC2 (",pct[2],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("LG ",chrom[i],sep=""),cex.main=cm)
	#text(pc$x[,1],pc$x[,2],ids[,1],cex=.7)
	xa<-min(pc$x[,1]) * .9
	ya<-pc$x[22,2] * .9
	if(ya > 0){
		legend(xa,ya,ids[,1],pch=rep(15:20,each=7),col=ids[,2],ncol=3,cex=.6)
	}
	plot(pc$x[,1],pc$x[,3],pch=rep(15:20,each=7),col=ids[,2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC3 (",pct[3],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("LG ",chrom[i],sep=""),cex.main=cm)
	#text(pc$x[,1],pc$x[,3],ids[,1],cex=.7)

}
dev.off()

## Fst
P<-vector("list",24)
n<-vector("list",24)
H<-vector("list",24)
for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n[[i]]<-a1+a2
	P[[i]]<-a2/(a1+a2) ## non-ref
	H[[i]]<-2*p*(1-p)
}
	
save(list=ls(),file="fst.rdat")

Npop<-35
Nx<-(Npop*(Npop-1))/2
fstGw<-matrix(NA,nrow=Nx,ncol=24)
fst90<-matrix(NA,nrow=Nx,ncol=24)
x<-1
for(i in 1:(Npop-1)){
	for(j in (i+1):Npop){
		for(k in 1:24){
			keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
			pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
			Ht<-2*pbar*(1-pbar)
			Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
			
			fstGw[x,k]<-mean(Ht-Hs)/mean(Ht)
			fst<-(Ht-Hs)/Ht
			fst90[x,k]<-quantile(fst,.9)
		}
		x<-x+1
	}
}

mord<-c(1,12,18:24,2:11,13:17)
pdf("FstLycPool.pdf",width=8,height=10)
par(mfrow=c(5,3))
par(mar=c(4,5,2.5,.5))
x<-1
for(i in 1:(Npop-1)){
	for(j in (i+1):Npop){
		plot(fstGw[x,mord],pch=19,ylim=c(0,1),xlab="Chromosome",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=.9)
		segments(x0=1:24,y0=fstGw[x,mord],x1=1:24,y1=fst90[x,mord])
		title(main=paste(ids[i,1]," x ",ids[j,1],sep=""),cex.main=1.2)
		mn<-round(mean(fstGw[x,mord]),2)
		text(5,.9,mn)
		x<-x+1
	}
}
dev.off()

### window examples

pdf("FstWinsABM20_SIN10.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-1;j<-28 
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsABM20_GNP17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-1;j<-18
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsMR20_YG20.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-23;j<-34
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsMR20_CP19.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-23;j<-13
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsBCR17_BTB17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-2;j<-9
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsBTB10_BTB17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-6;j<-9
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsBTB17_BTB17rep.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-9;j<-10
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsYG98_YG20.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-34;j<-35
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()
