## prelminary analyses examing patterns of change over time
library(data.table)

a1f<-list.files(pattern="ad1_filt")
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
N<-length(a1f)
ids<-read.table("IDS.txt",header=FALSE)
temp<-gsub("ad1_filt_lycpool_chrom","",a1f)
chrom<-gsub(".txt","",temp)

######################################################
### patterns of change ##############################
deg2rad <- function(deg) return(deg*pi/180)
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

coords<-read.table("coords.txt",header=FALSE)
lon<-deg2rad(coords[,3])
lat<-deg2rad(coords[,2])
d<-matrix(NA,nrow=N,ncol=N)
for(i in 1:N){for(j in 1:N){
        d[i,j]<-gcd.slc(lon[i],lat[i],lon[j],lat[j])
}}

taxon<-c("mel","jh","jh","jh","jh","sn","idas","jh","sn","jh","jh","jh","jh")
d_tax<-matrix(NA,nrow=N,ncol=N)
for(i in 1:N){for(j in 1:N){
        d_tax[i,j]<-as.numeric(taxon[i]==taxon[j])
}}



c_dp<-matrix(NA,nrow=24,ncol=2)
c_arcdp<-matrix(NA,nrow=24,ncol=2)

pdf("ChangeCorsDist.pdf",width=8,height=8)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

for(x in 1:24){
        a1<-as.matrix(fread(a1f[x],header=F))
        a2<-as.matrix(fread(a2f[x],header=F))
        n<-a1+a2
        p<-a2/(a1+a2) ## non-ref
        p[is.na(p)]<-0.001
        colnames(p)<-ids[,1]

        keep<-which(apply(n,1,min) > 40)
        combos<-matrix(c(2,1, ## abm 20, 12
          11,5, ## bcr 22, 13
          19,16, ## bld 22, 14
          24,20, ## bnp 22, 13
          35,26, ## btb 21, 13 ## btb 22 is funny
          44,41, ## cp 20, 12
          57,54, ## gnp 21, 13
          66,60, ## hnv 22, 13
          75,73, ## mr 20, 12
          81,78, ## psp 21, 13
          87,83, ## rnv 22, 13
          97,91, ## ski 22, 13
          106,101),## usl 22,13
        ncol=2,byrow=TRUE)

        N<-13
        L<-length(keep)
        dp<-matrix(NA,nrow=N,ncol=L)
        arcdp<-matrix(NA,nrow=N,ncol=L)
        for(i in 1:N){
                dp[i,]<-p[keep,combos[i,1]]-p[keep,combos[i,2]]
                arcdp[i,]<-asin(sqrt(p[keep,combos[i,1]]))-asin(sqrt(p[keep,combos[i,2]]))

        }

        cmat<-matrix(NA,nrow=N,ncol=N)
        for(i in 1:N){ for(j in 1:N){
                if(i==j){
                        cmat[i,j]<-mean(abs(dp[i,]))
                } else if(i < j){
                        cmat[i,j]<-cor(dp[i,],dp[j,])
                } else{
                        cmat[i,j]<-cor(arcdp[i,],arcdp[j,])
                }
        }}
        colnames(cmat)<-c("ABM","BCR","BLD","BNP","BTB","CP","GNP","HNV","MR","PSP","RNV","SKI","USL")
        rownames(cmat)<-c("ABM","BCR","BLD","BNP","BTB","CP","GNP","HNV","MR","PSP","RNV","SKI","USL")
        o<-cor.test(d[upper.tri(d)],cmat[upper.tri(cmat)])
        c_dp[x,1]<-o$estimate
        c_dp[x,2]<-o$p.value

        o<-cor.test(d[lower.tri(d)],cmat[lower.tri(cmat)])
        c_arcdp[x,1]<-o$estimate
        c_arcdp[x,2]<-o$p.value

        cs<-c("darkgray","firebrick")
        
        plot(d[upper.tri(d)],cmat[upper.tri(cmat)],pch=19,col=cs[1+d_tax[upper.tri(d_tax)]],xlab="Geographic distance",ylab="Correlation")
        o<-lm(cmat[upper.tri(cmat)]~d[upper.tri(d)])
        abline(o$coefficients)
        abline(h=0,lty=3)
        title(main=paste("Chromosome",chrom[x],sep=" "))

        plot(d[lower.tri(d)],cmat[lower.tri(cmat)],pch=19,col=cs[1+d_tax[lower.tri(d_tax)]],xlab="Geographic distance",ylab="Correlation")
        o<-lm(cmat[lower.tri(cmat)]~d[lower.tri(d)])
        abline(o$coefficients)
        abline(h=0,lty=3)
}
dev.off()  

########### scaling btb
btb<-25:36
yrs<-c(10,13,13:17,17:19,21,22)
N<-11
l_dp_scale<-vector("list",N)
l_arcdp_scale<-vector("list",N)
pdf("Scaling.pdf",width=8,height=8)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(x in 1:24){
        a1<-as.matrix(fread(a1f[x],header=F))
        a2<-as.matrix(fread(a2f[x],header=F))
        n<-a1+a2
        p<-a2/(a1+a2) ## non-ref
        p[is.na(p)]<-0.001
        colnames(p)<-ids[,1]

        k<-1
        dp_scale<-matrix(NA,nrow=(N*(N-1))/2,ncol=2)
        arcdp_scale<-matrix(NA,nrow=(N*(N-1))/2,ncol=2)
        for(i in 1:(N-1)){for(j in (i+1):N){
                a<-btb[i];b<-btb[j]
                keep<-which(n[,i] > 40 & n[,j] > 40)
                dp<-p[keep,b]-p[keep,a]
                dp_arc<-asin(sqrt(p[keep,b]))-asin(sqrt(p[keep,a]))
                dp_scale[k,2]<-mean(abs(dp))
                arcdp_scale[k,2]<-mean(abs(dp_arc))
                dp_scale[k,1]<-yrs[j]-yrs[i]
                arcdp_scale[k,1]<-yrs[j]-yrs[i]
                k<-k+1
        }}
        l_dp_scale[[x]]<-dp_scale
        l_arcdp_scale[[x]]<-arcdp_scale
        drop<-which(is.na(dp_scale[,1]) | dp_scale[,1]==0)
        xx<-log10(dp_scale[-drop,1])
        yy<-log10(dp_scale[-drop,2]/dp_scale[-drop,1])
        plot(xx,yy,xlab="log10 interval",ylab="log10 rate",pch=19)
        o<-lm(yy~xx)
        abline(o$coefficients)
        mtext(round(o$coefficients[2],3),side=3,line=-2,adj=.5)
        title(paste("chromosome",chrom[x]))
}

dev.off()


save(list=ls(),file="prelim.rdat") 

