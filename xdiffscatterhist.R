#The following script generates smoothscatter plots for female male replicates within each experimental group. Data for use in plotting has already been extracted by 'extractdata.R'

load('averages.RData') #average across male and female replicates (freqC) has been stored in averages.RData. Format of filenamess is favg1, mavg1, favg2, mavg2... for male and female averages in group 1 etc.

load('rawmf.RData') #raw freqC columns for male and female freqC in each experimental group. Format of filenames is f1, m1, f2, m2... etc. 

load('mfloc.RData') #load the extracted locations of each mf

N=1 #put the group number here to save time inputting each plot

pdf(paste('smoothscatterhistogram',N,'.pdf', sep="")) #If you are not using Xwindows and want pdf, uncomment this.

x=get(paste('favg',N,sep="")) #average of freqC column for females of group N; ie favgN
y=get(paste('mavg',N,sep="")) #average of freqC column for males of group N; ie mavgN
x1=get(paste('f',N,sep="")) #store raw data for females (x1), ie fN
y1=get(paste('m',N,sep="")) #store raw data for males (y1), ie mN

#For average densities with sd, first determine bin heights for histograms PER REPLICATE, then average the densities to get an average histogram.
histx1=hist(x1[,1],plot=FALSE)
histx2=hist(x1[,2],plot=FALSE)
histx3=hist(x1[,3],plot=FALSE)
histy1=hist(y1[,1],plot=FALSE)
histy2=hist(y1[,2],plot=FALSE)
histy3=hist(y1[,3],plot=FALSE)

#Build matrices with male and female histograms, and applies average
histxcounts=cbind(histx1$counts,histx2$counts,histx3$counts) #
histycounts=cbind(histy1$counts,histy2$counts,histy3$counts)
histxavg=apply(histxcounts,1,mean) #calculate average of histogram
histyavg=apply(histxcounts,1,mean) 
histxsd=apply(histycounts,1, sd)#calculate standard devation for error bars
histysd=apply(histycounts,1, sd)

#Setup plot layout
layout(mat = matrix(c(2, 1, 0, 3), 
                        nrow = 2, 
                        ncol = 2),
       heights = c(0.4, 2),    # Heights of the two rows
       widths = c(2, 0.4))     # Widths of the two columns

#Plot 1: Scatterplot
par(mar = c(5,4,0,0))
smoothScatter(x,y, nrpoints=0, xlab="Females", ylab="Males") #Create Smoothscatter
stats=cor.test(x1,y1,method="pearson") #Pearson correlation coefficient


abline(a=0.25,b=1) #Create ablines (positive and negative)
abline(a=-0.25,b=1)
outliershigh=cbind(x[x-y>0.25],y[x-y>0.25]) #Plot points above abline 0.25
outlierslow=cbind(x[y-x>0.25],y[y-x>0.25])  #Plot points below abline -0.25
outliersall=rbind(outliershigh,outlierslow) #Put together the columns for Outlier densities
points(outlierslow, cex=0.5, pch=".") #Plot the low outliers
points(outliershigh, cex=0.5, pch=".") #Plot the high outliers


# Plot 2: Top boxplot (females)
par(mar = c(0, 4, 2, 0))
plotTopx=max(histxavg+histxsd)
barcentersx = barplot(histxavg,axes=F,space=0,xaxs="i",yaxs="i",col="deepskyblue3", ylim=c(0,plotTopx))
segments(barcentersx, histxavg-histxsd, barcentersx, histxavg+histxsd)
#arrows(barcenters, histxavg-histxsd, barcenters, histxavg+histxsd, lwd=0.5,angle=90, code=3)
title(paste("Sex Differences Group", N, "PCC = ", round(stats$estimate,digits=3), "\n Proportion Outliers: ", length(outliersall),"/",length(x), " = ", round(length(outliersall)/length(x)*100,digits=1),"%")) #Added the title here. Title includes Pearson Correlation Coefficient value, and the proportion of outliers outside the ablines.

# Plot 3: Right boxplot (males) - Note the switch in the x and y to generate a horizontal plot
par(mar = c(5, 0, 0, 0))
plotTopy=max(histyavg+histysd)
barcentersy = barplot(histyavg,axes=F,space=0,horiz=T, xaxs="i",yaxs="i",col="deepskyblue3", xlim=c(0,plotTopy))
segments( histyavg-histysd, barcentersy, histyavg+histysd,barcentersy)

dev.off()




###OUTLIERS HISTOGRAMS (A)####################

N=1

x=get(paste('favg',N,sep="")) #average of freqC column for females of group N; ie favgN
y=get(paste('mavg',N,sep="")) #average of freqC column for males of group N; ie mavgN
x1=get(paste('f',N,sep="")) #store raw data for females (x1), ie fN
y1=get(paste('m',N,sep="")) #store raw data for males (y1), ie mN

#####


#calculate histogram parameters for female (and to generate ERROR BARS)
#collect outliers for replicates (female)
locateoutliers=c(which(x-y>0.25,1),which(y-x>0.25,1))

xoutliersrep1=c(x1[x-y>0.25,1],x1[y-x>0.25,1])
xoutliersrep2=c(x1[x-y>0.25,2],x1[y-x>0.25,2])
xoutliersrep3=c(x1[x-y>0.25,3],x1[y-x>0.25,3])
xoh1 = hist(xoutliersrep1, plot=FALSE)
xoh2 = hist(xoutliersrep2, plot=FALSE)
xoh3 = hist(xoutliersrep3, plot=FALSE)
xohall=cbind(xoh1$counts, xoh2$counts, xoh3$counts)
xohavg=apply(xohall, 1, mean)
xohsd=apply(xohall, 1, sd)
plotTop=max(xohavg+xohsd)

pdf(paste('Average of Histograms (Female)',N,'.pdf', sep=""))
barcentersx = barplot(xohavg,space=0,col="dark violet", main=paste("Group",N), ylim=c(0,plotTop))
segments(barcentersx, xohavg-xohsd, barcentersx, xohavg+xohsd)
dev.off()


#collect outliers for replicates (male)
youtliersrep1=c(y1[x-y>0.25,1],y1[y-x>0.25,1])
youtliersrep2=c(y1[x-y>0.25,2],y1[y-x>0.25,2])
youtliersrep3=c(y1[x-y>0.25,3],y1[y-x>0.25,3])
yoh1 = hist(youtliersrep1, plot=FALSE)
yoh2 = hist(youtliersrep2, plot=FALSE)
yoh3 = hist(youtliersrep3, plot=FALSE)
yohall=cbind(yoh1$counts, yoh2$counts, yoh3$counts)
yohavg=apply(yohall, 1, mean)
yohsd=apply(yohall, 1, sd)
plotTop=max(yohavg+yohsd)

pdf(paste('Average of Histograms (Male)',N,'.pdf', sep=""))
barcentersy = barplot(yohavg,space=0,col="dark green", main=paste("Group",N), ylim=c(0,plotTop))
segments(barcentersy, yohavg-yohsd, barcentersy, yohavg+yohsd)
dev.off()

#Run t-test on significant "differentially methylated sites"
pvalues=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #generate a vector to fill with p-values
for(i in 1:20){
ttest=t.test(xohall[i,],yohall[i,])
pvalues[i]=ttest$p.value
}

asterisks=pvalues<0.05 #check which pvalues are significant
asterisksint=as.integer(asterisks)
w=which(asterisksint>0)+1
upperbreaks=xoh1$breaks[c(w)] #Create a vector indicating upper and lower limits of each significant bar
lowerbreaks=xoh1$breaks[c(w-1)]
SigNum=length(upperbreaks) #number of significant points


#For this next section, I go through each significant bar in the differentially methylated sites, given the upperbreaks and lower breaks limits calculated in the previous segment. Note that if the lower bar is 0 or 1, I must specifically pull out those indicies who's values =0 or =1 (see the if statements below)
for(n in 1:SigNum){
significantrowsx=unique(c(which(xoutliersrep1>lowerbreaks[n] & xoutliersrep1<=upperbreaks[n]),which(xoutliersrep2>lowerbreaks[n] & xoutliersrep2<=upperbreaks[n]),which(xoutliersrep3>lowerbreaks[n] & xoutliersrep3<=upperbreaks[n])))
significantrowsy=unique(c(which(youtliersrep1>lowerbreaks[n] & youtliersrep1<=upperbreaks[n]),which(youtliersrep2>lowerbreaks[n] & youtliersrep2<=upperbreaks[n]),which(youtliersrep3>lowerbreaks[n] & youtliersrep3<=upperbreaks[n])))

if(lowerbreaks[n]==0){
significantrowsx0=unique(c(which(xoutliersrep1==0),which(xoutliersrep2==0),which(xoutliersrep3==0)))
significantrowsy0=unique(c(which(youtliersrep1==0),which(youtliersrep2==0),which(youtliersrep3==0)))
significantrowsx=unique(c(significantrowsx0,significantrowsx))
significantrowsy=unique(c(significantrowsy0,significantrowsy))
}

if(lowerbreaks[n]==1){
significantrowsx1=unique(c(which(xoutliersrep1==1),which(xoutliersrep2==1), which(xoutliersrep3==1)))
significantrowsy1=unique(c(which(youtliersrep1==1),which(youtliersrep2==1), which(youtliersrep3==1)))
significantrowsx=unique(c(significantrowsx1,significantrowsx))
significantrowsy=unique(c(significantrowsy1,significantrowsy))
}

significantrowsn=sort(unique(c(significantrowsx, significantrowsy))) #This pulls out all the significant indices from the outliers that match the significant bar given by 'n'

#now I want the vector to build outside the loop. so if n==1, then we start the first vector 'allsignificantrows'. If n>1, then we begin to build allsignificant rows for each additional value of n. When we leave the loop we have a very large vector with ALL of the significant values.
if(n==1){
allsignificantrows=significantrowsn
}

if(n>1){
allsignificantrows=sort(unique(c(allsignificantrows,significantrowsn)))
}

}


L=locateoutliers[allsignificantrows]
locations=get(paste("X",N,"loc",sep=""))[L,]
sigx=x[L]
sigy=y[L]
df = data.frame(locations,sigx,sigy) #Output dataframe is a collection of all the significant differentially methylated sites, with their location and methylation score. 




