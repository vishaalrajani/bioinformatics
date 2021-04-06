#The following script generates smoothscatter plots for female male replicates within each experimental group. Data for use in plotting has already been extracted by 'readdata.R'

load('averages.RData') #average across male and female replicates (freqC) has been stored in averages.RData. Format of filenamess is favg1, mavg1, favg2, mavg2... for male and female averages in group 1 etc.

load('rawmf.RData') #raw freqC columns for male and female freqC in each experimental group. Format of filenames is f1, m1, f2, m2... etc. 

N=3 #put the group number here to save time inputting each plot

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


pdf(paste('XdiffOutliers',N,'.pdf', sep=""))
plot(density(outliersall[,1]), col="dark violet", xlab="Female/Male FreqC", main=paste('Density of Outliers for Group',N))
lines(density(outliersall[,2]), col="dark green")
dev.off()

#






##########################OLD CODE#####################################################
#pdf('smoothscatterhistogram1.pdf')
#
#x=favg1
#y=mavg1
#
#layout(mat = matrix(c(2, 1, 0, 3), 
#                        nrow = 2, 
#                        ncol = 2),
#       heights = c(0.4, 2),    # Heights of the two rows
#       widths = c(2, 0.4))     # Widths of the two columns
#plot.new()
#
## Plot 1: Scatterplot
#par(mar = c(5, 4, 0, 0))
#smoothScatter(x,y, nrpoints=0, xlab="Females", ylab="Males")
#
## Plot 2: Top (height) boxplot
#par(mar = c(0, 4, 0, 0))
#xhist=hist(x, plot=FALSE)
#barplot(xhist$counts,axes=F,space=0,xaxs="i",yaxs="i",col="deepskyblue3")
#
#
## Plot 3: Right (weight) boxplot
#par(mar = c(5, 0, 0, 0))
#yhist=hist(y,plot=FALSE)
#barplot(yhist$counts,axes=F,space=0,horiz=T,xaxs="i",yaxs="i",col="deepskyblue3")
#
#mtext("Group 1: Females v. Males", side = 3, outer = TRUE)
#
#dev.off()
#
######
#
#require(ggplot2)
# x<-favg1
# y<-mavg1
# xy<-data.frame(x,y)
#     xhist <- qplot(x, geom="histogram") + scale_x_continuous(limits=c(min(x),max(x))) + opts(axis.text.x = theme_blank(), axis.title.x=theme_blank(), axis.ticks = theme_blank(), aspect.ratio = 5/16, axis.text.y = theme_blank(), axis.title.y=theme_blank(), background.colour="white")
#     yhist <- qplot(y, geom="histogram") + coord_flip() + opts(background.fill = "white", background.color ="black")
#
#     yhist <- yhist + scale_x_continuous(limits=c(min(x),max(x))) + opts(axis.text.x = theme_blank(), axis.title.x=theme_blank(), axis.ticks = theme_blank(), aspect.ratio = 16/5, axis.text.y = theme_blank(), axis.title.y=theme_blank() )
#
#
#     scatter <- smoothscatter(x,y, nrpoints=0)  + scale_x_continuous(limits=c(min(x),max(x))) + scale_y_continuous(limits=c(min(y),max(y)))
#none <- smoothscatter(x,y, nrpoints=0) + geom_blank()
#dev.off()
#
#grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
#
#
#
#par(mfrow=c(2,2))
#xhist
#
#plot(empty)
#
#scatter
#
#yhist
#dev.off()
#
#
#pdf('smoothscatterhist1.pdf')
#pdf('smoothscatterhist2.pdf')
#pdf('smoothscatterhist3.pdf')
#pdf('smoothscatterhist4.pdf')
#pdf('smoothscatterhist5.pdf')
#pdf('smoothscatterhist6.pdf')