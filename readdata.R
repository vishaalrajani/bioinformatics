#Updated: March 1, 2021, Vishaal Rajani
#This is an analysis file for Zia's data
#data is located in /benoukraf_lab/large_memory_share/proj/wgbs_rat/methRaw/
#this Rscript is saved in /research/project/shared/benoukraf_lab/users/vishaal

#format of data.frames are
#     chrBase  chr  base strand coverage freqC freqT
#1 chr1.25761 chr1 25761      F        7 71.43 28.57

#Explanation of groups:
#Group 1 - Early stress (S1-S6)
#Group 2 - Early enrichment(S7 - S12)
#Group 3 - Control E14 (S13 - S18)
#Group 4 - Control GFP (S19 - S24)
#Group 5 - Late stress (S25 - S30)
#Group 6 - Late enrichment (S31 - S36)
#NB: Females and males are first 3 and last 3 individuals of each group, respectively.

#~~~~~~~~~~~~~~~~~~~Parameters~~~~~~~~~~~~~~~~~~~~~~~~~~
depth=20 # indicate read depth 
bin=5 #optional input, havent included it below



#~~~~~~~~~~~~~~~~~~~Reading Data~~~~~~~~~~~~~~~~~~~~~~~~~~
#the following segment reads the data in a for loop, reading ONLY coverage and freqC columns
for(i in 1:36){
print(i)
x=read.table(paste('/benoukraf_lab/large_memory_share/proj/wgbs_rat/methRaw/S',i,'_CpG.txt', sep=""), colClasses = c(rep("NULL",4), rep("numeric",2), rep("NULL",1)), header=T, sep='\t')
xd = x[x$coverage >= depth,]
assign(paste('x',i,sep=""), hist(xd$freqC)$density) #Need to add bin into this histogram function if i want to change the bin size. For now, working with default
assign(paste('d',i, sep=""), density(xd$freqC))
}

#save values x1-36 which are the densities given by the 'hist' function for scatterplotting
#save the densities for each freqC distribution given by the 'density' function
#NB: I realize this is cumbersome, but after frustration i decided to just do this and save time so i dont have to re-read the data 
save(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36, file="rawfreq20.RData")
save(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d30,d31,d32,d33,d34,d35,d36, file="densities20.RData")



#~~~~~~~~~~~~~~~~~~~Scatterplots~~~~~~~~~~~~~~~~~~~~~~~~~~
#NB: can work off saved data from here forward, and commment the first 'read section if not using
load("rawfreq20.RData")

#stitch together densities for scatterplots, based on male/female groupings 1-6
#also cumbersome, but only have to do it once, and its fast.
female1= rbind(x1,x2,x3)
male1 = rbind(x4,x5,x6)
female2=rbind(x7,x8,x9)
male2=rbind(x10,x11,x12)
female3=rbind(x13,x14,x15)
male3=rbind(x16,x17,x18)
female4=rbind(x19,x20,x21)
male4=rbind(x22,x23,x24)
female5=rbind(x25,x26,x27)
male5=rbind(x28,x29,x30)
female6=rbind(x31,x32,x33)
male6=rbind(x34,x35,x36)


#calculate averages and standard deviations for each male/female grouping
for(j in 1:6){ #for each group, this will generate a vector for average and stdev for plotting
avgm=rep(NA,binnum) # clear/generate avg and stdev vectors
avgf=rep(NA,binnum)
sdm=rep(NA,binnum)
sdf=rep(NA,binnum)
assign('male', get(paste('male',j,sep=""))) #create a male column and female column containing the values that I want to use for this iteration
assign('female', get(paste('female',j,sep=""))) 

for(i in 1:binnum){
avgm[i]=mean(male[,i])#calculate average and standard deviation for each column.
avgf[i]=mean(female[,i])
sdm[i]=sd(male[,i])
sdf[i]=sd(female[,i])
}
#now for the j'th group, I have an avgm, avgf, sdm, and sdf vector, containing the average male, average female, standard devation of male, and standard devation of female
#here i save them for future use, as avgm'j',avgf'j',sdm'j', and sdf'j'
assign(paste('avgm', j, sep=""), avgm)
assign(paste('avgf',j, sep=""), avgf)
assign(paste('sdm', j, sep=""), sdm)
assign(paste('sdf',j, sep=""), sdf)

#Here i put everything in data frames for easier plotting later.
#Realizing this now, perhaps easier function to calculate the mean within a dataframe than doing a for loop!
df=data.frame(X = avgm, errX=sdm, Y=avgf, errY=sdf)
assign(paste('df',j,sep=""), df)
assign(paste('plot',j, sep=""), ggplot(data=df, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX)))
}

#plot each scatter. Tried to put this in the loop but the command dev.off() didnt work in the for loop..... so i just did it here manually. Now this can be updated using the par() command... will update. 
pdf('mfscatterplotgroup1.pdf')
ggplot(data=df1, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
dev.off()
pdf('mfscatterplotgroup2.pdf')
ggplot(data=df2, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
dev.off()
pdf('mfscatterplotgroup3.pdf')
ggplot(data=df3, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
dev.off()
pdf('mfscatterplotgroup4.pdf')
ggplot(data=df4, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
dev.off()
pdf('mfscatterplotgroup5.pdf')
ggplot(data=df5, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
dev.off()
pdf('mfscatterplotgroup6.pdf')
ggplot(data=df6, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
dev.off()



######parking lot 
#extra code that I didn't end up using, but maybe useful later

#
# df=data.frame(X = avgm1, errX=sdm1, Y=avgf1, errY=sdf1)
#assign(paste('df',1,sep=""), data.frame(X = avgm1, errX=sdm1, Y=avgf1, errY=sdf1))
#ggplot(data=df, aes(x=X, y=Y)) + geom_point() +
#+ geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) +
#+ geom_errorbar(aes(xmin=X-errX, xmax = X+errX))

#pdf(paste('mfscatterplotgroup',j,'.pdf', sep=""))
#print(df)
#ggplot(data=df, aes(x=X, y=Y)) + geom_point() + geom_errorbar(aes(ymin=Y-errY, ymax = Y+errY)) + geom_errorbar(aes(xmin=X-errX, xmax = X+errX))
#dev.off()
