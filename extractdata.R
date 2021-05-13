#This data unpacks the results of the UNITE function run by Geoff. 
#Data files are located in /benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/
#Each group (1-6) contains data file 'meth'. The relevant vectors can be called using getData function.

#Relevant data are in the format with columns:
#chr, start, end, strand, coverage1, numCs1, numTs1, coverage2, numCs2, numTs2, coverage3, numCs3, numTs3, coverage4, numCs4, numTs4, coverage5, numCs5, numTs5, coverage6, numCs6, numTs6

print(1) #Since this file takes a long time to run, these print commands add signposts to let me know how the program is running
load("/benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/group1.RData")
meth #This installs some packages so needs to be run before proceeding
X=getData(meth) #Store relevant data vector in X
print(1.5) #signpost
X1loc=data.frame(X$chr, X$start, X$end) #Store group1 locations
f1=cbind(X$numCs1/X$coverage1, X$numCs2/X$coverage2, X$numCs3/X$coverage3) #Calculate proportion methylation for individuals 1,2,3 (females)
m1=cbind(X$numCs4/X$coverage4, X$numCs5/X$coverage5, X$numCs6/X$coverage6) #Calculate proportion methylation for individuals 4,5,6 (males)
fsd1=apply(f1,1, sd)
msd1=apply(m1,1, sd)
favg1=apply(f1, 1, mean)
mavg1=apply(m1, 1, mean)

print(2)
load("/benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/group2.RData")
meth
X=getData(meth)
print(2.5)
X2loc=data.frame(X$chr, X$start, X$end)
f2=cbind(X$numCs1/X$coverage1, X$numCs2/X$coverage2, X$numCs3/X$coverage3)
m2=cbind(X$numCs4/X$coverage4, X$numCs5/X$coverage5, X$numCs6/X$coverage6)
fsd2=apply(f2,1, sd)
msd2=apply(m2,1, sd)
favg2=apply(f2, 1, mean)
mavg2=apply(m2, 1, mean)

print(3)
load("/benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/group3.RData")
meth
X=getData(meth)
print(3.5)
X3loc=data.frame(X$chr, X$start, X$end)
f3=cbind(X$numCs1/X$coverage1, X$numCs2/X$coverage2, X$numCs3/X$coverage3)
m3=cbind(X$numCs4/X$coverage4, X$numCs5/X$coverage5, X$numCs6/X$coverage6)
fsd3=apply(f3,1, sd)
msd3=apply(m3,1, sd)
favg3=apply(f3, 1, mean)
mavg3=apply(m3, 1, mean)

print(4)
load("/benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/group4.RData")
meth
X=getData(meth)
print(4.5)
X4loc=data.frame(X$chr, X$start, X$end)
f4=cbind(X$numCs1/X$coverage1, X$numCs2/X$coverage2, X$numCs3/X$coverage3)
m4=cbind(X$numCs4/X$coverage4, X$numCs5/X$coverage5, X$numCs6/X$coverage6)
fsd4=apply(f4,1, sd)
msd4=apply(m4,1, sd)
favg4=apply(f4, 1, mean)
mavg4=apply(m4, 1, mean)

print(5)
load("/benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/group5.RData")
meth
X=getData(meth)
print(5.5)
X5loc=data.frame(X$chr, X$start, X$end)
f5=cbind(X$numCs1/X$coverage1, X$numCs2/X$coverage2, X$numCs3/X$coverage3)
m5=cbind(X$numCs4/X$coverage4, X$numCs5/X$coverage5, X$numCs6/X$coverage6)
fsd5=apply(f5,1, sd)
msd5=apply(m5,1, sd)
favg5=apply(f5, 1, mean)
mavg5=apply(m5, 1, mean)

print(6)
load("/benoukraf_lab/large_memory_share/proj/wgbs_rat/RObjects/group6.RData")
meth
X=getData(meth)
print(6.5)
X6loc=data.frame(X$chr, X$start, X$end)
f6=cbind(X$numCs1/X$coverage1, X$numCs2/X$coverage2, X$numCs3/X$coverage3)
m6=cbind(X$numCs4/X$coverage4, X$numCs5/X$coverage5, X$numCs6/X$coverage6)
fsd6=apply(f6,1, sd)
msd6=apply(m6,1, sd)
favg6=apply(f6, 1, mean)
mavg6=apply(m6, 1, mean)

save(X1loc, X2loc, X3loc, X4loc, X5loc, X6loc, file="mfloc.RData") # save the locations for each CpG site
save(f1,f2,f3,f4,f5,f6,m1,m2,m3,m4,m5,m6, file="rawmf.RData")
save(fsd1,fsd2,fsd3,fsd4,fsd5,fsd6,msd1,msd2,msd3,msd4,msd5,msd6, file="stdev.RData")
save(favg1,favg2,favg3,favg4,favg5,favg6,mavg1,mavg2,mavg3,mavg4,mavg5,mavg6, file="averages.RData")
