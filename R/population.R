population <-
function(max.repeat, max.time, Na, N.loci, N.allele, n.recruit.fem, skew.recruit.fem, skew.recruit.mal, surv.ad) {


# expected generation length at start point (approx.)
s <- surv.ad
glength <- sum(s^c(0:100))

whole.rep <- as.data.frame(matrix(ncol=max.repeat*3,nrow=max.time))

## starting the repeat loop
for (r in 1:max.repeat) {


## initializing the population

# adult matrix

adults <- as.data.frame(matrix(NA, ncol=2*N.loci+2, nrow=Na))
seq.loc <- rep(1:N.loci, each = 2)
seq.all <- rep(c(1,2),N.loci)
print.loc <- paste("l",seq.loc,"all",seq.all, sep=".")
colnames (adults) <- c("sex", "age", print.loc)
rownames (adults) <- NULL

adults[,1] <- rep(c("M","F"),each=Na/2)
adults[,2] <- round(rtnorm(Na,1,glength,1,Inf),0)# adults are at least at age 1, sd = 1 generation
adults[,c(3:dim(adults)[2])] <- round(rtnorm((dim(adults)[2]-2)*Na,sum(1:N.allele)/N.allele,N.allele/5,1,N.allele),0)


# define single parameters and matices

n.adults <-  dim(adults)[1]
genotypes.adults <- adults[,3:dim(adults)[2]]
Het <- 1-sum(genotypes.adults[,seq(1,dim(genotypes.adults)[2],2)] == genotypes.adults[,seq(2,dim(genotypes.adults)[2],2)])/ (n.adults*N.loci)# whole Het
Gen.length <- round(mean(adults$age))# generation length
Nind <- n.adults# adult population size


### running the population

for (i in 2:max.time)

{

# sexes
females <- which(adults$sex=="F")
males <- which(adults$sex=="M")


## genotypes of the offspring

# offspring can have different fathers, random mating is assumed

# n.recruit offspring per female
# round(rtnorm(#juv (=mean*# females), mean, skew (=interindividual variance),1,Inf) for females, the reproductive skew is not that strong, most females get at least one offspring
dist.fem <- round(rtnorm(length(females),n.recruit.fem,sqrt(skew.recruit.fem),1,Inf))
genotypes.fem <- genotypes.adults[c(rep(females,dist.fem)),]

# n.recruit offspring per female per year including reproductive skew
n.juv <- sum(dist.fem)
juv <- data.frame(sex=rep("J",n.juv), age=0) # sex is determined later, age=0 per def.


# distribution of males: left skewed, between equal dist. and few males get most and many males get no offspring (rtnorm[0,Inf])
# using the probability as the reproductive skew

# males: sample males of data set randomly, amount=as much as offspring produced, with replace = true because of random mating

dist.mal <- round(rtnorm(length(males),n.juv/length(males),sqrt(skew.recruit.mal),0,Inf))
prob.dist.mal <- dist.mal/sum(dist.mal)
genotypes.mal <- genotypes.adults[males[sample(1:length(males),n.juv,replace=T,prob=prob.dist.mal)],]


# select alleles from females and build the genotypes of the offspring using bits

sel.alleles.fem <- matrix(sample(c(TRUE,FALSE),n.juv*N.loci, replace=T), ncol=N.loci)# TRUE/FALSE matrix Njuv ro x #loci col
sel.alleles.fem <- cbind(sel.alleles.fem,!sel.alleles.fem)# sel.allels and the opposite
sel.alleles.fem <- sel.alleles.fem[,rep(1:N.loci, each=2)+c(0,N.loci)]# sort cols by sel.alleles1,!sel.alleles1,sel.alleles2,!sel.alleles2...
sel.alleles.fem <- as.vector(t(sel.alleles.fem))# t() to convert the matrix, as.vector is writing columnwise
sel.alleles.fem <- as.bit(sel.alleles.fem)

genotypes.juv <-as.matrix(t(genotypes.fem))
gen.juv.vec <- as.vector(genotypes.juv)


# by now, the offspring have genotypes only from the females. We take that allele from the female, where the T/F Matrix is true. The other one has to be raplaced by an allele from the father. 
# That can be either the first or the second one. Therefore, I create a new T/F matrix for males.

sel.alleles.mal <- matrix(sample(c(TRUE,FALSE),n.juv*N.loci, replace=T), ncol=N.loci)# TRUE/FALSE matrix Njuv ro x #loci col
sel.alleles.mal <- cbind(sel.alleles.mal,!sel.alleles.mal)# sel.allels and the opposite
sel.alleles.mal <- sel.alleles.mal[,rep(1:N.loci, each=2)+c(0,N.loci)]# sort cols by sel.alleles1,!sel.alleles1,sel.alleles2,!sel.alleles2...
sel.alleles.mal <- as.vector(t(sel.alleles.mal))# t() to convert the matrix, as.vector is writing columnwise
sel.alleles.mal <- as.bit(sel.alleles.mal)

gen.mal.vec <- as.vector(as.matrix(t(genotypes.mal)))
gen.mal.sel <- gen.mal.vec[as.which(sel.alleles.mal)]

# replace always one allele by the one from the father
gen.juv.vec[as.which(!sel.alleles.fem)] <- gen.mal.sel

genotypes.juv <- as.data.frame(t(matrix(gen.juv.vec,nrow=2*N.loci)))
colnames(genotypes.juv) <- print.loc
juv <- cbind(juv, genotypes.juv)


## survival to the next year and status change

# after reproduction, some juv and adults will survive, others wont, whereby survival will be set as constant over the years

# new adults (seperatet by sex, because of the assumption of a constant sex ratio)
ad.mal <- adults[males[sample(1:length(males),(surv.ad*n.adults)/2)],]# surv.ad/2 of the adult males from the previous year
ad.fem <- adults[females[sample(1:length(females),(surv.ad*n.adults)/2)],]# surv.ad/2 of the adult females from the previous year
adults <- rbind(ad.mal,ad.fem)

# survived juveniles
n.juv <- Na-dim(adults)[1]
juv.sample <- juv[sample(1:dim(juv)[1],n.juv),]
juv.sample$sex <- rep(c("M","F"),n.juv/2)

# change status, former juveniles become adults, adults age 1 year
adults <- rbind(adults,juv.sample)# new adults = survived adults + survived juv from last year
n.adults <-  dim(adults)[1]# new # adults
adults$age <- adults$age+1# change status. age+1
rownames(adults) <- NULL


## calculating Heterozygosity
genotypes.adults <- adults[,3:dim(adults)[2]]

Het[i] <- 1-sum(genotypes.adults[,seq(1,dim(genotypes.adults)[2],2)] == genotypes.adults[,seq(2,dim(genotypes.adults)[2],2)])/(n.adults*N.loci)# whole Het
Gen.length[i] <- round(mean(adults$age))# observed generation length (control)
Nind[i] <- n.adults# adult population size (control)

}

# writing the files for each repeat
mat <- cbind(Het,Gen.length,Nind)

name.seq <- c(1:max.repeat)
name.mat <- paste(paste("mat",name.seq,"s",s,"b",n.recruit.fem, sep="_"),"txt",sep=".")
#write.table(mat,file = name.mat[r], col.names=TRUE, row.names=TRUE)

fill.whole <- seq(r,dim(whole.rep)[2],max.repeat)
whole.rep [,fill.whole] <- mat
colnames(whole.rep)[fill.whole] <- paste(colnames(mat),name.seq[r],sep="_")

message <- paste("I've done the repeat number",r,sep=" ","\n")
cat(message)


}

# extract the data to analyze
Het_pop <- whole.rep[,1:max.repeat]
Gen_pop <- whole.rep[,c((max.repeat+1):(2*max.repeat))]

## temporal reference of one year
mean.het.y <- rowMeans(Het_pop)
years <- c(1:max.time)
lm.het.y <- lm(log(mean.het.y[mean.het.y>0])~years[mean.het.y>0])
slope.y <- coefficients(lm.het.y)[2]
# calculate Ny from the simulation
Ny_sim <- -1/(2*slope.y)
# calculate Ny from the formula
Ny_cal <- (-n.recruit.fem^2*Na)/(1-n.recruit.fem^2+(-1.33+0.44*log(n.recruit.fem))*n.recruit.fem*surv.ad+(1+n.recruit.fem^2-0.001*n.recruit.fem^3)*surv.ad^2)

## temporal reference one generation
mean.gen <- rowMeans(Gen_pop)
mean.het.g <- mean.het.y[seq(1,max.time,round(mean(mean.gen)))]
years.g <- 1:length(mean.het.g)
gen.length <- mean(mean.gen)
lm.het.g <- lm(log(mean.het.g[mean.het.g>0])~years.g[mean.het.g>0])
slope.g <- coefficients(lm.het.g)[2]
# calculate Ne from the simulation
Ne_sim <- -1/(2*slope.g)
# calculate Ne from the formula
Ne_cal <- (-n.recruit.fem^2*Na)/((1-n.recruit.fem^2+(-1.33+0.44*log(n.recruit.fem))*n.recruit.fem*surv.ad+(1+n.recruit.fem^2-0.001*n.recruit.fem^3)*surv.ad^2)*(gen.length))


par(mfcol=c(2,2),mar = c(5,5,4,2) + 0.1)

# original heterozygosity data
plot(mean.het.y,
xlab="",
ylab="mean Heterozygosity",
type="b",pch=1,lwd=2,
cex.main=1.5,cex.lab=1.5,cex.axis=1.5)

# ln(heterozygosity) and the best fit of the linear model
plot(log(mean.het.y),
xlab="time in years",
ylab="ln(mean Heterozygosity)",
type="b",pch=1,lwd=2,
cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
lines(coefficients(lm.het.y)[1] + coefficients(lm.het.y)[2]*years, lwd=2, col="red")

# original heterozygosity per generation
plot(mean.het.g,
xlab="",
ylab="",
type="b",pch=1,lwd=2,
cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
# ln(heterozygosity) and the best fit of the linear model per generation
plot(log(mean.het.g),
xlab="time in generations",
ylab="",
type="b",pch=1,lwd=2,
cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
lines(coefficients(lm.het.g)[1] + coefficients(lm.het.g)[2]*years.g, lwd=2, col="red")

## results table
res_table <- matrix(ncol=2,nrow=7)
colnames(res_table) <- c("parameter","value")
res_table[,1] <- c("slope of heterozygosity loss per year","slope of heterozygosity loss per generation","mean generation length", "Ny[simulation]", "Ny[calc]", "Ne[simulation]","Ne[calc]")
res_table[,2] <- c(round(slope.y,4), round(slope.g,4), round(gen.length,4), round(Ny_sim,4), round(Ny_cal,4), round(Ne_sim,4), round(Ne_cal,4))
res_table <- as.data.frame(res_table)

message3 <- paste ("If you want to analyze another population, recall population(#repeats, time, Na, #loci, #alleles, bf, var(bf),var(bm), sa).\n",
"\n",
"\n",
"NEff 2013",
"\n",
"written by Annegret Grimm :-)\n",
"\n",sep="")

message2 <- paste ("Please switch to the Graphics window to have a look at the heterozygozity behavior over time.\n")
message1 <- paste ("Summary table of an optimal behaviour of your population:\n")

cat(paste("\n","\n"))
cat (message1)
print(res_table)
cat(paste("\n","\n"))
cat (message2)
cat(paste("\n","\n"))
cat (message3)
}
