data = read.csv("star_data.csv")
### get covariates and outcomes we want
var_name = c("GPA_year1","control","age","dad_edn","mom_edn", "sfsp","ssp","sfp", "female")
data1 = data[,var_name]
### delete units with missing value with respect to first year GPA
ind1 = which(is.na(data1$GPA_year1)==0)
data2 = data1[ind1,]
dim(data2)[1]
sum(data2$sfsp)

###############################
var_name = c("GPA_year1","GPA_year2","control", "sfsp", "ssp","sfp", "female")
data12 = data[,var_name] # data with GPA year 1 and year 2
### delete units with missing value with respect to any row
ind2 = which(rowSums(is.na(data12))==0)
data3 = data[ind2,var_name]
dim(data3)[1]
1656 - dim(data3)[1]
sum(data3$sfsp)


###############################
y1=data2$GPA_year1
y2=data3$GPA_year2

# Analysis of GPA_year1 (y1)

N=length(y1)
w=numeric(N)
w[data2$control==1]=1 # ssp=0, sfp=0
w[data2$sfp==1]=2 # ssp=0, sfp=1
w[data2$ssp==1]=3 # ssp=1, sfp=0
w[data2$sfsp==1]=4 # ssp=1, sfp=1


# Display of data

library(tidyverse)
StardataY1Ch6 <- data.frame("w"=w,"y1"=y1)
dim(StardataY1Ch6)
StardataY1Ch6[1:20,]

library(tidyverse)
N_j <- StardataY1Ch6 %>% count(w) %>% .$n %>% as.integer()
meanvec = c(mean(y1[w==1]),mean(y1[w==2]),mean(y1[w==3]),mean(y1[w==4]))
varvec = c(var(y1[w==1]),var(y1[w==2]),var(y1[w==3]),var(y1[w==4]))
rbind(N_j, meanvec, varvec) # TABLE 6.9


boxplot(y1 ~ w, xlab="Treatment combination", ylab=expression(Y[1])) # FIGURE

# ANOVA using lm()
n1=sum(w==1)
n2=sum(w==2)
n3=sum(w==3)
n4=sum(w==4)
anova(lm(y1 ~ factor(w)))
MSRes = ((n1-1)*varvec[1]+(n2-1)*varvec[2]+(n3-1)*varvec[3]+(n4-1)*varvec[4])/(N-4)
VNey=MSRes*(1/n1 + 1/n2 + 1/n3 + 1/n4)/4
SE = sqrt(VNey)

########################################################
#dev.new(width=15)

par(mfrow=c(1,3))

## ME SSP plot
ssp0=mean(c(y1[w==1],y1[w==2]))
ssp1=mean(c(y1[w==3],y1[w==4]))
plot(c(0,1),c(ssp0,ssp1),type="b",main=NULL,axes=FALSE,xlab="SSP",ylab="AVG GPA",ylim=c(1.75,1.90),pch=16, cex.lab=1.5, cex.axis=1.2)
axis(side=1, at=c(0:1))
axis(side=2, at=seq(1.75, 1.85, by=.01))


## ME SFP plot
sfp0=mean(c(y1[w==1],y1[w==3]))
sfp1=mean(c(y1[w==2],y1[w==4]))
plot(c(0,1),c(sfp0,sfp1),type="b",main=NULL,axes=FALSE,xlab="SFP",ylab="AVG GPA",ylim=c(1.75,1.90),pch=16, cex.lab=1.5, cex.axis=1.2)
axis(side=1, at=c(0:1))
axis(side=2, at=seq(1.75, 1.85, by=.01))


## Interaction plot
plot(c(0,1),c(meanvec[1],meanvec[2]),type="b",main=NULL,axes=FALSE,xlab="SFP",ylab="AVG GPA",ylim=c(1.75,1.90),pch=16, cex.lab=1.5, cex.axis=1.2)
points(c(0,1),c(meanvec[3],meanvec[4]),type="b",pch=16)
axis(side=1, at=c(0:1))
axis(side=2, at=seq(1.75, 1.85, by=.01))
text(.5,1.77,"SSP=0", cex=1.5)
text(.5,1.82,"SSP=1", cex=1.5)

###################################
## Factorial effects
###################################

library(AlgDesign)
D = gen.factorial(2,2) ## ME contrast columns in reverse order

D[, c(1:2)] = D[,c(2:1)] ## Reverse the order of the columns
colnames(D) = c("F1","F2")
L = model.matrix(~ F1 + F2 + F1*F2,D)

taubarhatadj = (t(L) %*% meanvec)/4
taubarhat = 2*taubarhatadj[2:4] # adjust factorial effects by 2
taubarhat
Tscaled=abs(taubarhat/SE)
pval=2*(1-pnorm(abs(Tscaled)))

cbind(taubarhat, rep(SE,3), Tscaled, pval)

#################################
## ANOVA from definition

(N_j) # vector of N_j's
(N = sum(N_j)) # Total population size
(y1bar = mean(y1))

# Calculation of SSTot
(SSTot = (N-1)*var(y1))
# Calculation of SSTre
SSTre = 0
for (j in 1:4) {
  SSTre = SSTre + N_j[j]*(meanvec[j] - y1bar)^2
}
SSTre

# Calculation of SSRes
(SSRes = SSTot - SSTre)

# Mean squared errors
(MSTre = SSTre/3)
(MSRes = SSRes/(N-4))

# Calculation of F
(FTre = MSTre/MSRes)


#####################################
# FISHER RANDOMIZATION TEST

library(AlgDesign)
D = gen.factorial(2,2) ## ME contrast columns in reverse order
D[, c(1:2)] = D[,c(2:1)] ## Reverse the order of the columns
colnames(D) = c("F1","F2")
L = model.matrix(~ F1 + F2 + F1*F2,D)
ITER=10000
Fobs=0.1954
Farray=NULL
set.seed(100)
for (count in 1:ITER) {
wnew = sample(c(rep(1,N_j[1]),rep(2,N_j[2]),
              rep(3,N_j[3]),rep(4,N_j[4])))

Fnew=anova(lm(y1 ~ factor(wnew)))$F[1]
Farray = c(Farray, Fnew)
}


hist(Farray, main=NULL, freq=FALSE, xlab="F")
abline(v=Fobs, lty=3, lwd=2)
pvalF = sum(Farray >= Fobs)/ITER
text(4,0.2,paste("p-value= ",round(pvalF,2)))
############################################

############################################
## ANALYSIS FOR FEMALES
############################################

data2f=data2[data2$female==1,]
y1f = data2f$GPA_year1
N=length(y1f)
w=numeric(N)
w[data2f$control==1]=1 # ssp=0, sfp=0
w[data2f$sfp==1]=2 # ssp=0, sfp=1
w[data2f$ssp==1]=3 # ssp=1, sfp=0
w[data2f$sfsp==1]=4 # ssp=1, sfp=1
Nz1z2f = c(sum(data2f$control),sum(data2f$sfp),sum(data2f$ssp),sum(data2f$sfsp))
meanvecf = c(mean(y1f[w==1]),mean(y1f[w==2]),mean(y1f[w==3]),mean(y1f[w==4]))
varvecf = c(var(y1f[w==1]),var(y1f[w==2]),var(y1f[w==3]),var(y1f[w==4]))
rbind(round(Nz1z2f), round(meanvecf,3), round(varvecf,3)) # TABLE 5.8
boxplot(y1f ~ w, xlab="Treatment combination", ylab=expression(Y[1])) # FIGURE

n1=sum(w==1)
n2=sum(w==2)
n3=sum(w==3)
n4=sum(w==4)
anova(lm(y1f ~ factor(w)))
MSRes = ((n1-1)*varvecf[1]+(n2-1)*varvecf[2]+(n3-1)*varvecf[3]+(n4-1)*varvecf[4])/(N-4)
VNeyf=MSRes*(1/n1 + 1/n2 + 1/n3 + 1/n4)/4
SEf = sqrt(VNeyf)

taubarhatadjf = (t(L) %*% meanvecf)/4
taubarhatf = 2*taubarhatadjf[2:4] # adjust factorial effects by 2
taubarhatf
Tscaledf=abs(taubarhatf/SEf)
pvalf=2*(1-pnorm(abs(Tscaledf)))

cbind(taubarhatf, rep(SEf,3), Tscaledf, pvalf)

# FISHER RANDOMIZATION TEST

library(AlgDesign)
D = gen.factorial(2,2) ## ME contrast columns in reverse order
D[, c(1:2)] = D[,c(2:1)] ## Reverse the order of the columns
colnames(D) = c("F1","F2")
L = model.matrix(~ F1 + F2 + F1*F2,D)
ITER=5000
Fobs=0.8495
Farray=NULL
for (count in 1:ITER) {
wnew = sample(c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4)))
Fnew=anova(lm(y1f ~ factor(wnew)))$F[1]
Farray = c(Farray, Fnew)
}

hist(Farray, main=NULL,xlab="F")
abline(v=Fobs, lty=3, lwd=2)
pvalF = sum(Farray >= Fobs)/ITER
text(4,1000,paste("p-value= ",pvalF))
