#setwd("C:/Users/td370/td work/Rutgers 2017 onward/TD 2025 spring/Stat 490/R codes")

PTSDdat <- matrix(scan("ptsd.txt"), ncol=2, byrow=T)
PTSDdat <- data.frame(treatment = PTSDdat[,1], PTSDscore = PTSDdat[,2])

y = PTSDdat$PTSDscore
w = PTSDdat$treatment

anova(lm(y ~ factor(w)))

Fobs <- anova(lm(y ~ factor(w)))[[4]][[1]]

Farray = NULL
for (i in 1:5000){
  wnew = sample(c(rep(1,14),rep(2,10),rep(3,11),rep(4,10)))
  Fnew = anova(lm(y ~ factor(wnew)))[[4]][[1]]
  Farray = c(Farray, Fnew)  
}  

hist(Farray)
abline(v=Fobs, lty=2, lwd=1.5)

pval = sum(Farray > Fobs)/5000


