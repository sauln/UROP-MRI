
#another interesting example on using a different
#hotelling t2 test
#http://www.uni-kiel.de/psychologie/rexrepos/posts/multHotelling.html




########
#One-sample Hotelling's T^2 test #
#########

#simulating data
set.seed(123)
library(mvtnorm)
Nj    <- c(15, 25)
Sigma <- matrix(c(16,-2, -2,9), byrow=TRUE, ncol=2)
mu1   <- c(-4, 4)
Y1    <- round(rmvnorm(Nj[1], mean=mu1, sigma=Sigma))


#use hotelling t2
muH0 <- c(-1, 2)
library(DescTools)
HotellingsT2Test(Y1, mu=muH0)


#using anova.mlm
Y1ctr  <- sweep(Y1, 2, muH0, "-")
a = (anRes <- anova(lm(Y1ctr ~ 1), test="Hotelling-Lawley"))




##############
#  Hotelling T2 test for two independent samples
##############

mu2 <- c(3, 3)
Y2  <- round(rmvnorm(Nj[2], mean=mu2, sigma=Sigma))
Y12 <- rbind(Y1, Y2)
IV  <- factor(rep(1:2, Nj))

library(DescTools)
a = HotellingsT2Test(Y12 ~ IV)
print(a)


a = anova(lm(Y12 ~ IV), test="Hotelling-Lawley")
print(a)

a = summary(manova(Y12 ~ IV), test="Hotelling-Lawley")
print(a)



##########
# Hotelling's T2-test for two dependent samples
##########

N    <- 20
P    <- 2
muJK <- c(90, 100, 85, 105)
Sig  <- 15
Y1t0 <- rnorm(N, mean=muJK[1], sd=Sig)
Y1t1 <- rnorm(N, mean=muJK[2], sd=Sig)
Y2t0 <- rnorm(N, mean=muJK[3], sd=Sig)
Y2t1 <- rnorm(N, mean=muJK[4], sd=Sig)
Ydf  <- data.frame(id=factor(rep(1:N, times=P)),
                   Y1=c(Y1t0, Y1t1),
                   Y2=c(Y2t0, Y2t1),
                   IV=factor(rep(1:P, each=N), labels=c("t0", "t1")))

dfDiff <- aggregate(cbind(Y1, Y2) ~ id, data=Ydf, FUN=diff)
DVdiff <- data.matrix(dfDiff[ , -1])
muH0   <- c(0, 0)

library(DescTools)
a = HotellingsT2Test(DVdiff, mu=muH0)
print(a)

a = anova(lm(DVdiff ~ 1), test="Hotelling-Lawley")
print(a)

try(detach(package:DescTools))
try(detach(package:mvtnorm))













print(a)




