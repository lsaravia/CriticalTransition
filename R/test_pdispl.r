require(poweRlaw)
source("power_fun.r")
#Create distribution object              
# 
m <- displ$new()
m$setXmin(1);m$setPars(2.51)
# 
#Generate random numbers            #
# 
testP <- dist_rand(m, 100)

m$setXmin(30);m$setPars(1.51)

testP <- c(testP,dist_rand(m, 100))

m <- displ$new(testP)

(e1 <- estimate_pars(m))
plot(m,pch=19, col="violet")
m$setPars(e1$pars)
lines(m,col=2)

(e2 = estimate_xmin(m))
m$setPars(e2$pars);m$setXmin(e2$xmin)
lines(m,col=4)


(est<-discpowerexp.fit(testP,1))

x <- sort(unique(testP))
x <- x[x>=est$threshold]
y<-pdiscpowerexp(x,est$exponent,est$rate,est$threshold)

lines(x,y,col=3)

cdfplot_displ_exp(testP,e1$pars,est$exponent,est$rate,1)

# Use contiuous power laws

m <- conpl$new(testP)

(e2 <- estimate_pars(m))
plot(m,pch=19, col="violet")
m$setPars(e1$pars)
lines(m,col=2)

ex <- conexp$new(testP)
(e3 <-estimate_pars(ex))
ex$setPars(e3$pars)
lines(ex,col=3)

(es1<-powerexp.fit(testP,1))
x <- sort(unique(testP))
x <- x[x>=es1$xmin]
y<-ppowerexp(x,es1$xmin,es1$exponent,es1$rate,lower.tail=F)

lines(x,y,col=1)
