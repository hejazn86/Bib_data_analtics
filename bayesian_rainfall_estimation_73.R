
rm(list=ls())
set.seed(2345)


# required libraries
require(ggplot2)
require(ggmap)
require(geoR)
require(geostats)
require(gstat)
library(dplyr)
library(sp)
library(fields)
library(spBayes)
library(MASS)
library(MBA)
library(R2jags)


# Loading the data
load("C:/Users/39351/Documents/spatial homeworks/datiAll.RData")
load("C:/Users/39351/Documents/spatial homeworks/quota32.RData")


## Preparing the data

# Removing the data of Puglia from the dataset since we are focing only on Basalicata and Calabria regions
datib <- dati[dati$Nome_Regione == "Basilicata", ]
datic <- dati[dati$Nome_Regione == "Calabria", ]
dati1 <- rbind(datib, datic)
dati1 <- dati1[, c(-1, -2)]         # removing the id columns for not offering any extra information
dati1$XUTM <- dati1$XUTM/1000
dati1$YUTM <- dati1$YUTM/1000

# The data of 1973
data1<-dati1[dati1$anno==1973,c('XUTM', 'YUTM', 'Quota', 'totanno')]
#coordinates(data1) <- ~ XUTM + YUTM
data_1 <- as.geodata(data1, coord.col= 1:2, data.col=4, covar.col = 3)


##visualizing the data
op <- par(no.readonly = TRUE) 
##some graphical adjustments
par(mfrow=c(1,2))
par(mar=c(1,1,2,2)) 
par(mgp = c(2,1,1))

plot(data_1)
par(mfrow=c(1, 3))
points(data_1, xlab = "Coord X", ylab = "Coord Y", pt.divide = "rank.prop")
points(data_1, xlab = "Coord X", ylab = "Coord Y", cex.max = 1.7, 
       col = gray(seq(1, 0.1, l=100)), pt.divide = "equal")
points(data_1, pt.divide = "quintile", xlab = "Coord X", ylab = "Coord Y")

par(op)

## Normality test
shapiro.test(data_1$data)

transf<-MASS::boxcox(data1$totanno~ data1$XUTM + data1$YUTM + data1$Quota ,data = data_1)
lambda<-transf$x[which(transf$y==max(transf$y))]    #very close to zero, so log transformation is the best choice 
## lambda is estimated to be approximately 0.38, so we are going to use this to transform our data

## Removing the id form the data of Basalica nad Calabria and transform the coordinates to km
Bas <- Bas[, -1]
Bas$XUTM <- Bas$XUTM/1000
Bas$YUTM <- Bas$YUTM/1000

Cal <- Cal[, -1]
Cal$XUTM <- Cal$XUTM/1000
Cal$YUTM <- Cal$YUTM/1000


## creating a grid for the predictions
quotacal<-Calabriaquota32[,-1]    # Calabria region
quotabas<-Basquota32[,-1]         # Basalica region

# harmonising the column names 
colnames(quotacal)<-colnames(quotabas)
quota<- rbind(quotabas,quotacal)
quota$XUTM<-quota$XUTM/1000
quota$YUTM<-quota$YUTM/1000
#coordinates(quota) <- ~ XUTM +YUTM

pred.grid1 <- data.frame(x=quota$X, y=quota$Y)


# Plotting the empirical variogram - 1973
var1 <- variog(data_1, trend = ~ data1$XUTM + data1$YUTM + data1$Quota, max.dist = 60, lambda = 0.38)

par(mfrow=c(1,1))
par(mar=c(2, 3, 3,3))
plot(var1, type="b", main="Empirical Variogram of 1973", col='blue')


var_dir1 <- variog4(data_1, trend = ~data1$XUTM + data1$YUTM + data1$Quota, lambda = 0.38, max.dist = 60)
par(mfrow=c(1,1))
par(mar=c(2, 2, 2,2))
plot(var_dir1, omnidirectional = T, legend=T) #isotropic


eyefit(var1)
#results of eyefit(var1)
##     phi    sigma    nugget(tausq)  max.distance   function 

##    19.22    14.12     1.77          57.67       exponential


## setting the prior and adding nugget
## the values for phi and tausq are always discrete

model1 <- model.control(trend.d= ~ data1$XUTM + data1$YUTM + data1$Quota,
                        trend.l= ~ pred.grid1$x + pred.grid1$y + quota$Quota,
                        cov.m="exponential", lambda = 0.38)

prior1 = prior.control(phi.disc=seq(0, 150, l=50),
                       phi.prior="uniform", #reciprocal
                       tausq.rel.discrete = seq(0,2,l=50),
                       tausq.rel.prior = "uniform")


bayes_1 <- krige.bayes(data1, coords=cbind(data1$XUTM + data1$YUTM), data=data1$totanno,
                       prior = prior1, model = model1,
                       output = output.control(1000, simulations.predictive = T, quantile = c(0.025,0.975))) 

# Plot histograms with samples from the posterior
par(mfrow=c(4,1))
par(mar=c(1,2,1,2))
hist(bayes_1)

par(mfrow=c(2,1))
par(mar=c(2.5,2.5,1,2))
plot(bayes_1)

prior11 = prior.control(phi.disc=seq(30, 100, l=10),
                       phi.prior="uniform", #reciprocal
                       tausq.rel.discrete = seq(0,1,l=10),
                       tausq.rel.prior = "uniform")

bayes_11 <- krige.bayes(data1, coords=cbind(data1$XUTM, data1$YUTM), data=data1$totanno,
                        loc = pred.grid1, prior = prior11, model = model1,
                        output = output.control(1000,simulations.predictive = T,  quantile = c(0.025,0.975)))



# Maps:

par(mfrow=c(1,1))
surface(mba.surf(cbind(pred.grid1,bayes_11$predictive$mean),100,100)$xyz.est)


###########################################
## Part 2 ##

# evaluate estimates precision using credibility intervals.

b_pred73 <- mba.surf(cbind(pred.grid1, bayes_11$predictive$mean.simulations),no.X=100,no.Y=100)$xyz.est
b_low73 <- mba.surf(cbind(pred.grid1,bayes_11$predictive$quantiles.simulations[,1]),no.X=100,no.Y=100)$xyz.est
b_up73 <- mba.surf(cbind(pred.grid1,bayes_11$predictive$quantiles.simulations[,2]),no.X=100,no.Y=100)$xyz.est

b_zlim <- range(b_low73$z,b_up73$z,na.rm=T)

par(mfrow=c(1,3))
par(mar=c(2, 2, 3, 1))
surface(b_low73,zlim=b_zlim,main ="Lower Bound")
points(Cal[,c(1,2)],col = "gray",pch = 20)
points(Bas[,c(1,2)],col = "gray",pch = 20)

surface(b_pred73,zlim=b_zlim,main ="Estimation")
points(Cal[,c(1,2)],col = "gray",pch = 20)
points(Bas[,c(1,2)],col = "gray",pch = 20)

surface(b_up73,zlim=b_zlim,main ="Upper Bound")
points(Cal[,c(1,2)],col = "gray",pch = 20)
points(Bas[,c(1,2)],col = "gray",pch = 20)
mtext("1973",cex=2, side = 3,line = - 2,outer = TRUE)


#############################################
## part 3 ##
# discuss your results and compare them to those obtained in homework 4

## ML krige using kirge.conv (homework-4 results)

v1<-variog(data1, coords = cbind(data1$XUTM,data1$YUTM), data = data1$totanno,
           trend = ~ data1$Quota + data1$XUTM +data1$YUTM,max.dist = 60, lambda=0.38)


par(mfrow=c(1,1))
plot(v1,type="b", col='blue')

# We estimate by hand the parameters using eyefit():
eyefit(v1)
#results of eyefit(v1)
##     phi    sigma    nugget(tausq)  max.distance   function

##    19.22    14.12      1.77          57.67       exponential


############  Kriging estimation

# Predict data points on the grid 
kc<-krige.control(trend.d= ~ data1$XUTM + data1$YUTM + data1$Quota,
                  trend.l= ~ pred.grid1$x + pred.grid1$y + quota$Quota,
                  cov.model="exponential",
                  cov.pars=c(14.12,19.22), lambda=0.38)


ml_krige <- krige.conv(data1, coords = cbind(data1$XUTM,data1$YUTM),
                       data = data1$totanno,
                       locations = pred.grid1, krige = kc, output = output.control(1000, quantile = c(0.025,0.975)))



### Comparing the results of Bayesian and ML kriging using the interval score index

# ML Kriging 
ml.pred.73 <- cbind(pred.grid1, ml_krige$mean.simulations)
ml.low.73 <- cbind(pred.grid1, ml_krige$quantiles.simulations[,1])
ml.up.73 <- cbind(pred.grid1, ml_krige$quantiles.simulations[,2])

ml.interval.score73 = (ml.up.73-ml.low.73)+(2/0.05)*(ml.low.73-ml.pred.73)*(ml.pred.73<ml.low.73)+(2/0.05)*(ml.pred.73-ml.up.73)*(ml.pred.73>ml.up.73)
#ml.interval.score73

# Bayesian Kriging 

b.pred.73 <- cbind(pred.grid1, bayes_11$predictive$mean.simulations)
b.low.73 <- cbind(pred.grid1,bayes_11$predictive$quantiles.simulations[,1])
b.up.73 <- cbind(pred.grid1,bayes_11$predictive$quantiles.simulations[,2])

b.interval.score73 = (b.up.73 - b.low.73)+(2/0.05)*(b.low.73 - b.pred.73) * (b.pred.73 < b.low.73) + (2/0.05) * (b.pred.73 - b.up.73) * (b.pred.73 > b.up.73)

## returning the lowest interval score to decide which model performed better
interval.Score.73 = c(mean(ml.interval.score73[,3]), mean(b.interval.score73[,3]))
interval.Score.73

## based on the interval score using the function krige.conv for ML kriging, the ML kriging still has the lower interval score compared to Bayesian kriging, so we can assume that ML kriging performed better on estimating the rainfall compared with bayesian kriging.


