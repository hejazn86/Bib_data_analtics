

rm(list=ls())

set.seed(2345)

# required libraries
require(ggplot2)
require(ggmap)
require(geoR)
library(rgdal)
require(geostats)
library(gstat)
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


# The data of 1975
data2<-dati1[dati1$anno==1975,c('XUTM', 'YUTM', 'Quota', 'totanno')]
#coordinates(data2) <- ~ XUTM + YUTM
data_2 <- as.geodata(data2, coord.col= 1:2, data.col=4, covar.col = 3)


##visualizing the data
op <- par(no.readonly = TRUE) 
##some graphical adjustments
par(mfrow=c(1,2))
par(mar=c(1,1,2,2)) 
par(mgp = c(2,1,1))

plot(data_2)

par(mfrow=c(1,3))
points(data_2, xlab = "Coord X", ylab = "Coord Y", pt.divide = "rank.prop")
points(data_2, xlab = "Coord X", ylab = "Coord Y", cex.max = 1.7, 
       col = gray(seq(1, 0.1, l=100)), pt.divide = "equal")
points(data_2, pt.divide = "quintile", xlab = "Coord X", ylab = "Coord Y")

par(op)

## Normality test
shapiro.test(data_2$data)

transf<-MASS::boxcox(data2$totanno~ data2$XUTM + data2$YUTM + data2$Quota ,data = data_2)
lambda<-transf$x[which(transf$y==max(transf$y))]    #very close to zero, so log transformation is the best choice 
## the lambda is estimated to be almost 0.62, but we will use lambda=0.5 to transform our data.

## Removing the id form the data of Basilica nad Calabria and transform the coordinates to km
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
pred.grid2 <- data.frame(x=quota$X, y=quota$Y)


# Plotting the empirical variogram - 1975
var2 <- variog(data_2, trend = ~ data2$XUTM + data2$YUTM + data2$Quota, max.dist = 50, lambda=0.5)


par(mfrow=c(1,1))
plot(var2, type="b", main="Empirical Variogram of 1975", col='blue')


var_dir2 <- variog4(data_2, trend = ~ data2$XUTM + data2$YUTM + data2$Quota, max.dist = 50, lambda = 0.5)
par(mfrow=c(1,1))
par(mar=c(2, 2, 3,3))
plot(var_dir2, omnidirectional = T, legend=T)  ## isotropic

eyefit(var2)
#results of eyefit(v2)
##     phi    sigma    nugget(tausq)  max.distance   function 

##    16.02    65.82     8.38         48.07       exponential


## setting the prior and adding nugget
## the values for phi and tausq are always discrete

model2 <- model.control(trend.d= ~ data2$XUTM + data2$YUTM + data2$Quota,
                        trend.l= ~ pred.grid2$x + pred.grid2$y + quota$Quota,
                        cov.m="exponential", lambda=0.5)

prior2 = prior.control(phi.disc=seq(0, 150, l=50),
                       phi.prior="uniform",
                       tausq.rel.discrete = seq(0,1,l=50),
                       tausq.rel.prior = "uniform")

## discrete priors for the range phi and the tau_rel parameters

bayes_2 <- krige.bayes(data2, coords= cbind(data2$XUTM, data2$YUTM), data=data2$totanno,
                       prior = prior2, model = model2,
                       output = output.control(1000,simulations.predictive = T, quantile = c(0.025,0.975))) 

# Plot histograms with smples from the posterior
par(mfrow=c(4,1))
par(mar=c(1,3,1,2))
hist(bayes_2)

par(mfrow=c(2,1))
par(mar=c(2,2,1,2.5))
plot(bayes_2)

prior22 = prior.control(phi.disc=seq(30, 100, l=10),
                        phi.prior="uniform", #reciprocal
                        tausq.rel.discrete = seq(0.2,1,l=10),
                        tausq.rel.prior = "uniform")

bayes_22 <- krige.bayes(data2, coords= cbind(data2$XUTM,data2$YUTM), data=data2$totanno,
                        loc = pred.grid2, prior = prior22, model = model2,
                        output = output.control(1000,simulations.predictive = T,  quantile = c(0.025,0.975)))


## map of the observed area
par(mfrow=c(1,1))
surface(mba.surf(cbind(pred.grid2,bayes_22$predictive$mean),100,100)$xyz.est)


#########################################################

##### PART 2 #####
# evaluate estimates precision using credibility intervals.

b_pred75 <- mba.surf(cbind(pred.grid2, bayes_22$predictive$mean.simulations),no.X=100,no.Y=100)$xyz.est     ## estimations of the model
b_low75 <- mba.surf(cbind(pred.grid2,bayes_22$predictive$quantiles.simulations[,1]),no.X=100,no.Y=100)$xyz.est   ## the lower bound of estimations
b_up75 <- mba.surf(cbind(pred.grid2,bayes_22$predictive$quantiles.simulations[,2]),no.X=100,no.Y=100)$xyz.est    ## the upper bound of the model

b_zlim2 <- range(b_low75$z,b_up75$z,na.rm=T)

par(mfrow=c(1,3))
par(mar=c(2, 2, 3, 1))
surface(b_low75,zlim=b_zlim2,main ="Lower Bound")
points(Cal[,c(1,2)],col = "white",pch = 20)
points(Bas[,c(1,2)],col = "white",pch = 20)

surface(b_pred75,zlim=b_zlim2,main ="Estimation")
points(Cal[,c(1,2)],col = "white",pch = 20)
points(Bas[,c(1,2)],col = "white",pch = 20)

surface(b_up75,zlim=b_zlim2,main ="Upper Bound")
points(Cal[,c(1,2)],col = "white",pch = 20)
points(Bas[,c(1,2)],col = "white",pch = 20)
mtext("1975",cex=2, side = 3, line = - 2, outer = TRUE)


#############################################
## part 3 ##
# discuss your results and compare them to those obtained in homework 4

## ML krige using kirge.conv (homework-4 results)

v2<-variog(data2, coords = cbind(data2$XUTM,data2$YUTM), data = data2$totanno,
           trend = ~ data2$Quota + data2$XUTM +data2$YUTM,max.dist = 50, lambda=0.5)


par(mfrow=c(1,1))
plot(v2,type="b", col='blue')

# We estimate by hand the parameters using eyefit():
eyefit(v2)
#results of eyefit(v2)
##     phi    sigma    nugget(tausq)  max.distance   function 

##    16.02    65.82     8.38         48.07       exponential


############  Kriging estimation

# Predict data points on the grid 
kc<-krige.control(trend.d= ~ data2$XUTM + data2$YUTM + data2$Quota,
                  trend.l= ~ pred.grid2$x + pred.grid2$y + quota$Quota,
                  cov.model="exponential",
                  cov.pars=c(65.82,16.02), lambda=0.5)

ml_krige <- krige.conv(data2, coords = cbind(data2$XUTM,data2$YUTM),
                       data = data2$totanno,
                       locations = pred.grid2, krige = kc, output = output.control(1000, quantile = c(0.025,0.975)))


### Comparing the results of Bayesian and ML kriging using the interval score index

# ML Kriging 
ml.pred.75 <- cbind(pred.grid2, ml_krige$mean.simulations)
ml.low.75 <- cbind(pred.grid2, ml_krige$quantiles.simulations[,1])
ml.up.75 <- cbind(pred.grid2, ml_krige$quantiles.simulations[,2])

ml.interval.score75 = (ml.up.75-ml.low.75)+(2/0.05)*(ml.low.75-ml.pred.75)*(ml.pred.75 < ml.low.75)+(2/0.05)*(ml.pred.75-ml.up.75)*(ml.pred.75>ml.up.75)
#ml.interval.score75

# Bayesian Kriging 

b.pred.75 <- cbind(pred.grid2, bayes_22$predictive$mean.simulations)
b.low.75 <- cbind(pred.grid2,bayes_22$predictive$quantiles.simulations[,1])
b.up.75 <- cbind(pred.grid2,bayes_22$predictive$quantiles.simulations[,2])

b.interval.score75 = (b.up.75 - b.low.75)+(2/0.05)*(b.low.75 - b.pred.75) * (b.pred.75 < b.low.75) + (2/0.05) * (b.pred.75 - b.up.75) * (b.pred.75 > b.up.75)

## returning the lowest interval score to decide which model performed better
interval.Score.75 = c(mean(ml.interval.score75[,3]), mean(b.interval.score75[,3]))
interval.Score.75

## the baysian kriging with lower interval score has better performance compared to ML kriging to estimate the rainfal data.

