# Data Analysis from Adsorption Experiments
#Calibration Curve
# Data for patrons must be organized in a .csv file,
  #first row adsorbate weight(g), second column volume of solution, 
  #third column Absorbance (a.u.) or just first column concentration(ppm) and second column (absorbance) 
cal.name <-"kath_calib"
kfile.name <- "kath_K1A_5mg_09-09"
dir.create(kfile.name)
datafile.cal<- read.table(file=paste0("Data/", cal.name,".csv"),
                          header=TRUE,sep=",",stringsAsFactors=FALSE) #calibration curve data
names(datafile.cal)
y.cal <- datafile.cal$absorbance
x.cal <- datafile.cal$concentration
plot(x.cal,y.cal,xlab="concentration(ppm)",ylab="absorption(a.u.)",pch=21,col="blue",bg="red") #plot of absorbance vs concentration
model.cal1 <- lm(y.cal ~x.cal) #regression linear model of absorbance and concentration
sink(file = paste0(kfile.name,"/","calibration.txt"))
print("calibration")
summary(model.cal1)
sink(file = NULL)
intercepto.cal <- coef(model.cal1)[1]
slope.cal <- coef(model.cal1)[2]
abline(model.cal1,col="red") #regression line
#function to calculate interval confidence lines
ci.lines <- function(model){ 
  xm <- sapply(model[[12]][2],mean)
  n <- sapply(model[[12]][2],length)
  ssx <- sum(model[[12]][2]^2)-sum(model[[12]][2])^2/n
  s.t <- qt(0.975,(n-2))
  xv <- seq(min(model[[12]][2]),max(model[[12]][2]),length=100)
  yv <- coef(model)[1]+coef(model)[2]*xv
  se <- sqrt(summary(model)[[6]]^2*(1/n+(xv-xm)^2/ssx))
  ci <- s.t*se
  uyv <- yv+ci
  lyv <- yv-ci
  lines(xv,uyv,lty=2,col="blue")
  lines(xv,lyv,lty=2,col="blue")
}
ci.lines(model.cal1) #plot confidence interval lines
intercepto.cal <- coef(model.cal1)[1]
slope.cal <- coef(model.cal1)[2]

################################################################
############## Function Print results #############
print.results <- function(kfile,kmodel,smodel){
  sink(file = paste0(kfile.name,"/",kfile,".txt"),append = TRUE)
  print("###################")
  print(paste("kinetic model = ",kmodel))
  print(summary(smodel))
  sink(file = NULL)}
###### Function Diagnostic of nonlinear models ##############
diag.model.nl <- function(model){
  plot(fitted(model), residuals(model), xlab="Fitted Values", ylab = "Residuals")
  abline(a=0,b=0)
  model.2 <- lm(qt ~as.factor(time.dat), data =data.f)
  Q <- 2*(logLik(model)-logLik(model.2))
  df.Q <- df.residual(model)-df.residual(model.2) 
  1-pchisq(Q,df.Q)
  standardRes <- residuals(model)/summary(model)$sigma
  qqnorm(standardRes, main="")
  abline(a=0,b=1)
  print(shapiro.test(standardRes))
  plot(profile(model))
}
##############################################################
##################### Data File from kinetic essay ###########
datafile.Kinetic<- read.table(file = paste0("Data/",kfile.name,".csv"),
                              header=TRUE,sep=",",stringsAsFactors=FALSE) #Read the kinetic data file
names(datafile.Kinetic) # Column names in kinetic data file
conc.estimate <- function(y,a,b){# function to estimate concentrations from calibration curve 
  (y-a)/b
}
conc.t <- conc.estimate(datafile.Kinetic$Absorbance,intercepto.cal,slope.cal)
t.infty <- length(datafile.Kinetic$Time) #number of data in data file 
conc.initial <- conc.estimate(datafile.Kinetic$Absorbance[1],intercepto.cal,slope.cal)
conc.e <-conc.estimate(datafile.Kinetic$Absorbance[t.infty],intercepto.cal,slope.cal) 
plot(datafile.Kinetic$Time, conc.t, main = "Kinetic", xlab="time(min)",
     ylab="concentration(ppm)",pch=21,col="blue",bg="red") #Kinetic plot 
ads.ads <- function(x,c0,v,w){ # function to estimate quantity of adsorbate adsorbed by adsorbent mass(mg/g) 
  (c0-x)*(w/v)
}
ms <- datafile.Kinetic$wadsorbate[1]
vol<- datafile.Kinetic$Volume[1]
q <- ads.ads(conc.t,conc.initial,vol,ms) # quantity of adsorbate adsorbed
data.f <- data.frame(time.dat = datafile.Kinetic$Time, ct =conc.t, qt = q)
plot(qt~time.dat, data = data.f, main = "Adsorption", xlab="time(min)",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
qe <- data.f$qt[t.infty]

##############  PFO linear model  #################
kmodel <- "PFO linear model"
interval <- c(qe>data.f$qt)
model.pfo.linear <- lm(log(qe-qt)~time.dat, subset(data.f,qt<qe))
summary(model.pfo.linear)
print.results(kfile.name, kmodel, model.pfo.linear)
qe.pfo.linear <- exp(coef(model.pfo.linear)[1])
k1 <-  coef(model.pfo.linear)[2]
svg(file = paste0(kfile.name,"/PFO linear.svg"))
plot(log(qe-qt)~time.dat, subset(data.f,qt<qe), main = "PFO Linear", 
     xlab="time(min)",ylab="ln(qe-qt)",pch=21,col="blue",bg="red") 
abline(model.pfo.linear,col="green")
dev.off()
##############  PFO nonlinear model  #################
kmodel <- "PFO non linear model"
nlpfoFct <- function(x,beta1,beta2){
  beta1*(1-exp(-beta2*x))
} 
model.pfo.nls <- nls(qt ~ nlpfoFct(time.dat,beta1,beta2), data = data.f, 
                   start = list(beta1 = qe.pfo.linear, beta2 = -k1))
summary(model.pfo.nls)
print.results(kfile.name,kmodel,model.pfo.nls)
unlist(summary(model.pfo.nls)[[1]][1])[1]
svg(file = paste0(kfile.name,"/PFO nonlinear.svg"))
plot(qt~time.dat, data = data.f, main = "PFO non Linear", xlab="time(min)",
     ylab="q(mg/g)",pch=21,col="blue",bg="red")
av <- seq(0,300,0.05)
bv <- predict(model.pfo.nls,list(time.dat=av))
lines(av,bv,col="green")
dev.off()
diag.model.nl(model.pfo.nls)
#anova(model.pfo.nls,model.pfo.2)
#par(mfrow=c(3,1))

############# second order  ##########################
#######  linear
kmodel <- "PSO linear model"
psox.Func <- function(x) {1/x}
model.pso.linear <- lm(1/qt~psox.Func(time.dat), data = data.f, 
                       subset = (time.dat > 2 & time.dat < 300))
summary(model.pso.linear)
print.results(kfile.name,kmodel,model.pso.linear)
qe.pso.linear <- 1/(coef(model.pso.linear)[1])
k2 <- 1/(coef(model.pso.linear)[2]*(qe.pso.linear)^2)
svg(file = paste0(kfile.name,"/PSO linear.svg"))
plot(1/qt~psox.Func(time.dat), data = data.f, 
     subset = (time.dat > 2 & time.dat < 300),
     main = "PSO Linear", xlab="1/time(min)",ylab="1/qt",pch=21,col="blue",bg="red")
abline(model.pso.linear,col="green") 
dev.off()
###############  nonlinear
kmodel <- "PSO non linear model"
nlpsoFct <- function(x,beta1,beta2){
  ((beta1^2)*beta2*x)/(1+beta1*beta2*x)
} 
model.pso.nls <- nls(qt ~ nlpsoFct(time.dat,beta1,beta2), data = data.f, 
                     subset = (time.dat > 2 & time.dat < 300),
                     start = list(beta1 = qe.pso.linear, beta2 = k2))
summary(model.pso.nls)
print.results(kfile.name, kmodel, model.pso.nls)
svg(file = paste0(kfile.name,"/PSO non linear.svg"))
plot(qt ~ time.dat, data = data.f, main= "PSO non Linear", xlab="time(min)",
     ylab="q(mg/g)",pch=21,col="blue",bg="red")
av <- seq(0,300,0.05)
bv <- predict(model.pso.nls,list(time.dat=av))
lines(av,bv,col="green")
dev.off()
diag.model.nl(model.pso.nls)
#################################
###############  Elovich Model ###############
kmodel <- "Elovich model"
elo.Func <- function(x,a,b){
  (1/b)*log(1+a*b*x)
}
plot(qt~time.dat, data = data.f, main = "Elovich non Linear", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
curve(elo.Func(x,a=1.65,b=0.022),col="blue",add = TRUE, lty =2) 
curve(elo.Func(x,a=1.62,b=0.022),col="red",add = TRUE, lty =3)
curve(elo.Func(x,a=70,b=0.10),col="green",add = TRUE, lty =4)
curve(elo.Func(x,a=80,b=0.10),col="orange",add = TRUE, lty =5)# Exploratory plot to find the initial set of parameters 
model.elo.nls <- nls(qt~elo.Func(time.dat,a,b),data=data.f,
                     start=list(a=1.65, b =0.022)) # start parameters from the exploratory plots 
summary(model.elo.nls)
print.results(kfile.name,kmodel,model.elo.nls)
svg(file = paste0(kfile.name,"/Elovich.svg"))
plot(qt~time.dat, data = data.f, main = "Elovich non Linear", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
av.elo <- seq(0,300,0.1)
bv.elo <- predict(model.elo.nls,list(time.dat=av.elo))
lines(av.elo, bv.elo, col="green")
dev.off()
diag.model.nl(model.elo.nls)
#####################################
#The mixer order model (MO model)
kmodel <- "MO model"
library(deSolve)
k1.mo <- coef(model.pfo.nls)[2]
qe.mo = coef(model.pfo.nls)[1]
#k2.mo <- 1
parameters <- c(qe.mo=86.65, k1.mo=0.02, k2.mo =0.000015)
times <- seq(0,300,length=300)
yini=0
derivs <- function(t,y,parms) {
  with(as.list(c(y,parameters)),{
    dy <- k1.mo*(qe.mo-y)+ k2.mo*(qe.mo-y)^2
  result <- dy
  list(result)
  })}
output <- ode(y = yini, time = times, func = derivs, parms = parameters)
head(output)
svg(file = paste0(kfile.name,"/MO.svg"))
plot(qt~time.dat, data = data.f, main = "MO", xlab="time(min)",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
lines(output[,1],output[,2],col="green")
dev.off()
sink(file = paste0(kfile.name,"/", kfile.name,".txt"),append = TRUE)
print("###################")
print(paste("kinetic model = ",kmodel))
print(parameters)
sink(file = NULL)

###################################
############################# External Diffusion Models #################
####### Boyd's external diffusion ###############
#linear approach
kmodel <- "Boyd External Diffusion linear model"
boyd.lin.Func <- function(x,beta1){log(1-(x/beta1))}
model.boyd.linear <- lm(boyd.lin.Func(qt,beta1 = 86.65)~time.dat, 
                        subset(data.f, time.dat > 2 & time.dat < 150))
summary(model.boyd.linear)
print.results(kfile.name,kmodel,model.boyd.linear)
b.boyd <- coef(model.boyd.linear)[2]
svg(file = paste0(kfile.name,"/Boyd's linear.svg"))
plot(boyd.lin.Func(qt,beta1 = 86.65)~time.dat, 
     subset(data.f, time.dat > 2 & time.dat < 150), 
     main = "Boyd's Linear", xlab="time(min))",
     ylab="ln(1-qt/qe)",pch=21,col="blue",bg="red")
abline(model.boyd.linear,col="green") 
dev.off()
#### Boyd's non linear #######
kmodel <- "Boyd External Diffusion non linear model"
boyds.Funct <- function(x,beta1,beta2) {
  beta1*(1-exp(beta2*x))
}
plot(qt~time.dat,data=data.f, main = "Boyd's non Linear", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
#curve(boyds.Funct(x,beta1=60,beta2=-0.025),add=TRUE,lty=2) # Exploratory plot to find the initial set of parameters 
#curve(boyds.Funct(x,beta1=86,beta2=-0.025),add=TRUE,lty=1)
model.boyd.nls <- nls(qt ~ boyds.Funct(time.dat, beta1,beta2),data=data.f,
                      start=list(beta1 = 86, beta2 =-0.025))
summary(model.boyd.nls)
print.results(kfile.name,kmodel,model.boyd.nls)
svg(file = paste0(kfile.name,"/Boyd's non linear.svg"))
plot(qt~time.dat,data=data.f, main = "Boyd's non Linear", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
av.boyd <- seq(0,350,1)
bv.boyd <- predict(model.boyd.nls,list(time.dat=av.boyd))
lines(av.boyd, bv.boyd, col="green")
dev.off()
diag.model.nl(model.boyd.nls)
############### Frusawa & Smith ##################
kmodel <- "Frusawa $ Smith model"
FS.Funct <- function(x,beta1,beta2) {
  w <- ms/vol
  alpha1 <- w*beta1
  alpha2 <- 1/(1+alpha1)
  (conc.initial/w)*(1-(alpha2)-(alpha1*alpha2)*exp((1/(alpha1*alpha2))*beta2*x)) 
}
plot(qt~time.dat,data=data.f, main = "F&S", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
#curve(FS.Funct(x, beta1=5, beta2=-0.01), col="blue", add=TRUE,lty=1)
#curve(FS.Funct(x, beta1=5, beta2=-0.015), col="red", add=TRUE,lty=2)
#curve(FS.Funct(x, beta1=5, beta2=-0.020), col="green", add=TRUE,lty=3)
#curve(FS.Funct(x, beta1=5, beta2=-0.025), col="orange", add=TRUE,lty=4)
# Exploratory plot to find the initial set of parameters 
model.FS.nls <- nls(qt ~ FS.Funct(time.dat, beta1,beta2),data=data.f,
                      start=list(beta1 = 5, beta2 =-0.015))
summary(model.FS.nls)
print.results(kfile.name, kmodel, model.FS.nls)
svg(file = paste0(kfile.name,"/Frusawa & Smith.svg"))
plot(qt~time.dat,data=data.f, main = "F&S", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
av.FS <- seq(0,350,1)
bv.FS <- predict(model.FS.nls,list(time.dat=av.FS))
lines(av.FS, bv.FS, col="green")
dev.off()
diag.model.nl(model.FS.nls)
################# Mathews & Weber #############
kmodel <- "Mathews & Weber model"
MW.Funct <- function(x,beta1) {
  w <- ms/vol
  (conc.initial/w)*(1-exp(-beta1*x)) 
}
plot(qt~time.dat,data=subset(data.f,time.dat<100), main = "M&W", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
#curve(MW.Funct(x,beta1=0.017),col="blue",add=TRUE,lty=2) # Exploratory plot to find the initial set of parameters 
curve(MW.Funct(x,beta1=0.016),col="blue",add=TRUE,lty=1)
model.MW.nls <- nls(qt ~ MW.Funct(time.dat, beta1),data=subset(data.f,time.dat<100),
                    start=list(beta1 = 0.016))
summary(model.MW.nls)
print.results(kfile.name, kmodel, model.MW.nls)
svg(file = paste0(kfile.name,"/Mathews & Weber.svg"))
plot(qt~time.dat,data=subset(data.f,time.dat<100), main = "M&W", xlab="time",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
av.MW <- seq(0,100,1)
bv.MW <- predict(model.MW.nls,list(time.dat=av.MW))
lines(av.MW, bv.MW, col="green")
dev.off()
diag.model.nl(model.MW.nls)
###################################
############################# Internal Diffusion Models #################
#########################################################################
######## Boyd’s intraparticle diffusion model ###########################
kmodel <- "Boyd’s intraparticle diffusion model"
boydIntra.Function <- function(x){-0.4977-log(1-(x/qe))}
newdata <- subset(data.f,(qt > 0 & qt < qe))
Bt <- boydIntra.Function(newdata$qt)
plot(newdata$time.dat,Bt, main="Boyd's intraparticle diffusion model",xlab="time(min)",
     ylab = "Bt",pch=21,col="blue",bg="red")
model.boyd.intra <- lm(Bt~time.dat, data=newdata)
print.results(kfile.name,kmodel,model.boyd.intra)
svg(file = paste0(kfile.name,"/Boyd intraparticle diffusion model.svg"))
plot(newdata$time.dat,Bt, main="Boyd's intraparticle diffusion model",xlab="time(min)",
     ylab = "Bt",pch=21,col="blue",bg="red")
av.boyd.intra <- seq(0,200,1)
bv.boyd.intra <- predict(model.boyd.intra,list(time.dat=av.boyd.intra))
lines(av.boyd.intra, bv.boyd.intra, col="green")
dev.off()
######## Weber & Morris intraparticle diffusion model ###########################
kmodel <- "Weber & Morris intraparticle diffusion model"
WM.Func <- function(x){x^(1/2)}
model.WM <- lm(qt ~ WM.Func(time.dat), data = data.f, subset = (time.dat<10 & time.dat<100))
#output.analysis(model.WM)
print.results(kfile.name,kmodel,model.WM)
svg(file = paste0(kfile.name,"/Weber Morris.svg"))
plot(qt ~ WM.Func(time.dat), data = data.f, main="W & M's intraparticle diffusion model",
     xlab="(time(min))^(1/2)", ylab = "qt(mg/(g L))",pch=21,col="blue",bg="red")
av.WM <- seq(0,200,1)
bv.WM <- predict(model.WM,list(time.dat=av.WM))
lines((av.WM)^(1/2), bv.WM, col="green")
dev.off()
############### Langmuir kinetics model ##########
kmodel <- "Langmuir kinetics model"
c0 <- conc.initial
w <- ms/vol
parameters <- c(ka=0.00035, qe=86.5, w, c0, kd =0.000000011)
times <- seq(0,350,length=350)
yini=0
derivs <- function(t,y,parms) {
 with(as.list(c(y,parameters)),{
  dy <- ka*(c0-y*w)*(qe-y)-kd*y
 result <- dy
  list(result)
 })}
output <- ode(y = yini, time = times, func = derivs, parms = parameters)
svg(file = paste0(kfile.name,"/Langmuir.svg"))
plot(qt~time.dat, data = data.f, main = "Langmuir", xlab="time(min)",ylab="q(mg/g)",
     pch=21,col="blue",bg="red")
lines(output[,1],output[,2],col="green")
dev.off()
sink(file = paste0(kfile.name,"/", kfile.name,".txt"),append = TRUE)
print("###################")
print(paste("kinetic model = ", kmodel))
print(parameters)
sink(file = NULL)