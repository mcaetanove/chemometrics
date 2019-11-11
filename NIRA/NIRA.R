# ***********NIR Spectra collection***************

#Create a list with the different directories (list.dir) 
list.dir <- c("09-27-2019","10-14-2019","10-15-2019","10-17-2019","10-18-2019","10-21-2019","10-24-2019","10-31-2019","11-05-2019","11-06-2019")
ndir <- length(list.dir)
#Create a list with the different .csv files (temp)
temp <- list()
for (i in 1:ndir) {
  temp <- append(temp, list.files(path = list.dir[i], pattern="*.csv", 
                                  full.names =TRUE))}
#Read and Create a list of data frames with each of the files (five columns: wavelength,%R,honey,adulterant,percentage)
myfiles <- lapply(temp, read.csv, header=TRUE, colClasses=c("integer","numeric","character","character","numeric"),col.names=c("nm","reflectance","honey","adulterant","percentage"))
nfiles <- length(myfiles)
#Loop, create a vector with the information of each spectrum 
#(vhoney: type of honey, type of adulterant and percentage of adulterant), 
#and create a matrix with the spectra trimed within a specified range of wavelengths  
maxlambda <- 2600
minlambda <- 1305
#A vector with the wavelengths evenly spaced 5nm (as was registered in the experiments)
lambda <- seq(maxlambda,minlambda,-5)
vhoney <- matrix(, nrow=nfiles, ncol=3) 
colnames(vhoney) <- c("honey","adulterant","percentage")
mhoney <- matrix(, nrow = nfiles, ncol = length(lambda)) #empty matrix with the number of columns and rows required

for (j in 1:nfiles) {
  x <-length(which(as.numeric(myfiles[[j]]$nm)> maxlambda)) #number of wavelenghts above maxlambda in jth spectrum  
  y <-length(which(as.numeric(myfiles[[j]]$nm) < minlambda)) #number of wavelenghts lower than minlambda in jth spectrum
  z <- nrow(myfiles[[j]]) #number of rows in each file
  vhoney[j,]<-c(myfiles[[j]]$honey[1],myfiles[[j]]$adulterant[1], myfiles[[j]]$percentage[1]) #vector with information of each spectrum
  mhoney[j,] <-t(myfiles[[j]][-c(1:x,(z-y+1):z), ]$reflectance) # matrix with trimed spectra
}

#Calculate absorbances from transflectance data as log(1/R) and store them in a dataframe Absorbances$NIRS
absorbances <- data.frame(vhoney) #create a data.frame and store the vector with information
tmphoney <-mhoney 
tmphoney[tmphoney<0] <- 0 # convert negative reflectances in zeros
absorbances$NIRS <- apply(tmphoney,1:2, function(x) {log10(100/x)}) 

#Plotting spectra of absorbance (honey, adulterant and percentages) 
matplot(lambda, t(absorbances$NIRS), type = "l", lty = 1, xlab = "Wavelength (nm)", 
        ylab = "Absorbances(log(1/R))")

#******************************************* Subsetting Spectra ****************************************
# Subset de spectra corn-sirup
cs <- absorbances[which(absorbances$adulterant == "cs"),]
matplot(lambda, t(cs$NIRS), type = "l", lty = 1, xlab = "Wavelength (nm)", 
        ylab = "Absorbances(log(1/R))")
#Subset spectra cane-sugar
canesugar <- absorbances[which(absorbances$adulterant == "p"),]
matplot(lambda, t(canesugar$NIRS), type = "l", lty = 1, xlab = "Wavelength (nm)", 
        ylab = "Absorbances(log(1/R))")

#Subset honeys, adulterants
ssha <- absorbances[which(c(absorbances$percentage == 0 | absorbances$percentage == 100)),]
honeysid <- levels(absorbances$honey)
colors <- c(rep(c("cyan","red"),5),"cyan","cyan")
colorlabels <- c(rep(1,24))


#**********************************************************Color Function **********************************************
colorfunc <- function(value,i) {
  colorlabels[value] <<- colors[i]
  cat("color",colors[i],"i",i,"\n")
}

for (i in 1:length(honeysid)) {
  temp <- which(ssha$honey == honeysid[i])
  if (length(temp) != 0 )
    lapply(temp, colorfunc, i)
}

#*********************************************** Plotting Spectra *********************************

matplot(lambda, t(ssha$NIRS), type = "l", lty = c(rep(1,10),2,rep(1,8),rep(2,3),1), 
        col=colorlabels, xlab = "Wavelength (nm)", 
        ylab = "Absorbances(log(1/R))")

ss.honey <- absorbances[which(c(absorbances$percentage == 0)),]
honeysid <- levels(absorbances$honey)
colors2 <- c("blue","red","orange","red","cyan","red","brown","red","gray","red","green","purple")
colorlabels2 <- c(rep(1,20))

colorfunc <- function(value,i) {
  colorlabels2[value] <<- colors2[i]
  cat("color",colors[i],"i",i,"\n")
}

for (i in 1:length(honeysid)) {
  temp <- which(ss.honey$honey == honeysid[i])
  if (length(temp) != 0 )
    lapply(temp, colorfunc, i)
}


honeysidA <- levels(absorbances$honey)
matplot(lambda, t(ss.honey$NIRS), type = "l", lty = 1, col = colorlabels2, 
        xlab = "Wavelength (nm)", 
        ylab = "Absorbances(log(1/R))",ylim = (c(0,0.3)))


#**********************************************************SPECTRA PRETREATMENT **************************

#NIR.mc <- scale(cs.ha$NIRS, scale = FALSE)
#TODO: Random Distribution; 


#¨****************************************************** Hierarchical Cluster *************


ahoney <- absorbances[which(absorbances$percentage != 100),]
ahoney <-ahoney[-c(grep("*p",ahoney$honey)),]
set.seed(113)
random.rows <- sample(nrow(ahoney))
ahoney.random <- ahoney[random.rows,]
ahoney.random.ad <- ahoney.random[which(ahoney.random$percentage != 0),]
ahoney.random.ad <- ahoney.random.ad[which(ahoney.random.ad$honey != "5"),]
subset <- sample(nrow(ahoney.random.ad), 20)
ahoney.random.adss <- ahoney.random.ad[subset,]
pure.honeys <- ahoney[which(ahoney$percentage == 0),]
pure.honeys <- pure.honeys[-c(2,4),] # exclude honey1 (85) and honey5 (2)
ahoney.random.adss <- rbind(ahoney.random.adss, pure.honeys)
vahoney <- data.frame(ahoney.random.adss$honey,ahoney.random.adss$adulterant,ahoney.random.adss$percentage)
colors3 <- c()
ahoneyLabels <- c()
apply(vahoney, 1, function(value){
  if (value[3] == 0) { #pure honey
    colors3 <<- append(colors3,"cyan")
    ahoneyLabels <<- append (ahoneyLabels, value[1])
  } else if (value[2] == "cs") { # corn syrup adulterated honey
    colors3 <<- append(colors3,"red")
    ahoneyLabels <<- append (ahoneyLabels, value[2])
    #ahoneyLabels <<- append (ahoneyLabels, paste(value[1], value[2], collapse = "-"))
  } else { # brown cane sugar adulterated honey
    colors3 <<- append(colors3,"green")
    ahoneyLabels <<- append (ahoneyLabels, value[2])
    #ahoneyLabels <<- append (ahoneyLabels, paste(value[1], value[2], collapse = "-"))
  }
})

NIR.mc.ahoney <- scale(ahoney.random.adss$NIRS, scale = FALSE)
honey.dist <- dist(NIR.mc.ahoney)
honey.hcsingle <- hclust(honey.dist, method = "single")
plot(honey.hcsingle, hang = -1, cex=0.8, labels = ahoneyLabels)


#***************PCA*******************

cs.ha <- absorbances#[which(absorbances$adulterant != "cs"|absorbances$percentage==0),]
NIR.mc <- scale(cs.ha$NIRS, scale = FALSE)
NIR.sweep<- sweep(cs.ha$NIRS, MARGIN = 1, apply(cs.ha$NIRS, 1, max),
         FUN = "/")
matplot(lambda, t(NIR.sweep),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "%R (mean-centered)", lty = 1, col = 1)
NIR.sweep.scaled <- scale(NIR.sweep, scale = FALSE)
matplot(lambda, t(NIR.sweep.scaled),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "%R (mean-centered)", lty = 1, col = 1)

matplot(lambda, t(NIR.mc),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "%R (mean-centered)", lty = 1, col = 1)

nir.prcomp <- prcomp(NIR.mc)
summary(nir.prcomp)
plot(nir.prcomp)

nir.prcomp <- prcomp(cs.ha$NIRS)
summary(nir.prcomp)
plot(nir.prcomp)
nir.loadings <- nir.prcomp$rotation[,1:5]
#par(mfrow = c(1,2), pty = "s")
offset <- c(0, 0.09)
plot(nir.loadings[,1:2], type = "l",
     xlim = range(nir.loadings[,1]) + offset,
     xlab = "PC 1 (98.37%)", ylab = "PC 2 (1.75%)")
#par(mfrow = c(1,1))
vhoney2 <- data.frame(cs.ha$honey,cs.ha$adulterant,cs.ha$percentage)
label <- apply(vhoney2, 1, paste, collapse = "-")
biplot(nir.prcomp, xlabs = label, ylabs = lambda, cex=0.6)

## Test custom biplot
library(ggfortify)
colors4 <- c()
apply(vhoney2,1,function(value){
  if (value[3] == 0) { #pure honey
    colors4 <<- append(colors4,1)
  } else if (value[2] == "cs") { # corn syrup adulterated honey
    colors4 <<- append(colors4,2)
  } else { # brown cane sugar adulterated honey
    colors4 <<- append(colors4,3)
  }
})
autoplot(nir.prcomp, loadings=TRUE, loadings.label=TRUE, loading.label.color = 'blue')
## ---

NIR.mc <- scale(honeylog$NIRS, scale = FALSE)
matplot(lambda, t(NIR.mc),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "%R (mean-centered)", lty = 1, col = 1)
nir.prcomp <- prcomp(honeylog$NIRS)
summary(nir.prcomp)
plot(nir.prcomp)
nir.loadings <- nir.prcomp$rotation[,1:4]
#par(mfrow = c(1,2), pty = "s")
offset <- c(0, 0.09)
plot(nir.loadings[,1:2], type = "p",
     xlim = range(nir.loadings[,1]) + offset,
     xlab = "PC 1 (86.37%)", ylab = "PC 2 (10.75%)")
#par(mfrow = c(1,1))

biplot(nir.prcomp)
biplot(nir.prcomp, xlabs = label, ylabs = lambda)

extremes <- c(14,1,2,3,4)
Xextr <- scale(honeydataframe$NIRS, scale = FALSE)[extremes,]
matplot(lambda, t(Xextr),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "%R (mean-centered)", lty = c(1,4,3,2,2,2,2), col = c(1,2,3,4,4,4,4))





#*******************************************************++++++*****************  PLSR ******************
library(pls)
set.seed(62)
random.rows <- sample(nrow(absorbances))
absorbances.random <- absorbances[random.rows,]
absorbances.random <- absorbances.random[which(absorbances.random$adulterant == "cs"),]
#selected.honey <- absorbances.random[which(absorbances.random$percentage == "0"),]
#absorbances.random <- absorbances.random[which(absorbances.random$percentage != "0"),]
#absorbances.random <- rbind(absorbances.random,selected.honey[which(selected.honey$honey == "4"),])
#absorbances.random <-absorbances.random[-c(grep("*p",absorbances.random$honey)),]
#absorbances.random <- absorbances.random[-11,]
absorbances.random$NIRS <- scale(absorbances.random$NIRS, scale = FALSE)
absorbances.training <- absorbances.random
#absorbances.test <- absorbances.random[1:15,]
concentrations <- as.numeric(levels(absorbances.random$percentage))
adulterant.percentage <- c()
lapply(as.numeric(absorbances.random$percentage),function(value){
  adulterant.percentage <<- append(adulterant.percentage, concentrations[value])})
honey.pls <- plsr(adulterant.percentage ~ NIRS, 
                  data = absorbances.training, ncomp = 4, validation = "LOO")
summary(honey.pls)
plot(RMSEP(honey.pls), legendpos = "topright")
plot(honey.pls, ncomp = 2, asp = 1, line = TRUE)
plot(honey.pls, plottype = "scores", comps = 1:3)
plot(honey.pls, "loadings", comps = 1:2, legendpos = "bottomright", xlab = "nm")
abline(h = 0)
#*************** TEST ************************

absorbances.random$NIRS <- scale(absorbances.random$NIRS, scale = FALSE)
absorbances.training <- absorbances.random[1:80,]
absorbances.test <- absorbances.random[81:109,]
adulterant.percentage <-as.numeric(absorbances.random$percentage)
honey.pls <- plsr(adulterant.percentage[1:80] ~ NIRS, 
                  data = absorbances.training, ncomp = 4, validation = "LOO")
summary(honey.pls)
plot(RMSEP(honey.pls), legendpos = "topright")
plot(honey.pls, ncomp = 2, asp = 1, line = TRUE)
  
  
  library(chemometrics)  
pcr_dcv<-mvr_dcv(vhoney[,3]~.,ncomp = 50,data = honeylog,method = "svdpc")
pcr_plot2 <- plotcompmvr(pcr_dcv)  
rpls = prm_cv(vhoney[,3],honeylog, a = 50, trim = 0.2, plot.opt=TRUE)



