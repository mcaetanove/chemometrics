 ###########BAND-GAP###########


## Import spectral data --------------------------------------------------------

wd <- "~/Code/BGF"
setwd(wd)
datafile <- "data/UV-Vis.csv"
uvvis.data <- read.csv(datafile, header = TRUE)

# Number of features in the UV-Vis datafile (Files, thickness & intensities)
n.features <- 3
wavel.index <- n.features:ncol(uvvis.data)
n.records <- nrow(uvvis.data)

min.wavel <- as.numeric(sub(".", "", colnames(uvvis.data)[n.features]))
max.wavel <- as.numeric(sub(".", "", tail(colnames(uvvis.data), n = 1)))
lambda <- seq(from = min.wavel, to = max.wavel, by = 1)

## Smooth the data with an improved Savitzky Golay transform -------------------

# Apply Adaptive Degree Polynomial Filter
#
# features      - extra features in the dataframe
# x             - a numeric data.frame, matrix or vector to transform
# diff.order    - differentiation order
# max.order     - maximum polynomial order
# window.size   - window size (must be odd)
#
# Returns original UV-Vis dataframe with smoothed data

adpf <- function(features, x, diff.order, max.order, window.size) {
    library(ADPF)
    func = function(x) ADPF(x, diff.order, max.order, window.size)[,3]
    matrix.smooth <- apply(x, MARGIN = 1, FUN = func)
    uvvis.smooth.data <- cbind(features, absorbance = I(t(matrix.smooth)))
    return(uvvis.smooth.data)
}

diff.order <- 0
max.order <- 4
window.size <- 11
x <- uvvis.data[, wavel.index]
features <- uvvis.data[1:n.features-1]

uvvis.smooth.data <- adpf(features, x, diff.order, max.order, window.size)

## Plot absorbance markers -----------------------------------------------------

# pdf(file = "output/absorbance_plot_markers.pdf", paper = "a4")
cairo_ps("output/absorbance_plot_markers.eps")

par(mar = c(5,6,1,1)+.1)
plot(1,
     type = "n",
     xlim = c(min.wavel, max.wavel),
     ylim = c(0,2),
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = "wavenumber [nm]", 
     ylab = "Absorbance [u.a.]")

col.set <- c("black","red","blue", "violetred")
pch.set <- c(0,1,3,4)
legend <- c();

for (i in 1:n.records) {
    lines(lambda,
          t(uvvis.smooth.data$absorbance[i,]),
          type = "o",
          col = col.set[i],
          lty = i,
          pch = (pch.set[i]),
          cex = 1)

    legend <- c(legend, as.character(uvvis.smooth.data[i,1]))
}

legend(x = 450,
       y = 2.0,
       legend = legend,
       col = col.set,
       lty = 1:n.records,
       pch = c(0:4),
       cex = 1)

invisible(dev.off())
#dev.off()
# browseURL("output/absorbance_plot_markers.eps")

## Calculate absortivity coefficients ------------------------------------------

# Absorption Coefficient Alpha of the film
#
# features      - extra features in the dataframe
# n.records     - number of records
# lambda.length - number of wavelengths
# t             - thickness of the films
# A             - the absorbance
#
# Returns UV-Vis dataframe with absorption coefficients

acaf <- function(features, n.records, lambda.length, t, A) {
    m.alpha <- matrix(, nrow = 0, ncol = length(lambda))
    for (i in 1:n.records) {
        m <- mapply(function(t, A) (2.303 / t) * A, t, A)
        m.alpha <- rbind(m.alpha, m)
    }
    alpha.data <- cbind(features, alpha = I(m.alpha))
    return(alpha.data)
}

lambda.length <- length(lambda)
t <- uvvis.smooth.data$thickness[i]
A <- uvvis.smooth.data$absorbance[i,]

alpha.data <- acaf(features, n.records, lambda.length, t, A)

# Plot absorption coefficient vs energy in eV ----------------------------------

# pdf(file = "output/absorbance_coefficient_plot_markers.pdf", paper = "a4")
cairo_ps("output/absorbance_coefficients_plot_markers.eps")

par(mar=c(5,6,2,1) + .1)

# Vector of energy in eV
energy <- 1240 / lambda
col.set <- c("black","red","green","blue", "violetred")
pch.set <- c(0,1,2,3,4)

matplot(energy,
        t(I(alpha.data$alpha))/10E+3,
        type = "o",
        lty = c(1:n.records),
        pch = pch.set,
        col = col.set,
        cex = 0.75,
        cex.axis = 1.5,
        cex.lab = 1.5,
        xlab = expression(paste("Photon energy(h",nu,")  [eV]")),
        ylab = expression(paste(alpha %*% 10^3, "  [m"^-1,"]")))

legend(x = 2,
       y = 15,
       legend = legend,
       col = col.set,
       lty= 1:n.records,
       pch = pch.set,
       cex=1)

abline(-158.1846, 28.99253, col="red")
abline(-158.1981, 29.2276, col="red")
abline(-104.0234, 19.68363, col="red")
abline(-104.9601, 19.76176, col="red")
abline(-89.01537, 17.1056, col="red")
abline(h=0)

#dev.off()
invisible(dev.off())

# TODO
finite.differences <- function(x, y) {
    n <- length(x)
    fdx <- vector(length = n)
    for (i in 2:n) {
        fdx[i - 1] <- (y[i - 1] - y[i]) / (x[i - 1] - x[i])
    }
    fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
    return(fdx)
}

for (i in 1:n.records) {
cairo_ps(paste("output/absorbance_coefficient_plot_markers", i, ".eps", sep =""))

    par(mar = c(5,6,2,1) + .1)

    plot(1,
        type = "n",
        xlim = c(1,6.5),
        ylim = c(0,15),
        cex.axis = 1.5,
        cex.lab = 1.5,
        xlab = expression(paste("Photon energy [h",nu,"]")),
        ylab = expression(paste(alpha %*% 10^3, " [m" ^-1,"]")))

    lines(energy,
        t(alpha.data$alpha[i,]) / 10E+3,
        type = "o",
        col = col.set[i],
        pch = pch.set[i],
        lty = i)

    ac <- data.frame(x = energy, y = as.vector(alpha.data$alpha[i,] / 10E+3))
    slopes <- diff(ac[,2]) / diff(ac[,1])
    slopes.k <- kmeans(slopes,3)
    print(slopes.k$centers)
    slopes.max <- max(abs(slopes.k$centers))
    print(slopes.max)
    slopes.v <- which(slopes.k$cluster == match(slopes.max, slopes.k$centers))

    lp <- ac[slopes.v,]
    points(lp, pch=15)
    abline(coef(lm(y~x,lp)), col="red")

    #min <- 5.6
    #max <- 5.85
    #print(min)
    #print(max)
    #abline(v = c(1240 /min, 1240 /max))

    #min
    #max
    #lim.inf <- which(abs(energy - min) == min(abs(energy - min))) #find the index in the energy array corresponding
    #lim.sup <- which(abs(energy - max) == min(abs(energy - max)))

    #model <- lm((alpha.data$alpha[i,lim.inf:lim.sup]) / 10E+3~energy[lim.inf:lim.sup])
    #print(summary(model))
    #edge <- - model$coefficients[1] / model$coefficients[2]# Computing the band gap
    abline(h = 0) # trace a horizontal line at y = 0
   # abline(model$coefficients[1], model$coefficients[2],col = "red") # trace the fitted model, red line
    #print(edge)
    #print(model$coefficients[1])
    #print(model$coefficients[2])

    dev.off()
}
