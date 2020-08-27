#BAND-GAP
# The symbol # is a commentary
# Specify your working directory (where the files with the data are) to R .  
#For example "~/Dropbox/Fernando/",equivalent to "/home/Dropbox/Fernando"  
setwd( "")
#A continuation we create a list, temp, with the name of the different .csv files, wich are 
#in the working directory 
temp <- list()
temp <- append(temp, list.files(pattern="*.csv", 
                                  full.names =TRUE))
#The next function read the information in the first row of each file and store it in a string, (sample <-). 
sample <- lapply(temp, read.csv, nrows=1, header=FALSE, sep=","
                   ,stringsAsFactors = FALSE)

# The instruction grep search for matches to argument pattern within each element of a 
#character vector. For example, b_vector <- grep("B",sample) or h_vector <- grep("H",sample)
#The idea is to group the files according to sample nature. The name of the files with 
# hydroxyapatites begin with h, oxides of Ba by B, etc.

#b_vector <-grep("B",sample) #Returns a vector containing a list of all the files whose 
# names have a B
hc_vector <- grep("H","C",sample)
#z_vector <- grep("Z",sample)

#Read and Create a list of data frames with each of the files. Skip = 1, skip the first row 
# where the information about the sample and conditions are
#(two columns: wavelength,%R)
myfiles <- lapply(temp, read.csv, skip = 1, header=TRUE, 
                  colClasses=c("numeric","numeric"),
                  col.names=c("nm","reflectance"))

# plot(1, type = n ...) return an empty plot whitin the coordinates xlim and ylim
plot(1,type = "n", xlim = c(190,700),ylim = c(0,100), xlab = "Wavelength (nm)"
     ,ylab = "%Reflectance")

vcolor <- c("black","red","blue","green") # return a vector, vcolor, with different 
#type of colours. The number of colours will be equal to the number of files in the 
#respective group, for example if length(h_vector) = 4, we need 4 colours
v_position <- c(57.5,60,62.5,65) #return a vector with the legend coordinates

# The next function add lines of different colour for each diffuse reflection spectrum within the groups 
#created with the instruction "grep", for example h_vector[i]  
lapply(1:length(h_vector), function(i) 
  {
  lines(myfiles[[h_vector[i]]]$nm ,myfiles[[h_vector[i]]]$reflectance, 
        type = "l", lty = i) #, col = vcolor[i])
  legend(x=500, y=v_position[i], col = i, lty = i, bty = "n", legend = sample[[h_vector[i]]][2])
}
  )


library(signal) #load signal package

# The next functions calculate the Kubelka-Munk parameters, R, K, S and K-M   
R_list <- lapply(1:length(myfiles), function(i) myfiles[[i]]$reflectance/100)
K_list <- lapply(1:length(myfiles), function(i) (1-R_list[[i]])^2)
S_list <- lapply(1:length(myfiles), function(i) (2*R_list[[i]]))
K_M <- lapply(1:length(myfiles), function(i) (K_list[[i]])/S_list[[i]])


plot(1,type = "n", xlim = c(190,700),ylim = c(0,120), xlab = "Wavelength (nm)"
     ,ylab = "KM_F") #Return an empty plot whitin the coordinates xlim and ylim

v_position2 <- c(20, 25, 30,35) # legend position

#The next function add lines of different colour for each absorption spectrum within the groups 
#created with the instruction "grep", for example h_vector[i]  
lapply(1:length(z_vector), function(i) 
{
  lines(myfiles[[z_vector[i]]]$nm, K_M[[z_vector[i]]], 
        type = "l", lty = i, col = vcolor[i])
  legend(x=450, y=v_position2[i], legend = sample[[z_vector[i]]][2], 
         col = i, lty = i, bty = "n")
}
)


energy <- lapply(1:length(myfiles), function(i) 1240/(myfiles[[i]]$nm)) # Return energy
square.K.M.Energy <- lapply(1:length(myfiles), 
                            function(i) (K_M[[i]]*energy[[i]])^2)

plot

v_position3 <- c(20, 25, 30,35)

lapply(1:length(z_vector), function(i) 
{
  lines(energy[[z_vector[i]]], square.K.M.Energy[[z_vector[i]]], 
        type = "l", lty = i, col = vcolor[i])
  legend(x=450, y=v_position2[i], legend = sample[[z_vector[i]]][2], 
         col = i, lty = i, bty = "n")
}
)
