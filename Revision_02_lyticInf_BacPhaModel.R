# title: "lytic infection"
# ODE model

rm(list = ls())

##### 'libraries' #####
library(odeintr)
library(EMD)       #number of extrema; local minima/maxima
library(ggplot2)

##### 'start values' #####
Resource <- 2
Bacteria <- 0.83
Bacteria2 <- 0.83
Infected <- 0.17
Infected2 <- 0.17
Phage <- 0.00031
Phage2 <- 0.00031
xStart <- c(Resource, Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)

##### 'parameter' #####

##### 'densities' 
# all densities are normalized by half. sat density 5.2333e06 µm^3 mL^-03'
S <- Resource  # supply concentration
K <- 1         # halfsaturation density

##### 'rates' 
# normalized by metabolic rate 0.002365714285714 hr^-1
D <- 1/24    # turnover rate, mean nutrient input
m <- 0.0206        # metabolic rate
i <- 5.23        # infection/ adsorption rate of fast-growing bacteria
i1 <- 5.23      # infection/ adsorption rate of slow-growing bacteria
s <- 0.330     # bacterial lysis
d <- 0.0866     # phage decay

##### 'conversion rates'
y <- 7.5          # metabolic scaling constant for fast-growing bacteria
y1 <- 4         # metabolic scaling constant for slow-growing bacteria
c <- 3.14e-04     # ratio of bacteria-phage biovolume
b <- 0.02         # phage burst biovolume

a <- 0     # amplitude, variation of nutrient fluctuation

TL <- 168     # turnover period
p <- 2*pi/TL  # frequency
# 24
# 168
# 720
# 8760

##### 'Model of phage-bacteria-interaction' #####
BasicPhage.sys <- '
double XX0 = x[0];
double XX1 = x[1];
double XX2 = x[2];
double XX3 = x[3];
double XX4 = x[4];
double XX5 = x[5];
double XX6 = x[6];

if (XX1 <= pow(10,-7)){
double XX1 = 0.0;
}
if (XX2 <= pow(10,-7)){
double XX2 = 0.0;
}
if (XX3 <= pow(10,-7)){
double XX3 = 0.0;
}
if (XX4 <= pow(10,-7)){
double XX4 = 0.0;
}

dxdt[0] = (a*D*sin(p*t)+D)*(S-XX0) - y*m*XX0/(K+XX0)*XX1 - y1*m*XX0/(K+XX0)*XX2 + (1-b)*s*XX3 + (1-b)*s*XX4;

if(XX1 <= pow(10,-7)){
  dxdt[1] = {0.0};
}
else{
  dxdt[1] = y*m*XX0/(K+XX0)*XX1 - m*XX1 - i*XX1*XX5;
}

if(XX2 <= pow(10,-7)){
  dxdt[2] = {0.0};
}
else{
  dxdt[2] = y1*m*XX0/(K+XX0)*XX2 - m*XX2 - i1*XX2*XX6;
}

if(XX3 <= pow(10,-7)){
  dxdt[3] = {0.0};
}
else{
  dxdt[3] = i*(1+c)*XX1*XX5 - m*XX3 - s*XX3;
}

if(XX4 <= pow(10,-7)){
  dxdt[4] = {0.0};
}
else{
  dxdt[4] = i1*(1+c)*XX2*XX6 - m*XX4 - s*XX4;
}

dxdt[5] = b*s*XX3 - c*i*XX1*XX5 - d*XX5;

dxdt[6] = b*s*XX4 - c*i1*XX2*XX6 - d*XX6;

'
compile_sys(name = "BasicPhage",
            sys = BasicPhage.sys,
            pars = c("a", "p", "D", "S", "y", "y1", "K", "m",
                     "i", "i1", "c", "b", "s", "d"),
            method = "rk54_a")

BasicPhage_set_params(a = 0, p = p, D = D, S = S, y = y, y1 = y1, m = m, K = K,
                      i = i, i1 = i1, c = c, b = b,  s = s, d = d)

TS.01 <- BasicPhage_adap(xStart, 2000)

TS.01[, 3:8][TS.01[, 3:8] < 10e-7] = 0.0
##### 'Plot Data' #####
par(mfrow = c(1,2))

plot(c(TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time), 
     log10(c(TS.01$X1, TS.01$X2, TS.01$X3, TS.01$X4, TS.01$X5, TS.01$X6, TS.01$X7)),
     type = "n",
     xlab = "Time",
     ylab = "log10 norm. biovolume",
     cex.lab=1.4)
lines(TS.01$Time, log10(TS.01$X1), col = "navy", lwd = 2)
lines(TS.01$Time, log10(TS.01$X2), col = "green2", lwd = 2)
lines(TS.01$Time, log10(TS.01$X3), col = "goldenrod1", lwd = 2)
lines(TS.01$Time, log10(TS.01$X4), col = "green4", lwd = 2)
lines(TS.01$Time, log10(TS.01$X5), col = "darkgoldenrod3", lwd = 2)
lines(TS.01$Time, log10(TS.01$X6), col = "brown2", lwd = 2)
lines(TS.01$Time, log10(TS.01$X7), col = "firebrick", lwd = 2)

plot(c(TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time, TS.01$Time), 
     c(TS.01$X1, TS.01$X2, TS.01$X3, TS.01$X4, TS.01$X5, TS.01$X6, TS.01$X7),
     type = "n",
     xlab = "Time",
     ylab = "norm. biovolume")
lines(TS.01$Time, TS.01$X1, col = "navy", lwd = 2)
lines(TS.01$Time, TS.01$X2, col = "green2", lwd = 2)
lines(TS.01$Time, TS.01$X3, col = "goldenrod1", lwd = 2)
lines(TS.01$Time, TS.01$X4, col = "green4", lwd = 2)
lines(TS.01$Time, TS.01$X5, col = "darkgoldenrod3", lwd = 2)
lines(TS.01$Time, TS.01$X6, col = "brown2", lwd = 2)
lines(TS.01$Time, TS.01$X7, col = "firebrick", lwd = 2)

legend("bottomright", 
       legend=paste(c("Nutrients", "fast bacteria", "slow bacteria", "lysogenic bacteria (fast)", "lysogenic bacteria (slow)", "Phage (fast)", "Phage (slow)")), 
       col= c("navy", "green2",  "goldenrod1", "green4", "darkgoldenrod3", "brown2", "firebrick"), 
       lty=1, ncol=1, cex=0.5, xpd = T)

##### 'Nutrient Bifurcation' #####
TSlength <- 200000
TSanalyzed <- 0.999
log10.NutCon <- seq(log10(0.6), log10(30), length = 1000)
NutCon <- 10^log10.NutCon

Nut.min <- c()
Nut.max <- c()
Bac.min <- c()
Bac.max <- c()
iBac.min <- c()
iBac.max <- c()
Pha.min <- c()
Pha.max <- c()
Bac2.min <- c()
Bac2.max <- c()
iBac2.min <- c()
iBac2.max <- c()
Pha2.min <- c()
Pha2.max <- c()

for(j in 1:length(NutCon)){
  xStart <- c(NutCon[j], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)
  BasicPhage_set_params(a = 0, p = p, D = D, S = NutCon[j], y = y, y1 = y1, m = m, K = K,
                        i = i, i1 = i1, c = c, b = b,  s = s, d = d)
  
  Sim.01 <- BasicPhage_adap(xStart, TSlength)
  Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
  Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
  
  if(extrema(Sim.01.Sub$X1)$nextreme == 0){
    Nut.min <- c(Nut.min, Sim.01.Sub$X1[1:10]) 
    Nut.max <- c(Nut.max, Sim.01.Sub$X1[1:10])}else{
      Nut.min <- c(Nut.min, Sim.01.Sub$X1[extrema(Sim.01.Sub$X1)$minindex[1:10]])
      Nut.max <- c(Nut.max, Sim.01.Sub$X1[extrema(Sim.01.Sub$X1)$maxindex[1:10]])
    }
  
  if(extrema(Sim.01.Sub$X2)$nextreme == 0){
    Bac.min <- c(Bac.min, Sim.01.Sub$X2[1:10]) 
    Bac.max <- c(Bac.max, Sim.01.Sub$X2[1:10])}else{
      Bac.min <- c(Bac.min, Sim.01.Sub$X2[extrema(Sim.01.Sub$X2)$minindex[1:10]])
      Bac.max <- c(Bac.max, Sim.01.Sub$X2[extrema(Sim.01.Sub$X2)$maxindex[1:10]])
    }
  
  if(extrema(Sim.01.Sub$X3)$nextreme == 0){
    Bac2.min <- c(Bac2.min, Sim.01.Sub$X3[1:10]) 
    Bac2.max <- c(Bac2.max, Sim.01.Sub$X3[1:10])}else{
      Bac2.min <- c(Bac2.min, Sim.01.Sub$X3[extrema(Sim.01.Sub$X3)$minindex[1:10]])
      Bac2.max <- c(Bac2.max, Sim.01.Sub$X3[extrema(Sim.01.Sub$X3)$maxindex[1:10]])
    }
  
  if(extrema(Sim.01.Sub$X4)$nextreme == 0){
    iBac.min <- c(iBac.min, Sim.01.Sub$X4[1:10]) 
    iBac.max <- c(iBac.max, Sim.01.Sub$X4[1:10])}else{
      iBac.min <- c(iBac.min, Sim.01.Sub$X4[extrema(Sim.01.Sub$X4)$minindex[1:10]])
      iBac.max <- c(iBac.max, Sim.01.Sub$X4[extrema(Sim.01.Sub$X4)$maxindex[1:10]])
    }
  
  if(extrema(Sim.01.Sub$X5)$nextreme == 0){
    iBac2.min <- c(iBac2.min, Sim.01.Sub$X5[1:10]) 
    iBac2.max <- c(iBac2.max, Sim.01.Sub$X5[1:10])}else{
      iBac2.min <- c(iBac2.min, Sim.01.Sub$X5[extrema(Sim.01.Sub$X5)$minindex[1:10]])
      iBac2.max <- c(iBac2.max, Sim.01.Sub$X5[extrema(Sim.01.Sub$X5)$maxindex[1:10]])
    }
  
  if(extrema(Sim.01.Sub$X6)$nextreme == 0){
    Pha.min <- c(Pha.min, Sim.01.Sub$X6[1:10]) 
    Pha.max <- c(Pha.max, Sim.01.Sub$X6[1:10])}else{
      Pha.min <- c(Pha.min, Sim.01.Sub$X6[extrema(Sim.01.Sub$X6)$minindex[1:10]])
      Pha.max <- c(Pha.max, Sim.01.Sub$X6[extrema(Sim.01.Sub$X6)$maxindex[1:10]])
    }
  
  if(extrema(Sim.01.Sub$X7)$nextreme == 0){
    Pha2.min <- c(Pha2.min, Sim.01.Sub$X7[1:10]) 
    Pha2.max <- c(Pha2.max, Sim.01.Sub$X7[1:10])}else{
      Pha2.min <- c(Pha2.min, Sim.01.Sub$X7[extrema(Sim.01.Sub$X7)$minindex[1:10]])
      Pha2.max <- c(Pha2.max, Sim.01.Sub$X7[extrema(Sim.01.Sub$X7)$maxindex[1:10]])
    }
  
  #print(j)
}

NutCon2 <- sort(rep(NutCon, 10))

par(mfrow = c(1,1))
plot(log10(rep(NutCon2, 14)), 
     log10(c(Nut.min, Nut.max, Bac.min, Bac.max, Bac2.min, Bac2.max, iBac.min, iBac.max, iBac2.min, iBac2.max, Pha.min, Pha.max, Pha2.min, Pha2.max)),
     type = "n",
     ylim = c(-10, 3),
     xlab = "log10 resource supply (N))",
     ylab = "log10 normalized Biovolume")

points(log10(c(NutCon2, NutCon2)), log10(c(Nut.min, Nut.max)), col = "navy", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Bac.min, Bac.max)), col = "green2", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Bac2.min, Bac2.max)), col = "goldenrod1", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(iBac.min, iBac.max)), col = "green4", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(iBac2.min, iBac2.max)), col = "darkgoldenrod3", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Pha.min, Pha.max)), col = "brown2", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Pha2.min, Pha2.max)), col = "firebrick", cex = 0.25, pch = 20)
abline(v = c(-0.25, 0, 0.25, 0.50, 0.75))

legend("bottomright", 
       legend=paste(c("Nutrients", "fast bacteria", "slow bacteria", "lysogenic bacteria (fast)", "lysogenic bacteria (slow)", "Phage (fast)", "Phage (slow)")), 
       col= c("navy", "green2",  "goldenrod1", "green4", "darkgoldenrod3", "brown2", "firebrick"), 
       lty=1, ncol=1, cex=0.5, xpd = T)

##### 'Nutrients and Amplitude' #####
TSlength <- 200000
TSanalyzed <- 0.999
log10.NutCon <- seq(log10(0.60), log10(8), length = 400)
NutCon <- 10^log10.NutCon
Amplitude <- seq(0, 0.99, length = 400)

TL <- 720
p <- 2*pi/TL
      
results1 <- data.frame(Nutrients = c(),
                      Amplitude = c(),
                      SpecNum = c(),
                      Bac = c(),
                      Bac2 = c(),
                      iBac = c(),
                      iBac2 = c(),
                      Pha = c(),
                      Pha2 = c())
      
for (k in 1:length(NutCon)) {
  for(j in 1:length(Amplitude)){
    
    xStart <- c(NutCon[k], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)
    
    BasicPhage_set_params(a = Amplitude[j], p = p, D = D, S = NutCon[k], y = y, y1 = y1, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b,  s = s, d = d)
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Surviver = c()
    Speccount <- 0
    
    for(l in 3:8){
      if(Sim.01.Sub[[l]] > 0.0){
        Surviver <- c(Surviver, 1)
      }else{
        Surviver <- c(Surviver, 0)  
      }
    }
    
    Speccount <- sum(Surviver)
    
    results.loop <- data.frame(Nutrients = NutCon[k],
                               Amplitude = Amplitude[j],
                               SpecNum = Speccount,
                               Bac = Surviver[1],
                               Bac2 = Surviver[2],
                               iBac = Surviver[3],
                               iBac2 = Surviver[4],
                               Pha = Surviver[5],
                               Pha2 = Surviver[6])
    
    results1 <- rbind(results1, results.loop)
    
    #print(paste("Nutrients = ",k,"; Amplitude = ",j))
  }
}
str(results1)

ggplot(results1, aes(Nutrients, Amplitude, fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = rainbow(7),
                       breaks = seq(0, 6, 1),
                       limits = c(0, 6)) +
  scale_x_continuous(trans = 'log10')

##### 'Nutrients and Frequency' #####
TSlength <- 200000
TSanalyzed <- 0.999
log10.NutCon <- seq(log10(0.60), log10(8), length = 400)
NutCon <- 10^log10.NutCon
log10.Frequency <- seq(log10(2*pi/(24)), log10(2*pi/(8760)), length = 400)
Fre <- 10^log10.Frequency

# Different metabolic scaling constants
y <- 7.5  # 3.75 # 15        #fast-growing bacteria
y1 <- 4   # 2 # 8            #slow-growing bacteria
  
results2 <- data.frame(Nutrients = c(),
                      Frequency = c(),
                      SpecNum = c(),
                      Bac = c(),
                      Bac2 = c(),
                      iBac = c(),
                      iBac2 = c(),
                      Pha = c(),
                      Pha2 = c())
  
for (k in 1:length(NutCon)) {
  for(j in 1:length(Fre)){
      
    xStart <- c(NutCon[k], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)
      
    BasicPhage_set_params(a = 0.9, p = Fre[j], D = D, S = NutCon[k], y = y, y1 = y1, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b,  s = s, d = d)
      
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Surviver = c()
    Speccount <- 0
    
    for(l in 3:8){
      if(Sim.01.Sub[[l]] > 0.0){
        Surviver <- c(Surviver, 1)
      }else{
        Surviver <- c(Surviver, 0)  
      }
    }
    
    Speccount <- sum(Surviver)
      
    results.loop <- data.frame(Nutrients = NutCon[k],
                               Frequency = Fre[j],
                               SpecNum = Speccount,
                               Bac = Surviver[1],
                               Bac2 = Surviver[2],
                               iBac = Surviver[3],
                               iBac2 = Surviver[4],
                               Pha = Surviver[5],
                               Pha2 = Surviver[6])
      
    results2 <- rbind(results2, results.loop)
      
    #print(paste("Nutrients = ",k,"; Frequency = ",j))
   }
}
str(results2)
 
ggplot(results2, aes(Nutrients, (2*pi/Frequency/0.0206), fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = rainbow(7),
                       breaks = seq(0, 6, 1),
                       limits = c(0, 6)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10')