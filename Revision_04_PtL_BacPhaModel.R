# title: "phage-bacteria interaction"
# ODE model

rm(list = ls())

##### 'libraries' #####
library(odeintr)
library(EMD)       #number of extrema; local minima/maxima
library(ggplot2)
library(viridis)

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
D <- 1/24       # turnover rate, mean nutrient input
m <- 0.0206     # metabolic rate
i <- 5.23       # infection/ adsorption rate of fast-growing bacteria
i1 <- 5.23      # infection/ adsorption rate of slow-growing bacteria
s <- 0.33       # bacterial lysis
SH <- 4
d <- 0.087      # phage decay

##### 'conversion rates'
y <- 7.5          # metabolic scaling constant for fast-growing bacteria
y1 <- 4           # metabolic scaling constant for slow-growing bacteria
y2 <- 7.5         # metabolic scaling constant for lysogenic fast-growing bacteria
y3 <- 4           # metabolic scaling constant for lysogenic fast-growing bacteria
c <- 3.14e-04     # ratio of bacteria-phage biovolume
b <- 0.02         # phage burst biovolume

a <- 0        # amplitude, variation of nutrient fluctuation
TL <- 168     # period
p <- 2*pi/TL  # oscillation frequency
# 24 hr
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

dxdt[0] = (a*D*sin(p*t)+D)*(S-XX0) - y*m*XX0/(K+XX0)*XX1 - y1*m*XX0/(K+XX0)*XX2 - y2*m*XX0/(K+XX0)*XX3 - y3*m*XX0/(K+XX0)*XX4 + (1-b)*s*pow((XX1+XX3),2)/(pow((XX1+XX3),2)+SH)*XX3 + (1-b)*s*pow((XX2+XX4),2)/(pow((XX2+XX4),2)+SH)*XX4;

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
else if(XX1+XX3 >= WP){
  dxdt[3] = i*(1+c)*XX1*XX5 - m*XX3 - s*pow((XX1+XX3),2)/(pow((XX1+XX3),2)+SH)*XX3;
}
else{
  dxdt[3] = i*(1+c)*XX1*XX5 + y2*m*XX0/(K+XX0)*XX3 - m*XX3 - s*pow((XX1+XX3),2)/(pow((XX1+XX3),2)+SH)*XX3;
}

if(XX4 <= pow(10,-7)){
  dxdt[4] = {0.0};
}
else if(XX2+XX4 >= WP){
  dxdt[4] = i1*(1+c)*XX2*XX6 - m*XX4 - s*pow((XX2+XX4),2)/(pow((XX2+XX4),2)+SH)*XX4;
}
else{
  dxdt[4] = i1*(1+c)*XX2*XX6 + y3*m*XX0/(K+XX0)*XX4 - m*XX4 - s*pow((XX2+XX4),2)/(pow((XX2+XX4),2)+SH)*XX4;
}

dxdt[5] = b*s*pow((XX1+XX3),2)/(pow((XX1+XX3),2)+SH)*XX3 - c*i*XX1*XX5 - d*XX5;

dxdt[6] = b*s*pow((XX2+XX4),2)/(pow((XX2+XX4),2)+SH)*XX4 - c*i1*XX2*XX6 - d*XX6;
' 
compile_sys(name = "BasicPhage",
            sys = BasicPhage.sys,
            pars = c("a", "p", "D", "S", "y", "y1", "y2", "y3", "K", "m",
                     "i", "i1", "c", "b", "s", "d", "SH", "WP"),
            method = "rk54_a")

BasicPhage_set_params(a = 0, p = p, D = D, S = S, y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                      i = i, i1 = i1, c = c, b = b,  s = s, d = d, SH = SH, WP = 0.2010076)

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
abline(v=250)

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

##### 'Nutrientsupply and Switch' #####
TSlength <- 10000
TSanalyzed <- 0.99

TL <- 720
p <- 2*pi/TL

log10.SW <- seq(log10(0.05), log10(40), length = 100)
SwitchPoint <- 10^log10.SW
log10.NutCon <- seq(log10(0.1), log10(10), length = 100)
NutCon <- 10^log10.NutCon

results <- data.frame(NutrientSupply = c(),
                      SwitchL = c(),
                      concNut = c(),
                      concBac = c(),
                      concBac2 = c(),
                      conciBac = c(),
                      conciBac2 = c(),
                      concPha = c(),
                      concPha2 = c(),
                      SpecNum = c()
)

for(j in 1:length(NutCon)){
  for (k in 1:length(SwitchPoint)){
  xStart <- c(NutCon[j], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)  
  BasicPhage_set_params(a = 0, p = p, D = D, S = NutCon[j], y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                        i = i, i1 = i1, c = c, b = b, s = s, d = d, SH = SH, WP = SwitchPoint[k])
  
  Sim.01 <- BasicPhage_adap(xStart, TSlength)
  Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
  Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
  
  Speccount <- 0
  
  for(l in 3:8){
    if(tail(Sim.01.Sub[,l],1) > 0.0){
      Speccount <- Speccount+1
    }
  }
  
  results.loop <- data.frame(NutrientSupply = NutCon[j],
                             SwitchL = SwitchPoint[k],
                             concNut = mean(Sim.01.Sub$X1),
                             concBac = mean(Sim.01.Sub$X2),
                             concBac2 = mean(Sim.01.Sub$X3),
                             conciBac = mean(Sim.01.Sub$X4),
                             conciBac2 = mean(Sim.01.Sub$X5),
                             concPha = mean(Sim.01.Sub$X6),
                             concPha2 = mean(Sim.01.Sub$X7),
                             SpecNum = Speccount
  )
  results <- rbind(results, results.loop)
  print(paste("Nut = ",j,"; Switch = ",k))
  }
}
SH <- 4

ggplot(results, aes(NutrientSupply, SwitchL, fill = SpecNum))+
  geom_raster(interpolate = T) +
  geom_hline(yintercept = c((sqrt(4/(0.33/0.165-1))), (sqrt(4/(0.33/0.0033-1))), (sqrt(4/(0.33/0.001-1))))) +
  scale_fill_gradientn(name = "", colors = magma(7)) +
  scale_x_continuous(trans = "log10", expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", expand = c(0, 0))

##### 'Lysis Rate Bifurcation' #####
TSlength <- 10000
TSanalyzed <- 0.99
log10.Lysis <- seq(log10(0.1), log10(1), length = 1000)
Lysis <- 10^log10.Lysis

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

for(j in 1:length(Lysis)){
  BasicPhage_set_params(a = 0, p = p, D = D, S = S, y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                        i = i, i1 = i1, c = c, b = b, s = Lysis[j], d = d, SH = SH, WP = 0.2010076)
  
  Sim.01 <- BasicPhage_adap(xStart, TSlength)
  Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
  Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
  
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

Lysis2 <- sort(rep(Lysis, 10))

par(mfrow = c(1,1))
plot(log10(rep(Lysis2, 12)), 
     log10(c(Bac.min, Bac.max, Bac2.min, Bac2.max, iBac.min, iBac.max, iBac2.min, iBac2.max, Pha.min, Pha.max, Pha2.min, Pha2.max)),
     type = "n",
     ylim = c(-10, 2),
     xlab = "log10(Lysis (s))",
     ylab = "normalized Biovolume")

points(log10(c(Lysis2, Lysis2)), log10(c(Bac.min, Bac.max)), col = "green2", cex = 0.25, pch = 20)
points(log10(c(Lysis2, Lysis2)), log10(c(Bac2.min, Bac2.max)), col = "goldenrod1", cex = 0.25, pch = 20)
points(log10(c(Lysis2, Lysis2)), log10(c(iBac.min, iBac.max)), col = "green4", cex = 0.25, pch = 20)
points(log10(c(Lysis2, Lysis2)), log10(c(iBac2.min, iBac2.max)), col = "darkgoldenrod3", cex = 0.25, pch = 20)
points(log10(c(Lysis2, Lysis2)), log10(c(Pha.min, Pha.max)), col = "brown2", cex = 0.25, pch = 20)
points(log10(c(Lysis2, Lysis2)), log10(c(Pha2.min, Pha2.max)), col = "firebrick", cex = 0.25, pch = 20)

legend("bottomright", 
       legend=paste(c("fast bacteria", "slow bacteria", "lysogenic bacteria (fast)", "lysogenic bacteria (slow)", "Phage (fast)", "Phage (slow)")), 
       col= c("green2",  "goldenrod1", "green4", "darkgoldenrod3", "brown2", "firebrick"), 
       lty=1, ncol=1, cex=0.5, xpd = T)

##### 'Resource supply and Lysis' #####
TSlength <- 10000
TSanalyzed <- 0.99
TL <- 720
p <- 2*pi/TL 
log10.Lysis <- seq(log10(0.1), log10(1), length = 100)
Lysis <- 10^log10.Lysis
log10.NutCon <- seq(log10(0.6), log10(8), length = 100)
NutCon <- 10^log10.NutCon

results <- data.frame(NutrientSupply = c(),
                      LysisR = c(),
                      concNut = c(),
                      concBac = c(),
                      concBac2 = c(),
                      conciBac = c(),
                      conciBac2 = c(),
                      concPha = c(),
                      concPha2 = c(),
                      SpecNum = c()
)

for(j in 1:length(NutCon)){
  for (k in 1:length(Lysis)){
    xStart <- c(NutCon[j], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)  
    BasicPhage_set_params(a = 0, p = p, D = D, S = NutCon[j], y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b, s = Lysis[k], d = d, SH = SH, WP = 19.89975)
    
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Speccount <- 0
    
    for(l in 3:8){
      if(tail(Sim.01.Sub[,l],1) > 0.0){
        Speccount <- Speccount+1
      }
    }
    
    results.loop <- data.frame(NutrientSupply = NutCon[j],
                               LysisR = Lysis[k]*0.0206,
                               concNut = mean(Sim.01.Sub$X1),
                               concBac = mean(Sim.01.Sub$X2),
                               concBac2 = mean(Sim.01.Sub$X3),
                               conciBac = mean(Sim.01.Sub$X4),
                               conciBac2 = mean(Sim.01.Sub$X5),
                               concPha = mean(Sim.01.Sub$X6),
                               concPha2 = mean(Sim.01.Sub$X7),
                               SpecNum = Speccount
    )
    results <- rbind(results, results.loop)
    #print(paste("Nut = ",j,"; Switch = ",k))
  }
}

ggplot(results, aes(NutrientSupply, LysisR, fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = magma(7)) +
  scale_x_continuous(trans = "log10", expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", expand = c(0, 0))

##### 'IndHalfSat Bifurcation' #####
TSlength <- 10000
TSanalyzed <- 0.99
IndHalfSat <- seq(1, 5, length = 1000)

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

for(j in 1:length(IndHalfSat)){
  BasicPhage_set_params(a = 0, p = p, D = D, S = S, y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                        i = i, i1 = i1, c = c, b = b, s = s, d = d, SH = IndHalfSat[j], WP = 0.2010076)
  
  Sim.01 <- BasicPhage_adap(xStart, TSlength)
  Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
  Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
  
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
IndHalfSat2 <- sort(rep(IndHalfSat, 10))

par(mfrow = c(1,1))
plot(rep(IndHalfSat2, 12), 
     log10(c(Bac.min, Bac.max, Bac2.min, Bac2.max, iBac.min, iBac.max, iBac2.min, iBac2.max, Pha.min, Pha.max, Pha2.min, Pha2.max)),
     type = "n",
     ylim = c(-10, 2),
     xlab = "log10(IndHalfSat (s))",
     ylab = "normalized Biovolume")

points(c(IndHalfSat2, IndHalfSat2), log10(c(Bac.min, Bac.max)), col = "green2", cex = 0.25, pch = 20)
points(c(IndHalfSat2, IndHalfSat2), log10(c(Bac2.min, Bac2.max)), col = "goldenrod1", cex = 0.25, pch = 20)
points(c(IndHalfSat2, IndHalfSat2), log10(c(iBac.min, iBac.max)), col = "green4", cex = 0.25, pch = 20)
points(c(IndHalfSat2, IndHalfSat2), log10(c(iBac2.min, iBac2.max)), col = "darkgoldenrod3", cex = 0.25, pch = 20)
points(c(IndHalfSat2, IndHalfSat2), log10(c(Pha.min, Pha.max)), col = "brown2", cex = 0.25, pch = 20)
points(c(IndHalfSat2, IndHalfSat2), log10(c(Pha2.min, Pha2.max)), col = "firebrick", cex = 0.25, pch = 20)

legend("bottomright", 
       legend=paste(c("fast bacteria", "slow bacteria", "lysogenic bacteria (fast)", "lysogenic bacteria (slow)", "Phage (fast)", "Phage (slow)")), 
       col= c("green2",  "goldenrod1", "green4", "darkgoldenrod3", "brown2", "firebrick"), 
       lty=1, ncol=1, cex=0.5, xpd = T)

##### 'Resource supply and IndHalfSat' #####
TSlength <- 10000
TSanalyzed <- 0.99

TL <- 720
p <- 2*pi/TL

IndHalfSat <- seq(1, 5, length = 100)
log10.NutCon <- seq(log10(0.6), log10(8), length = 100)
NutCon <- 10^log10.NutCon

results <- data.frame(NutrientSupply = c(),
                      HalfSaturation = c(),
                      concNut = c(),
                      concBac = c(),
                      concBac2 = c(),
                      conciBac = c(),
                      conciBac2 = c(),
                      concPha = c(),
                      concPha2 = c(),
                      SpecNum = c()
)

for(j in 1:length(NutCon)){
  for (k in 1:length(IndHalfSat)){
    xStart <- c(NutCon[j], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)  
    BasicPhage_set_params(a = 0, p = p, D = D, S = NutCon[j], y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b, s = s, d = d, SH = IndHalfSat[k], WP = 0.2010076)
    
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Speccount <- 0
    
    for(l in 3:8){
      if(tail(Sim.01.Sub[,l],1) > 0.0){
        Speccount <- Speccount+1
      }
    }
    
    results.loop <- data.frame(NutrientSupply = NutCon[j],
                               HalfSaturation = IndHalfSat[k],
                               concNut = mean(Sim.01.Sub$X1),
                               concBac = mean(Sim.01.Sub$X2),
                               concBac2 = mean(Sim.01.Sub$X3),
                               conciBac = mean(Sim.01.Sub$X4),
                               conciBac2 = mean(Sim.01.Sub$X5),
                               concPha = mean(Sim.01.Sub$X6),
                               concPha2 = mean(Sim.01.Sub$X7),
                               SpecNum = Speccount
    )
    results <- rbind(results, results.loop)
    #print(paste("Nut = ",j,"; Switch = ",k))
  }
}
tail(results)

ggplot(results, aes(NutrientSupply, HalfSaturation, fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = magma(7)) +
  scale_x_continuous(trans = "log10", expand = c(0, 0))

##### 'IndHalfSat & SW' #####
TSlength <- 10000
TSanalyzed <- 0.99
log10.SW <- seq(log10(0.05), log10(40), length = 50)
SW <- 10^log10.SW
log10.IndHalfSat <- seq(log10(0.5), log10(5), length = 50)
IndHalfSat <- 10^log10.IndHalfSat

results <- data.frame(SH = c(),
                      Switch = c(),
                      concNut = c(),
                      concBac = c(),
                      concBac2 = c(),
                      conciBac = c(),
                      conciBac2 = c(),
                      concPha = c(),
                      concPha2 = c(),
                      SpecNum = c()
)

for(j in 1:length(IndHalfSat)){
  for (k in 1:length(SW)){
    BasicPhage_set_params(a = 0, p = p, D = D, S = S, y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b, s = s, d = d, SH = IndHalfSat[j], WP = SW[k])
    
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Speccount <- 0
    
    for(l in 3:8){
      if(tail(Sim.01.Sub[,l],1) > 0.0){
        Speccount <- Speccount+1
      }
    }
    
    results.loop <- data.frame(SH = IndHalfSat[j],
                               Switch = SW[k],
                               concNut = mean(Sim.01.Sub$X1),
                               concBac = mean(Sim.01.Sub$X2),
                               concBac2 = mean(Sim.01.Sub$X3),
                               conciBac = mean(Sim.01.Sub$X4),
                               conciBac2 = mean(Sim.01.Sub$X5),
                               concPha = mean(Sim.01.Sub$X6),
                               concPha2 = mean(Sim.01.Sub$X7),
                               SpecNum = Speccount
    )
    results <- rbind(results, results.loop)
    #print(paste("Nut = ",j,"; Switch = ",k))
  }
}
tail(results)

ggplot(results, aes(SH, Switch, fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = magma(7), breaks=seq(0, 6, 1),
                       limits=c(0, 6)) +
  scale_x_continuous(trans = "log10", expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", expand = c(0, 0))

##### 'Nutrient Bifurcation' #####
TSlength <- 10000
TSanalyzed <- 0.99
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
  
  BasicPhage_set_params(a = 0, p = p, D = D, S = NutCon[j], y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                        i = i, i1 = i1, c = c, b = b,  s = s, d = d, SH = SH, WP = 0.1102636)
  
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
  
  print(j)
}
tail(Pha.min)
NutCon2 <- sort(rep(NutCon, 10))
par(mfrow = c(1,1))
plot(log10(rep(NutCon2, 14)), 
     log10(c(Nut.min, Nut.max, Bac.min, Bac.max, Bac2.min, Bac2.max, iBac.min, iBac.max, iBac2.min, iBac2.max, Pha.min, Pha.max, Pha2.min, Pha2.max)),
     type = "n",
     ylim = c(-10, 3),
     xlab = "log10 resource supply (N)",
     ylab = "log 10 normalized Biovolume")

points(log10(c(NutCon2, NutCon2)), log10(c(Nut.min, Nut.max)), col = "navy", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Bac.min, Bac.max)), col = "green2", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Bac2.min, Bac2.max)), col = "goldenrod1", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(iBac.min, iBac.max)), col = "green4", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(iBac2.min, iBac2.max)), col = "darkgoldenrod3", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Pha.min, Pha.max)), col = "brown2", cex = 0.25, pch = 20)
points(log10(c(NutCon2, NutCon2)), log10(c(Pha2.min, Pha2.max)), col = "firebrick", cex = 0.25, pch = 20)
abline(v = c(-0.25, 0, 0.25, 0.50, 0.75), h = log10(0.2010076))

legend("bottomright", 
       legend=paste(c("Nutrients", "Bacteria", "Bacteria2", "Infected", "Infected2", "Phage", "Phage2")), 
       col= c("navy", "green2",  "goldenrod1", "green4", "darkgoldenrod3", "brown2", "firebrick"), 
       lty=1, ncol=1, cex=0.5, xpd = T)

##### 'Amplitude and Nutrients' #####
TSlength <- 20000
TSanalyzed <- 0.99
log10.NutCon <- seq(log10(0.60), log10(8), length = 400)
NutCon <- 10^log10.NutCon
Amplitude <- seq(0, 0.99, length = 400)

TL <- 720 # 30 Tage
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
    
    BasicPhage_set_params(a = Amplitude[j], p = p, D = D, S = NutCon[k], y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b,  s = s, d = d, SH = SH, WP = 0.2010076)
    
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Surviver = c()
    Speccount <- 0
    
    for(l in 3:8){
      if(tail(Sim.01.Sub[,l],1) > 0.0){
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
#str(results1)

ggplot(results1, aes(Nutrients, Amplitude, fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = rainbow(7),
                       breaks = seq(0, 6, 1),
                       limits = c(0, 6)) +
  scale_x_continuous(trans = 'log10')

##### 'Nutrients and Period' #####
TSlength <- 10000
TSanalyzed <- 0.99
log10.NutCon <- seq(log10(0.60), log10(8), length = 350)
NutCon <- 10^log10.NutCon
log10.Period <- seq(log10(2*pi/(24)), log10(2*pi/(8760)), length = 350)
Period <- 10^log10.Period

y = 7.5   #3.75 #15
y1 = 4    #2 #8
y2 = 7.5  #3.75 #15
y3 = 4    #2 #8

results2 <- data.frame(Nutrients = c(),
                      Period = c(),
                      SpecNum = c(),
                      Bac = c(),
                      Bac2 = c(),
                      iBac = c(),
                      iBac2 = c(),
                      Pha = c(),
                      Pha2 = c())

for (k in 1:length(NutCon)) {
  for(j in 1:length(Period)){
    
    xStart <- c(NutCon[k], Bacteria, Bacteria2, Infected, Infected2, Phage, Phage2)
    
    BasicPhage_set_params(a = 0.9, p = Period[j], D = D, S = NutCon[k], y = y, y1 = y1, y2 = y2, y3 = y3, m = m, K = K,
                          i = i, i1 = i1, c = c, b = b,  s = s, d = d, SH = SH, WP = 0.2010076)
    
    Sim.01 <- BasicPhage_adap(xStart, TSlength)
    Sim.01.Sub <- subset(Sim.01, Time > TSlength*TSanalyzed)
    Sim.01.Sub[, 3:8][Sim.01.Sub[, 3:8] < 10e-7] = 0.0
    
    Surviver = c()
    Speccount <- 0
    
    for(l in 3:8){
      if(tail(Sim.01.Sub[,l],1) > 0.0){
        Surviver <- c(Surviver, 1)
      }else{
        Surviver <- c(Surviver, 0)  
      }
    }
    
    Speccount <- sum(Surviver)
    
    results.loop <- data.frame(Nutrients = NutCon[k],
                               Period = Period[j],
                               SpecNum = Speccount,
                               Bac = Surviver[1],
                               Bac2 = Surviver[2],
                               iBac = Surviver[3],
                               iBac2 = Surviver[4],
                               Pha = Surviver[5],
                               Pha2 = Surviver[6])
    
    results2 <- rbind(results2, results.loop)
    
    #print(paste("Nutrients = ",k,"; Period = ",j))
  }
}

ggplot(results2, aes(Nutrients, (2*pi/Period/24), fill = SpecNum))+
  geom_raster(interpolate = T) +
  scale_fill_gradientn(name = "", colors = rainbow(7),
                       breaks = seq(0, 6, 1),
                       limits = c(0, 6)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10')
