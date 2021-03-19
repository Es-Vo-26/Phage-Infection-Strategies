# title: "Bacteria only, without phage infection"
# ODE model

rm(list = ls())

##### 'libraries' #####
library(odeintr)

##### 'start values' #####
Resource <- 1.3
Bacteria <- 0.83
Bacteria2 <- 0.83
xStart <- c(Resource, Bacteria, Bacteria2)

##### 'parameter' #####

##### 'densities' 
# all densities are normalized by half. sat density 5.2333e06 µm^3 mL^-03'
S <- Resource  # supply concentration
K <- 1         # halfsaturation density

##### 'rates' 
# normalized by metabolic rate 0.002365714285714 hr^-1
D <- 17.36    # turnover rate, mean nutrient input
m <- 1        # metabolic rate

##### 'conversion rates'
y <- 3.5        # metabolic scaling constant for bacterial growth
y1 <- 2.5       # metabolic scaling constant for inf. bacterial growth
a <- 0 #0.5     # amplitude, variation of nutrient fluctuation

TL <- 168       # time of one resource turnover
# 24 hr
# 168
# 720
# 8760

p <- 2*pi/TL    # oscillation period

##### 'Model of phage-bacteria-interaction' #####
BasicBacteria.sys <- '
double XX0 = x[0];
double XX1 = x[1];
double XX2 = x[2];

if (XX1 <= pow(10,-7)){
double XX1 = 0.0;
}
if (XX2 <= pow(10,-7)){
double XX2 = 0.0;
}

dxdt[0] = (a*D*sin(p*t)+D)*(S-XX0) - y*m*XX0/(K+XX0)*XX1 - y1*m*XX0/(K+XX0)*XX2;

if(XX1 <= pow(10, -7)){
  dxdt[1] = 0.0;
}
else{
  dxdt[1] = y*m*XX0/(K+XX0)*XX1 - m*XX1;
}

if(XX2 <= pow(10, -7)){
  dxdt[2] = 0.0;
}
else{
  dxdt[2] = y1*m*XX0/(K+XX0)*XX2 - m*XX2;
}
' 
compile_sys(name = "BasicBacteria",
            sys = BasicBacteria.sys,
            pars = c("a", "p", "D", "S", "y", "y1", 
                     "K", "m"),
            method = "rk54_a")

BasicBacteria_set_params(a = a, p = p, D = D, S = S, y = y, y1 = y1, K = K, m = m)

TS.01 <- BasicBacteria_adap(xStart, 500)
TS.01[TS.01 < 10e-7] = 0.0

##### 'Plot Data' #####
par(mfrow = c(1,2))
 
plot(c(TS.01$Time, TS.01$Time, TS.01$Time), 
     log10(c(TS.01$X1, TS.01$X2, TS.01$X3)),
     type = "n",
     xlab = "Time",
     ylab = "Log10 norm. biovolume")
lines(TS.01$Time, log10(TS.01$X1), col = "navy", lwd = 2)
lines(TS.01$Time, log10(TS.01$X2), col = "green2", lwd = 2)
lines(TS.01$Time, log10(TS.01$X3), col = "goldenrod1", lwd = 2)

plot(c(TS.01$Time, TS.01$Time, TS.01$Time), 
     c(TS.01$X1, TS.01$X2, TS.01$X3),
     type = "n",
     xlab = "Time",
     ylab = "norm. biovolume")
lines(TS.01$Time, TS.01$X1, col = "navy", lwd = 2)
lines(TS.01$Time, TS.01$X2, col = "green2", lwd = 2)
lines(TS.01$Time, TS.01$X3, col = "goldenrod1", lwd = 2)

legend("bottomright", 
       legend=paste(c("Nutrients", "Bacteria", "Bacteria2")), 
       col= c("blue", "green2", "green4"), 
       lty=1, ncol=1, cex=0.5, xpd = T)