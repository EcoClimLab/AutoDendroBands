########################################
# Purpose: Runing a wavelet anaylsis
# Data: TreeHugger raw data: object called x with at least Date, Time, mm, Ch2 and Ch4 columns
# Date: 8/7/2017
# Author: Valentine Herrmann, HerrmannV@si.edu
# Reference paper : Herrmann et al. 2016. Tree Circumference Dynamics in Four Forests Characterized Using Automated Dendrometer Bands. PlosOne
#########################################


# Load libraries ####
library(WaveletComp)
library(circular)


# Load raw data ####
load("data_example.RData")

# Create TimeStamp ####
x$TimeStamp <- as.POSIXct(paste(x$Date, x$Time, sep = " "), format = "%d.%m.%Y %H:%M:%S")
x$Date <- as.Date(as.character(x$Date), format =  "%d.%m.%Y")

# Calculate Band Temperature (NB: In our paper we used an aggregated temerature record, using measurements from all TreeHuggers at one site)

x$Band.Temp <- (1 / 298.15+( 1 / 3974 ) * log((x$Ch2*10000/(x$Ch4-x$Ch2)) / 10000 ))^ (-1)  -273.15

# Calculate Spline and Residuals ####
lin.int <- as.data.frame(approx(x = x$TimeStamp, y = x$mm, xout = x$TimeStamp, rule = 2, method = "linear"))
colnames(lin.int) <- c("TimeStamp", "mm.int")


smoo <- smooth.spline(x = lin.int$TimeStamp, y = lin.int$mm.int, df = length(unique(x$Date)))

x$Spline <- ifelse(is.na(x$mm), NA, predict(smoo, as.numeric(x$TimeStamp))$y)
x$Residuals <- x$mm - x$Spline



# Wavelet analysis ####

## NB: In our paper, we removed rainy days for the wavelet analysis.

## On all days ####

x <- x[x$Date %in% unique(x$Date)[tapply(x$Residuals, x$Date, function(x) sum(!is.na(x)) == 96)],] # Filter for complete days
  
if (any(!is.na(x$Residuals))){
  ts.x <- ts(x$Residuals, star = 1, frequency = (24 * 60/15)) # 96 measurements per day, 1 unit = 1 day
  
  
  ts.y <- ts(x$Band.Temp, star = 1, frequency = (24 * 60/15)) # 96 measurements per day, 1 unit = 1 day but NA every 30 minutes because really we have one measurement per half-hour. 
  
  
  my.date = x$TimeStamp
  my.data = data.frame(date = my.date, x = as.numeric(ts.x), y = as.numeric(ts.y))
  
  
  my.wc = analyze.coherency(my.data, my.pair = c("x","y"),
                            loess.span = 1,
                            dt = 1/(24 * 60/15), # dj = 1/20,
                            lowerPeriod = 24/24,
                            upperPeriod = 25/24,
                            make.pval = T, n.sim = 100)
  
  
  Period <- which((my.wc$Period - 1) == min(my.wc$Period - 1))
  Coherence <- my.wc$Coherence[Period,]
  Angle <- my.wc$Angle[Period,]
  
  Coherence.70 <- Coherence > 0.70
  Mean.Phase.Angle <- mean.circular(mean.circular(as.circular(Angle[Coherence.70], type= "angles", units = "radians", modulo = 'asis', zero = 0,  rotation = 'counter', template = 'none')))
  Mean.Coherence <- mean(Coherence)
  Proportion.days.coh.70 <- sum(tapply(Coherence, x$Date, mean) > 0.70) / length(unique(x$Date))
  
  }


Mean.Phase.Angle <- as.circular(Mean.Phase.Angle, type= "angles", units = "radians", modulo = 'asis', zero = 0,  rotation = 'counter', template = 'none')



### circular plot

par(oma = c(0,0,0,0), mar = c(0,0,0,0))
res <- rose.diag(circular(0), bins = 4, prop = 0.75, col = "indianred1", cex = 1.8, shrink = 1.1)
res <- rose.diag(circular(pi/2), bins = 4, prop = 0.75, col = "indianred3", add = T, axes =  F)
res <- rose.diag(circular(pi), bins = 4, prop = 0.75, col = "palegreen2", add = T, axes = F)
res <- rose.diag(circular(3*pi/2), bins = 4, prop = 0.75, col = "palegreen3", add = T, axes = F)


text(0.3,0.3, "In phase\ntree leading", cex = 1.5)
text(-0.3,-0.3, "Out of phase\ntree leading", cex = 1.5)
text(-0.3,0.3, "Out of phase\ntree lagging", cex = 1.5)
text(0.3,-0.3, "In phase\ntree lagging", cex = 1.5)

points.circular(plot.info = res, Mean.Phase.Angle,  stack=FALSE, col = "chocolate", pch = 20, cex = 2)




## looking at each day separately ####

x <- x[x$Date %in% unique(x$Date)[tapply(x$Residuals,x$Date, function(x) sum(!is.na(x)) >= 72)],]  #Filter for 75% complete days

if(nrow(x) >0){
  m <- data.frame(Date = as.character(unique(x$Date)), Mean.Phase.Angle = NA, Mean.Coherence = NA, Proportion.coh.70 = NA)
  
  for(d in as.character(unique(x$Date))){
    
    print(d)
    
    x1 <- x[x$Date == d,]
    
    if (any(!is.na(x$Residuals))){
      
      ts.x <- ts(x1$Residuals, star = 1, frequency = (24 * 60/15)) # 96 measurements per day, 1 unit = 1 day
      
      
      ts.y <- ts(x1$Band.Temp, star = 1, frequency = (24 * 60/15)) # 96 measurements per day, 1 unit = 1 day but NA every 30 minutes because really we have one measurement per half-hour. 
      
      
      my.date = x1$TimeStamp
      my.data = data.frame(date = my.date, x = as.numeric(ts.x), y = as.numeric(ts.y))
      
      
      my.wc = analyze.coherency(my.data, my.pair = c("x","y"),
                                loess.span = 1,
                                dt = 1/(24 * 60/15), # dj = 1/20,
                                lowerPeriod = 24/24,
                                upperPeriod = 25/24,
                                make.pval = T, n.sim = 100)
      
      
      #       my.wc$Angle " This equals the difference of individual phases, Phase.x 􀀀 Phase.y, when converted to an angle in the interval [-pi; pi]. An absolute value less (larger) than  =2 indicates that the two series move in phase (anti-phase, respectively) referring to the instantaneous time as time origin and at the frequency (or period) in question, while the sign of the phase difference shows which series is the leading one in this relationship. Figure 2 (in the style of a diagram by Aguiar-Conraria and Soares [2]) illustrates the range of possible phase differences and their interpretation. In line with this style, phase differences are displayed as arrows in the image plot of cross-wavelet power."
      
      Period <- which((my.wc$Period - 1) == min(my.wc$Period - 1))
      Coherence <- my.wc$Coherence[Period,]
      Angle <- my.wc$Angle[Period,]
      
      Coherence.70 <- all(Coherence > 0.70)
      Mean.Phase.Angle <- ifelse(Coherence.70, mean.circular(mean.circular(as.circular(Angle, type= "angles", units = "radians", modulo = 'asis', zero = 0,  rotation = 'counter', template = 'none'))), NA)
      Mean.Coherence <- mean(Coherence)
      Proportion.coh.70 <- sum(Coherence > 0.70) / length(Coherence)
      
      m[m$Date == d, ]$Mean.Phase.Angle <- Mean.Phase.Angle
      m[m$Date == d, ]$Mean.Coherence <- Mean.Coherence
      m[m$Date == d, ]$Proportion.coh.70 <- Proportion.coh.70
    }else{
      m[m$Date == d, ]$Mean.Phase.Angle <- NA
      m[m$Date == d, ]$Mean.Coherence <- NA
      m[m$Date == d, ]$Proportion.coh.70 <- NA
    }
    
    }
}  



m$Mean.Phase.Angle <- as.circular(m$Mean.Phase.Angle, type= "angles", units = "radians", modulo = 'asis', zero = 0,  rotation = 'counter', template = 'none')


### circular histogram
plot(m$Mean.Phase.Angle,  stack = FALSE, pch = 19, axes = F, ticks = F, bins = 9, cex = 2)
rose.diag(m$Mean.Phase.Angle, bins = 9, tcl.text = .3, cex = 1.5, prop = .9, ticks = F, add = T)
legend("topleft", paste("n =", sum(!is.na(m$Mean.Phase.Angle)), "days"), bty = "n")
