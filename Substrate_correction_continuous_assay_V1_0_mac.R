# ***Script by Rory R. Koenen, October 2024, version V1.0
# ***Correction method and models by Jan Rosing, DOI: 10.1160/TH15-04-0354
# Perform substrate correction and kinetic fitting on parabolic curves of FVIIa/TF - Xa generation and inactivation by TFPI

library(here)
require(stats)

#*******************************************
#*** section 1, data import and cleaning ***
#*******************************************
rm(list = ls()) # clear workspace
setwd(here())
date.1<-format(Sys.time(), "%Y-%m-%d")
time.1<-format(Sys.time(), "%H-%M")
filepath<-here(date.1)
dir.create(filepath)
Xa.df<- read.table(pipe("pbpaste"), sep="\t", header=T) #use this when running on mac
#Xa.df<- read.table(file="clipboard", sep="\t", header=T) #use this when running on windows
dim.1<-dim(Xa.df) #dimensions of the data
dim.2<-dim.1[2]-1 # take dimension of data frame minus the time vector.
names(Xa.df)<-make.names(colnames(Xa.df), unique = TRUE) #make names syntactically correct 
names.old <- names(Xa.df) #store the old names
names.1<-paste("XaCurve", 1:dim.2, sep="_") #make a nice list of universal names
names(Xa.df)<- c("time", names.1) #set correct labels as time and "XaCurve_n"
Xa.df<-as.data.frame(Xa.df, row.names = NULL) # convert into a dataframe
#the following lines check for empty cells, rows, columns:
if (any(is.na(Xa.df))){
Xa.df<-Xa.df[, colSums(is.na(Xa.df)) != nrow(Xa.df)] # remove columns only containing NA
Xa.df <- Xa.df[-which(!complete.cases(Xa.df)),] # remove rows containing NA
}
Xa.df<-apply(Xa.df, 2, function(x) {ifelse(x < 0, 0, x)})  # replace any negative values with zeros
dim.1<-dim(Xa.df)#update dimensions as columns might have been removed
dim.2<-dim.1[2]-1
names.old <- names.old[1:dim.1[2]] #shorten the names.old vector if columns were removed.
Xa.df<-as.data.frame(Xa.df)  # convert data back into a data frame ("apply" turned it into a matrix)
Xa_cor<-as.data.frame(Xa.df[,1:2])#take time and control curve, without inhibitor for substrate correction

#***************************************
#*** section 2, substrate correction ***
#***************************************

#The data is parabolic, so fit initial time points with parabola
initial_slope <-seq(0,0, length.out = dim.2)
x<-Xa_cor[3:13,1]
y<-Xa_cor[3:13,2] 
initial_slope<-lm(y ~ poly(x,2, raw=TRUE)) # we only take the first part (10 points) of the curve and ignore the first 3 time points because they sometimes cause trouble.

blank<- initial_slope$coefficients[1] #blank
slope.1<-initial_slope$coefficients[2] #slope
slope.2<-initial_slope$coefficients[3] #square
initial_slope$coefficients #show coefficients

#plot curve and initial fit, you can readily see the effect of substrate consumption 
Xa_cor$linfit<- (blank+(Xa_cor[,1]*slope.1)+(Xa_cor[,1]^2*slope.2))
x <- Xa_cor$time
y <- Xa_cor$XaCurve_1
l <- Xa_cor$linfit
plot(x,y,main = "Xa_curve + slope_fit", pch="°", xlab = "time", ylab = "mOD", xlim = c(0,max(x)))
lines(x,l, col="red")

n<- length(Xa_cor[,1]) # set length of vector

#fit the original curve over ca. 50% of its length using a 6-degree polynomial 
l<-as.integer(n*0.5) #this value, which should be 1 (entire curve) or less is critical. It should neither be too high or too low and can be empirically determined.
x<-Xa_cor$time[1:l] #set x values
y<-Xa_cor$XaCurve_1[1:l] #set y values
polyfit.1<-lm(y ~ poly(x,6, raw=TRUE))
sum.fit.1<-summary(polyfit.1)
params.1<-as.numeric(c(coefficients(polyfit.1),0,0)) #store model parameters in a vector to use for the calculation below
coef.diff<-abs(slope.2-params.1[2])
coef.diff
#plot(polyfit.1) #unhash to plot model details
params.1

#predict using the parameters from the 6-degree polynomial fit (polyfit.1)
#below is the function "poly_cal" of a 6th degree polynomial to be fitted:
poly_cal<-function(x){
  a<-params.1[1]
  b<-params.1[2]
  c<-params.1[3]
  d<-params.1[4]
  e<-params.1[5]
  f<-params.1[6]
  g<-params.1[7]
  
  y=a+(b*x)+(c*x^2)+(d*x^3)+(e*x^4)+(f*x^5)+(g*x^6)
  return(y)
}

Xa_fit.df<-as.data.frame(cbind(Xa_cor$time[1:l],Xa_cor$XaCurve_1[1:l]))
times<- Xa.df$time[1:l]
datafit<-sapply(times, poly_cal) #the sapply command runs the polynomial function over all values
Xa_fit.df$datafit<-datafit
names(Xa_fit.df)<- c("t", "data", "fit")

#make a plot of the fitted uncorrected Xa-curve, the 2 curves should be completely overlapping
Xa_curve<-as.data.frame(Xa_fit.df)
names(Xa_curve)<- c("t", "data", "fit")
x<-Xa_curve$t
y<-Xa_curve$data
f<-Xa_curve$fit
plot(x,y,main = "Xa_curve + 6th degree fit", pch="°", xlab = "time", ylab = "mOD", xlim = c(0,20), ylim = c(0,0.9))
lines(x,f, col="red")

#now calculate the residuals from the linear fit of the initial slope and the fitted progress curve and plot the result
Xa_fit.df$linfit<- (Xa_cor$linfit[1:l])
Xa_fit.df$res<-Xa_fit.df$linfit-Xa_fit.df$fit

x<-Xa_fit.df$data
y<-Xa_fit.df$res
ylim.1<-1.1*max(Xa_fit.df$res)
ylim.1b<- 1.1*min(Xa_fit.df$res)
xlim.1<-1.2*max(Xa_fit.df$data)
plot(x,y,main = "residuals vs corrected mOD", pch="°", xlab = "mOD", ylab = "residuals", xlim = c(0,xlim.1), ylim = c(ylim.1b,ylim.1)) #remove the bracket and the hash to set plot limits
#as one can see, the residuals increase as mOD increases, due to substrate consumption.

#fit the curve residuals versus mOD using a 6th degree polynomial
x<-Xa_fit.df$data[1:l]
y<-Xa_fit.df$res[1:l]
polyfit.2<-lm(y ~ poly(x,6, raw=TRUE))
sum.fit.2<-summary(polyfit.2)
# #unhash to plot model details

#predict residuals using a defined 6-degree polynomial function
params.2<-as.numeric(coefficients(polyfit.2))
params.2
sixt_cal<-function(x){
  a<-0
  b<-params.2[2]
  c<-params.2[3]
  d<-params.2[4]
  e<-params.2[5]
  f<-params.2[6]
  g<-params.2[7]
  y=a+(b*x)+(c*x^2)+(d*x^3)+(e*x^4)+(f*x^5)+(g*x^6)
  return(y)
}

#predict data using 6-degree polynomial
corr_mOD<- Xa_fit.df$data
Xa_fit_val<-sapply(corr_mOD,sixt_cal) #the sapply command runs the polynomial function over all values. Note: "predict()" returns erratic results
Xa_fit.df$val<-Xa_fit_val
Xa_fit.df$cor<-Xa_fit.df$val+Xa_fit.df$data
res.cor<-Xa_fit.df$cor-Xa_fit.df$linfit
sumsq_rescor<-sum( (res.cor - mean(res.cor) )^2 )

x<-Xa_fit.df$data
y<-Xa_fit.df$res
f<-Xa_fit.df$val
ylim.1<-0.5*max(Xa_fit.df$res)
ylim.1b<- -min(Xa_fit.df$res)
xlim.1<-0.8*max(Xa_fit.df$data)
plot(x,y,main = "residuals  + 6th degree fit", pch="°", xlab = "mOD", ylab = "residuals")#, xlim = c(0,xlim.1), ylim = c(-1,ylim.1))
lines(x,f, col="blue")

#plot to see if the corrected curve now lies over the ideal fit

x<-Xa_fit.df$t
y<-Xa_fit.df$cor
f<-Xa_fit.df$linfit
ylim.1<-1.1*max(Xa_fit.df$linfit)
xlim.1<-1.1*max(Xa_fit.df$t)
plot(x,y,main = "cons. corrected mOD + initial fit", pch="°", xlab = "time", ylab = "mOD", xlim = c(0,xlim.1), ylim = c(0,ylim.1))
lines(x,f, col="blue")

params<-rbind(params.1[1:7],params.2)
colnames(params)<-c("intercept","X1","X2","X3","X4","X5","X6")
filename<-paste(filepath,"/","5_model_parameters",time.1,".csv", sep = "")
write.csv(params, file = filename) #saving a csv with the fitted parameters
rm(Xa_cor) #clean up obsolete dataframes
rm(Xa_curve) #clean up obsolete dataframes

#************************************************
#*** section 3, correcting experimental data  ***
#************************************************

#import the data to be corrected for substrate consumption from a new Excel file
TFPI_raw.df<-Xa.df

#fit all data in the imported dataset
TFPI_cor.df<-as.data.frame(TFPI_raw.df[,1])
for (i in 2:dim.1[2]){
corr_mOD<- TFPI_raw.df[,i]-blank
corr_mOD<-sapply(corr_mOD, function(x) {ifelse(x < 0, 0, x)}) #make any negatives zero
TFPI_cor.1<-ifelse(corr_mOD>blank,sapply(corr_mOD,sixt_cal),0) #only correct values that are over the blank
TFPI_cor.df[,i]<-corr_mOD+TFPI_cor.1
}

names.3<-paste("TFPI_cor", 1:dim.2, sep="_") #make a nice list of universal names
names(TFPI_cor.df)<- c("time", names.3) #set new labels as time and "TFPI_fit_n"

for(i in 2:dim.1[2]){
  x<-TFPI_raw.df[,1]
  y<-TFPI_raw.df[,i]
  xc<-TFPI_cor.df[,1]
  yc<-TFPI_cor.df[,i]
  ylim.1<-1.2*max(TFPI_cor.df[,i])
  plot(x,y,main = paste("corrected vs raw curves","\n",names.old[i]), xlab = "time", ylab = "mOD", ylim = c(0,ylim.1))
  par(new=TRUE) 
  lines(xc,yc,col="purple")
}

#fit the corrected curves with a 2nd degree polynomial

dim.1<-dim(TFPI_cor.df) #dimensions of the data
dim.2<-dim.1[2]-1 # take dimension of data frame minus the time vector.
names.2<-paste("TFPI_raw", 1:dim.2, sep="_") #make a nice list of universal names
names(TFPI_cor.df)<- c("time", names.2) #set correct labels as time and "TFPI_curve_n"
test_raw.df<-as.data.frame(TFPI_cor.df, row.names = NULL) # convert into a dataframe
test_raw.df[1,]<-0

square <-seq(0,0, length.out = dim.2)
slope <-seq(0,0, length.out = dim.2)
intercept<-seq(0,0, length.out = dim.2)
corcoefs<-NULL
for(i in 2:dim.1[2]){
  x<-test_raw.df[2:22,1]#the quality of the fit depends on the part you take in, here it is ca. 5 minutes
  y<-test_raw.df[2:22,i] 
  linearfit<-lm(y ~ poly(x,2, raw=TRUE)) #here is the 2nd degree polynomial
  coefs.1<-coefficients(linearfit)
  square[i]<-coefs.1[3]
  slope[i]<-coefs.1[2]
  intercept[i]<-coefs.1[1]
}

intercept<-intercept[-1]#remove the zeros from the parameters
slope<-slope[-1]
square<-square[-1]
TFPI_fit.df<-seq(0,0, length.out = dim.1[1])
TFPI_fit.df<-(data.frame(cbind(TFPI_fit.df)))
TFPI_fit.df[,1]<-test_raw.df[,1]
x<-test_raw.df[,1]
for(i in 2:dim.1[2]){
TFPI_fit.df[,i]<-intercept[i-1]+ (x*slope[i-1])+(x^2*square[i-1])
}

#the resulting curves are those that you would get if there was no substrate consumption effect (and assuming that there is no further slow-tight rearrangement, which is probably the case)
for(i in 2:dim.1[2]){
  x<-TFPI_cor.df[,1]
  y<-TFPI_cor.df[,i]
  xc<-TFPI_fit.df[,1]
  yc<-TFPI_fit.df[,i]
  ylim.1<-1.2*max(TFPI_cor.df[,i])
  plot(x,y,main = paste("fitted vs corrected curves","\n",names.old[i]), xlab = "time", ylab = "mOD", ylim = c(0,ylim.1))
  par(new=TRUE) 
  lines(xc,yc,col="green")
}

#define a simple numeric 1st derivative function
first.deriv <- function(x,y){
  d<-(diff(y)/diff(x))
  return(d)
}

#take first derivative of fitted inactivation curves
#first test for "double zeros" in the time vector, which makes first.deriv() return Inf. If "TRUE" rows are removed.
if(min(diff(TFPI_fit.df[,1]))==0){
  TFPI_fit.df <- TFPI_fit.df[-which(diff(TFPI_fit.df[,1])==0),]
}

TFPI_fitderiv.df<-as.data.frame(TFPI_fit.df[-1,1])
for (i in 2:dim.1[2]){
  y<-TFPI_fit.df[,i]
  x<-TFPI_fit.df[,1]
  TFPI.fitderiv.1<- first.deriv(x,y)
  TFPI_fitderiv.df[,i]<-TFPI.fitderiv.1
  color.1<-as.hexmode(5*i)
  plot(TFPI_fitderiv.df[,1],TFPI_fitderiv.df[,i], main = "Xa activity",type = "l", col=color.1 ,xlab = "time", ylab = "%Xa",xlim= c(0,20), ylim = c(0,0.5))
  par(new=TRUE)
}

#calculate the slopes of the Xa curves
dim.3 <- dim(TFPI_fitderiv.df)[2]
slopes_Xa <- seq(0,0,length.out = dim.3)
for (i in 1:dim.3){
y<-TFPI_fitderiv.df[,i]
x<-TFPI_fitderiv.df[,1]
lincof<-lm(y~x)
slopes_Xa[i] <- lincof$coefficients[2]
}

slopes_Xa <- slopes_Xa[-1]
percentage_Xa <- 100*(slopes_Xa/slopes_Xa[1])#calculate as percentage of FXa w/o TFPI

summary_Xa.df <- as.data.frame(cbind(slopes_Xa, percentage_Xa))

#save the datafiles
names<-gsub("[X]","",names.old)#remove the "X" before the column titles, may yield unwanted resuls
names(TFPI_cor.df)<-names
filename.1<-paste(filepath,"/","6_corrected FVIIa curves",time.1,".csv", sep = "")#datafile with substrate-corrected curves
write.csv(TFPI_cor.df, file=filename.1, row.names = FALSE)
rm(TFPI_raw.df) #clean up obsolete dataframes

names(TFPI_fit.df)<-names
filename.1<-paste(filepath,"/","7_Fitted FVIIa curves",time.1,".csv", sep = "")#datafile with 2nd degree poly-fitted curves
write.csv(TFPI_fit.df, file=filename.1, row.names = FALSE)

names(TFPI_fitderiv.df)<-names
filename.1<-paste(filepath,"/","8_derivatives FVIIa curves",time.1,".csv", sep = "")#datafile with 1st derivates of fitted curves, convertible to nM Xa 
write.csv(TFPI_fitderiv.df, file=filename.1, row.names = FALSE)

rownames(summary_Xa.df)<-names[-1]
filename.1<-paste(filepath,"/","9_slopes of Xa curves",time.1,".csv", sep = "")#datafile with the slopes, OD/min, and percentage of initial Xa activity 
write.csv(summary_Xa.df, file=filename.1, row.names = TRUE)






