# Bioscreen functions

#===== Function =======
# Spline function (traditional approach) ----------------------------------

# "nderiv" is a function that takes a model "fit" and a vector of time "x" and gives the derivative of the model at each point x
# it calculates the derivative by taking a very small step (eps, 1e-5) in both directions from each point
# and taking the difference in predicted value for the two points (x + eps, and x - eps)
# it returns a vector with derivatives for all n elements in vector x
nderiv <- function(fit, x, eps=1e-5){
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)} #

# "spline.slope" is a function that takes in a time vector, OD measurement, a count n (nr of time points), a distance eps, and span.
# it uses its input to form a model ("fit"), which is the argument for the nderiv function.
# the model is a polynomial surface, loess(), of the growthcurve, log(OD)
# span is a parameter for the smoothness of the surface (higher value = fewer points = rougher surface)
# spline.slope calls nderiv and calculates the maximum derivative.

# Regarding the span parameter
# In Gerstein Biology Letters from 2012, they used a span of 0.1 for 48 hours with measurements every 30 minutes (a total of 96 data points).
# Basically, for each point, 10 % of the data was used, which means 9.6 data points, which will encompass 9.6/2 = 4.8 hours.
# So basically, just divide 4.8 with the number of hours, to standardize it according to gerstein procedure

# Regarding the n parameter
# In Gerstein Biology Letters from 2012, they used n = 101. 48 hours / 101 is a little more than two points per hour. 
# I'd say to get the number of points by multiplying the hours by 2.

spline.slope <- function(time, OD, n=max(time)*2, eps=1e-5, span=4.8/max(time)){ 
  max(nderiv(loess(log(OD) ~ time, degree=1, span=span), seq(min(time), max(time), length=n)), na.rm=TRUE)
}

spline <-function(time, OD, n=length(time), eps=1e-5, span=0.2){
  span = 4.8/max(time)
  n = max(time)*2
  maxr_time <- time[which(nderiv(loess(log(OD) ~ time, degree=1, span=span), seq(min(time), max(time), length=n)) == 
                            max(nderiv(loess(log(OD) ~ time, degree=1, span=span), seq(min(time), max(time), length=n)), na.rm = T))*2]
  return(maxr_time)     
}
