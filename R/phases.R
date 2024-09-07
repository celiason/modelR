#THIS FUNCTION TAKES REFRACTIVE INDICES OF THREE DIFFERENT MATERIALS, THICKNESS OF FILM, WAVELENGTH, AND ANGLE AND GIVES REFLECTANCE VERSUS THETA AND PHASE OF REFLECTANCE
phases <- function(n1, n2, n3, d, lambda, theta1=seq(0,90,length=91)){

theta1 <- theta1/180*pi
theta2 <- asin(n1*sin(theta1)/n2)
theta3 <- asin(n2*sin(theta2)/n3)

r1 <- (n1*cos(theta2)-n2*cos(theta1))/(n1*cos(theta2)+n2*cos(theta1))

r2 <- (n2*cos(theta3)-n3*cos(theta2))/(n2*cos(theta3)+n3*cos(theta2))

phi <- 4*pi*n2*d*cos(theta2)/lambda

phase <- atan(r2*sin(phi)/(r1+r2*cos(phi))) - atan(r1*r2*sin(phi)/(1+r1*r2*cos(phi)))

phase <- cumsum(abs(c(0,diff(phase))))

rt <- abs(sqrt((r1^2+r2^2+2*r1*r2*cos(phi))/(1+r1^2*r2^2+2*r1*r2*cos(phi))))

par(mfrow=c(2,2))
plot(r1~I(theta1/pi*180), type='l')
plot(r2~I(theta1/pi*180), type='l')
plot(rt~I(theta1/pi*180), type='l')
plot(phase~I(theta1/pi*180), type='l')

}