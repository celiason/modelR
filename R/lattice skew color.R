#snippet to calculate multiply Bragg diffracted peaks in hexagonal (or centered rectangular) lattices (i.e. a1!=a2) given lattice constants a1 and a2 (a1 parallel to surface)
foo <- function(a1=500, a2=500) {

#a1=535  # lattice constant parallel to surface
#a2=632  # lattice constant 2

d1=a1*sin(acos(a1/2/a2))  # distance b/w Bragg planes 
d2=sqrt(a2^2-(a1/2)^2)

theta2=acos(d1/a1)
theta3=acos(d1/a2)

#theta2
#theta3

n.eff=1.73

lambda1=2*n.eff*cos(theta2)*d1/2
lambda2=2*n.eff*cos(theta3)*d2/2

return(lambda1)
#lambda2

}
