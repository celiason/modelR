# calculate thin film thickness given a set of interference fringes starting 
# at the leftmost peak (l1) and ending with rightmost peak (l2). dn is the #
# of peaks from l1 to l2 (not counting l1, but counting l2). n is the refr.
# index of the film
d <- function(l1, l2, n=1.56, dn, th=0) {
  (dn*l1*l2)/(2*n*(l2-l1)*cos(th))
}

# d(525,589, dn=1)
