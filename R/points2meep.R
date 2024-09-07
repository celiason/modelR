points2meep <- function(x, y, filename="points.txt", diam=NULL) {

x <- x - min(x)
y <- y - min(y)

#center on x,y = 0
x <- x - sum(range(x))/2
y <- y - sum(range(y))/2

#export numbers to string
if (is.null(diam)) 
  
  write(c("xx", gsub(",", "", toString(x)), "yy", gsub(",", "", toString(y)), filename))

  write(c("xx", gsub(",", "", toString(x)), "yy", gsub(",", "", toString(y)), "rr", gsub(",", "", toString(diam/2))), filename)

list(range(x), range(y))

}
