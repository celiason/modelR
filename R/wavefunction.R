# wave equation
amp <- function(lambda=500, x=-500:500, t=0, amp0=1, n=1) {
	f <- 1/lambda
	T <- lambda
	w <- 2*pi/T
	c <- lambda/T
	k <- (2*pi*n)/(lambda*c)
	amp <- amp0 * cos(k*x-c*t)
	amp
}

# effect of time
plot(amp(lambda=500, x=-500:500, t=0), type='l')
lines(amp(lambda=500, x=-500:500, t=1), type='l', col="red")
lines(amp(lambda=500, x=-500:500, t=2), type='l', col="blue")

# effect of time (animated version)
times <- seq(0, 20, length=100)
for (i in seq_along(times)) {
	plot(amp(lambda=500, x=-500:500, t=times[i]), type='l')
}

# effective of changing refractive index (wave slows down)
# low RI
plot(amp(lambda=500, x=-500:500, t=0)~c(-500:500), type='l', col="darkblue", xlim=c(-100, 100))
lines(amp(lambda=500, x=-500:500, t=1)~c(-500:500), col="lightblue")
# hi RI
plot(amp(lambda=500, x=-500:500, t=0, n=10)~c(-500:500), type='l', col="darkblue", xlim=c(-100, 100))
lines(amp(lambda=500, x=-500:500, t=1, n=10)~c(-500:500), col="lightblue")
