; resonant mode of hollow tube

(use-output-directory "output")

(define-param diam .136) ; diameter
(define-param cor .248) ; cortex
(define-param n_rod 2) ; index of rod
(define-param n_cor 1.56) ; index of cortex

(define-param pad 2)

(define-param dpml 1)

(define sx (+ diam cor pad (* 2 dpml)))
(define sy (+ diam (* 2 dpml)))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(define thetas (do ((vec (make-vector 91))
     (i 0 (+ i 1)))
    ((= i 91) vec)
  (vector-set! vec i i)))


; (set! default-material (make dielectric (index m)))


(define-param thickness 1) ; graded stack thickness
(define-param rug 100) ; number of slices to cut stack into

(define thetas )
theta <- seq(90, 0, length=100)/180*pi

x <- r*cos(theta) + r
y <- r*sin(theta)
# plot(y~x, type='l')
thick <- 2*y
# plot(thick~x)
mel <- thick
air <- d-thick
plot(mel/d, type='l')
lines(air/d)
lines((mel+air)/d)

# pick RIs to use
A <- d^2
A1 <- pi*r^2
A2 <- A-pi*r^2
n1 <- 2
# neff <- (A1*n1+A2*n2)/A
neff <- 1.785398
n1 <- seq(2, 1.785398, length=10)
n2 <- (neff*A-A1*n1)/A2
n1 <- round(n1, 3)
n2 <- round(n2, 3)




(let loop ((n 0))
(if (<= n 100)
(begin
    (set! geometry (append geometry (list
  (make block
     (size 50 0.1 )
        (center -6 (/ n 10) )
     (material (make dielectric (index (- LN (/ n 100) ) ) ) ) )
  (make block
        (size 50 0.1 ) 
        (center -6 (- 0 (/ n 10) ) )
     (material (make dielectric (index (- LN (/ n 100) ) ) ) ) )
   ) ) ); end geometry 
(loop (+ n 1) )
) ;end begin
   ) ; end if
) ; end let




(set! geometry (list

    (if (= 0 cor)
        (make block (center (+ (* -0.5 sx) (/ (+ dpml pad cor) 2)) 0) 
                    (size (+ dpml cor pad) sy) 
                    (material (make dielectric (index n_cor))))
				(make block (center (+ (* -0.5 sx) dpml pad (/ cor 2)) 0) 
										(size cor sy)
										(material (make dielectric (index n_cor)))))
        (make cylinder (center (+ (* -0.5 sx) dpml pad cor (/ diam 2)) 0)
        							 (height infinity)
        							 (radius (/ diam 2))
        							 (material (make dielectric (index n_rod))))))

(set! pml-layers (list (make pml (thickness dpml))))

(define-param res 30)

(set-param! resolution res)

(define-param fcen 2) ; pulse center frequency (2 = 500 nm)

(define-param df .1) ; pulse width (in frequency)

(set! sources (list
			(make source
					(src (make gaussian-src (frequency fcen) (fwidth df)))
						(component Hz)
							; (center (+ (* -0.5 sx) dpml pad (/ cor 2)) 0)
							(center (+ (* -0.5 sx) dpml 1) 0)
							(size 0 sy)
							)))

; (run-sources+ 10)
; 		(at-beginning output-epsilon)
; 			(after-sources (harminv Hz (vector3 (+ (* -0.5 sx) dpml pad (/ cor 2)) 0) fcen df)))

(run-until (/ 10 fcen) 
	(at-beginning output-epsilon)
	(at-every (/ 10 fcen 20) output-hfield-z))
