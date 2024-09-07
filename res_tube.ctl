; resonant mode of hollow tube

(define-param n 2) ; index of shell = melanin
(define-param w .135) ; width of melanin shell in um
(define-param r .135) ; radius of air space in center

(define-param pad 4)
(define-param dpml 1)

(define sxy (* 2 (+ r w pad dpml))) ; size of cell
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

(set! geometry (list
				(make cylinder (center 0 0) (height infinity)
					(radius (+ r w)) (material (make dielectric (index n))))
				(make cylinder (center 0 0) (height infinity)
					(radius r) (material air))))
(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 40)

(define-param fcen 2.4) ; pulse center frequency
(define-param df 2) ; pulse width (in frequency)
(set! sources (list
				(make source
					(src (make gaussian-src (frequency fcen) (fwidth df)))
					(component Ez) (center (+ r 0.05) 0))))

(run-sources+ 300
		(at-beginning output-epsilon)
		(after-sources (harminv Ez (vector3 (+ r 0.05)) fcen df)))
