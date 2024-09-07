; resonant mode of hollow tube

(use-output-directory "output")

(define-param n 2) ; index of shell = melanin
(define-param w 1) ; width of melanin shell in um
(define-param r 1) ; radius of air space in center
(define-param mat 3.1) ; index of matrix

(define-param pad 1)
(define-param dpml 1)

(define sxy (* 2 (+ r w pad dpml))) ; size of cell
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

; (set! default-material (make dielectric (index m)))

(set! geometry (list
				(make block (center 0 0) (size sxy sxy infinity) (material (make dielectric (index mat))))
        (make cylinder (center 0 0) (height infinity)
					(radius (+ r w)) (material (make dielectric (index n))))
				(make cylinder (center 0 0) (height infinity)
					(radius r) (material air))))
(set! pml-layers (list (make pml (thickness dpml))))

(set-param! resolution 30)

(define-param fcen .15) ; pulse center frequency
(define-param df .1) ; pulse width (in frequency)
(set! sources (list
				(make source
					(src (make gaussian-src (frequency fcen) (fwidth df)))
					(component Ez) (center (+ r 0.1) 0))))

(run-sources+ 300
		(at-beginning output-epsilon)
		(after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))

; (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-efield-z))
