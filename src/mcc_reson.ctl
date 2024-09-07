; resonant mode of hollow tube

(use-output-directory "output")

(define-param diam .136) ; diameter
(define-param cor .248) ; cortex
(define-param n_rod 2) ; index of rod
(define-param n_cor 1.56) ; index of cortex

(define-param pad 2)

(define-param dpml 1)

(define sx (+ diam cor pad (* 2 dpml)))
; (define sy (+ (* 3 diam) (* 2 dpml)))
(define sy (+ (* 3 diam) (* 2 dpml)))

(set! geometry-lattice (make lattice (size sx sy no-size)))

; (set! default-material (make dielectric (index m)))

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
        							 (material (make dielectric (index n_rod))))
        ; (make cylinder (center (+ (* -0.5 sx) dpml pad cor (/ diam 2)) (* -1 diam))
        ;                (height infinity)
        ;                (radius (/ diam 2))
        ;                (material (make dielectric (index n_rod))))
        ; (make cylinder (center (+ (* -0.5 sx) dpml pad cor (/ diam 2)) diam)
        ;                (height infinity)
        ;                (radius (/ diam 2))
        ;                (material (make dielectric (index n_rod))))
        ))

(set! pml-layers (list (make pml (thickness dpml))))
; (set! pml-layers (list (make pml (thickness dpml) (direction X))))

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

(run-sources+ 100
		(at-beginning output-epsilon)
			(after-sources (harminv Hz (vector3 (+ (* -0.5 sx) dpml pad (/ cor 2)) 0) fcen df)))

(run-until (/ 20 fcen) 
	(at-beginning output-epsilon)
	(at-every (/ 10 fcen 20) output-hfield-z))
