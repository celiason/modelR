(define-param no-struc? true)
(use-output-directory "output")

; geometrical input params
(define-param a .300)  ; lattice parameter
(define-param r1 .5)  ; outer radius / lattice parameter
(define-param r2 .5)  ; inner radius / lattice parameter
(define-param r3 .5)
(define r1 (* r1 a))  ; radius size
(define r2 (* r1 r2))
(define r3 (* r3 a))
(define-param nx 13)  ; number of supercells in x direction
(define-param ny 7)  ; number of supercells in y direction
(define d (* (/ a 2) (sqrt 3)))
(define dx a)  ; width of supercell
(define dy (* d 2))  ; height of supercell
(define pad 2)
(define dpml 1)
(define sx0 (* nx dx))  ; width of lattice
(define sy0 (* ny dy))  ; height of lattice

; source characteristics
(define-param fcen 2.4)
(define-param df 4)
(define-param nfreq 100)

; material properties
(define-param eps1 2.43)  ; permittivity of low index material
(define-param eps2 4)  ; permittivity of high index material
(define-param k1 0)  ; extinction coeff mat 1
(define-param k2 0)  ; extincion coeff mat 2
(define mat1 (make medium (epsilon eps1) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps1) k1)) eps1))))
(define mat2 (make medium (epsilon eps2) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps2) k2)) eps2))))

; make geometric lattice
(set! geometry-lattice (make lattice (size sx0 sy0 no-size)))

; define geometry
(if no-struc?
	(set! geometry 
   (list 
		(make block (center 0 0) 
                (size (+ sx0 (* 2 dpml)) (+ sy0 (* 2 dpml) (* 2 pad)) .01)
                (material mat1))))
	(set! geometry 
			(append
				(list
					(make block (center 0 0) (size (+ sx0 (* 2 dpml)) (+ sy0 (* 2 pad) (* 2 dpml))) (material mat1)))
				(geometric-objects-lattice-duplicates
					(list 
						(make cylinder (center 0 0) (radius r1) (height infinity) (material mat2))
						(make cylinder (center 0 0) (radius r2) (height infinity) (material air))
						(make cylinder (center (/ dx 2) (/ dy 2)) (radius r1) (height infinity) (material mat2))
						(make cylinder (center (/ dx 2) (/ dy 2)) (radius r2) (height infinity) (material air))
						(make cylinder (center (/ dx 2) (/ dy -2)) (radius r1) (height infinity) (material mat2))
						(make cylinder (center (/ dx 2) (/ dy -2)) (radius r2) (height infinity) (material air))
						(make cylinder (center (/ dx -2) (/ dy 2)) (radius r1) (height infinity) (material mat2))
						(make cylinder (center (/ dx -2) (/ dy 2)) (radius r2) (height infinity) (material air))
						(make cylinder (center (/ dx -2) (/ dy -2)) (radius r1) (height infinity) (material mat2))
						(make cylinder (center (/ dx -2) (/ dy -2)) (radius r2) (height infinity) (material air))
						)
					dx dy)
				; point defect in center
				(list 
					(make cylinder (center 0 0) (radius r3) (material mat1) (height infinity)))
				; line defect in center (missing row of melanosomes)
				; (geometric-object-duplicates (vector3 a 0) 0 12
	   ;      (make cylinder (center (+ (* -0.5 sx0) (/ dx 2)) 0) (radius r1) (height infinity)
	   ;            (material mat1)))
		)
	)
)

; expand geometry, add dpml
(set! geometry-lattice (make lattice (size (+ sx0 (* 2 dpml)) (+ sy0 (* dpml 2)) no-size)))


; source for transmission spectrum
; (set! sources (list
; 	 (make source
; 	  (src (make gaussian-src (frequency fcen) (fwidth df)))
; 	    (component Hz)
; 			(center 0 (- (* 0.5 sy0) dpml 1))
; 			(size sx0 0))))

; source for resonant mode
(set! sources (list
				(make source
					(src (make gaussian-src (frequency fcen) (fwidth df)))
					(component Hz) (center 0 0))))

; (set! pml-layers (list (make pml (thickness dpml) (direction Y))))
(set! pml-layers (list (make pml (thickness dpml))))


(define-param res 20)
(set! resolution res)
; (set! k-point (vector3 0 0 0))
; (set! ensure-periodicity true)

; transmitted flux
; (define trans
; 	(add-flux fcen df nfreq
; 	 (make flux-region
; 	 (center 0 (+ (* -0.5 sy0) dpml 1)) (size (* sx0 2) 0))))

; symmetries for Hz (TE) mode, resonance calculation
(set! symmetries
      (list (make mirror-sym (direction Y) (phase -1))
            (make mirror-sym (direction X) (phase -1))))

; run source (for transmission spectrum)
; (run-sources+
; 	(stop-when-fields-decayed 50 Hz
; 		(vector3 0 (- (* 0.5 sy0) dpml 1))
; 		1e-2)
; 	(at-beginning output-epsilon)
; 	(at-end (output-png Hz "-Zc bluered")))
; ; print out the flux spectrum
; (display-fluxes trans)

; run source (for resonant mode calculation)
(run-sources+ 400
		(at-beginning output-epsilon)
		(after-sources (harminv Hz (vector3 0) fcen df)))

; for producing visualization of field in cavity
(run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-hfield-z))

; use this code to convert fields inside cavity to a movie:
; h5topng -RZc dkbluered -C output/eps-000000.00.h5 output/hz-*.h5
; convert output/hz-*.png movie-hz.gif
