; script to compute reflectance spectrum of a lattice at a given rotation
; non-periodic boundary conditions w/dpml layers
; in progress

(define-param no-struc? true)
(use-output-directory "output")
(set! eps-averaging? true)

; geometrical input params
(define-param rot 0)  ; angle to rotate lattice in degrees
(define rot (* (/ rot 180) pi))  ; rotation angle in radians
(define-param a1 .100)  ; lattice constant 1 (parallel to surface)
(define-param a2 .100)  ; lattice constant 2
(define-param ra .5)  ; radius / lattice constant (a1)
(define-param r2 .5) ; radius of air core / radius of particle
(define r (* ra a1))  ; radius size
(define r2 (* r2 r))
(define-param nx 3)  ; number of supercells in x direction
(define-param ny 3)  ; number of supercells in y direction
(define d (sqrt (- (* a2 a2) (* (/ a1 2) (/ a1 2)))))  ; spacing b/w Bragg planes
(define dx a1)  ; width of supercell
(define dy (* d 2))  ; height of supercell
(define sx0 (* nx dx))  ; width of lattice
(define sy0 (* ny dy))  ; height of lattice
(define dpml 1)
;(define sx (+ 2 sx0))
;(define sy (+ 7 sy0))

; source characteristics
(define fcen 2.4)
(define df 4)
(define nfreq 100)

; material properties
(define-param eps1 2.43)  ; permittivity of low index material
(define-param eps2 4.00)  ; permittivity of high index material
(define-param k1 0.01)  ; extinction coeff mat 1
(define-param k2 0.01)  ; extincion coeff mat 2
(define mat1 (make medium (epsilon eps1) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps1) k1)) eps1))))
(define mat2 (make medium (epsilon eps2) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps2) k2)) eps2))))

; make geometric lattice
(set! geometry-lattice (make lattice (size sx0 sy0 no-size)))

(if no-struc? 
	(set! geometry 
		(list 
			(make block (center 0 0) (size sx0 sy0 infinity) (material mat1))))
	(set! geometry
		(append
			(list 
				(make block (center 0 0) (size sx0 sy0 infinity) (material mat1)))
			(geometric-objects-lattice-duplicates
				(list 
					(make cylinder (center 0 0) (radius r) (height infinity) (material mat2))
					(make cylinder (center 0 0) (radius r2) (height infinity) (material air))
					(make cylinder (center (/ dx 2) (/ dy 2)) (radius r) (height infinity) (material mat2))
					(make cylinder (center (/ dx 2) (/ dy 2)) (radius r2) (height infinity) (material air))
					(make cylinder (center (/ dx 2) (/ dy -2)) (radius r) (height infinity) (material mat2))
					(make cylinder (center (/ dx 2) (/ dy -2)) (radius r2) (height infinity) (material air))
					(make cylinder (center (/ dx -2) (/ dy 2)) (radius r) (height infinity) (material mat2))
					(make cylinder (center (/ dx -2) (/ dy 2)) (radius r2) (height infinity) (material air))
					(make cylinder (center (/ dx -2) (/ dy -2)) (radius r) (height infinity) (material mat2))
					(make cylinder (center (/ dx -2) (/ dy -2)) (radius r2) (height infinity) (material air)))
				dx dy))))

(set! geometry-lattice (make lattice (size (+ sx0 (* dpml 2)) (+ sy0 (* dpml 2)) no-size)))

; excite sources
(set! sources (list
	 (make source
	  (src (make gaussian-src (frequency fcen) (fwidth df)))
	    (component Hz)
		(center 0 (- (* 0.5 sy0) 1))
		(size sx0 0))))

(set! pml-layers (list (make pml (thickness 1.0))))
; (set! pml-layers (list (make pml (direction Y) (thickness 1.0)))) ; only on y extremes

(define-param res 50)
(set! resolution res)
(set! k-point (vector3 0 0 0))
(set! ensure-periodicity true)

; define flux regions (planes)
(define refl ; reflected flux
	(add-flux fcen df nfreq
	 (make flux-region
	  (center 0 (- (* 0.5 sy0) 2)) (size sx0 0))))
(define trans ; transmitted flux
	(add-flux fcen df nfreq
	 (make flux-region
	 (center 0 (+ (* -0.5 sy0) 1)) (size sx0 0)))) ; changed from (* sy2)

; run sources
(if (not no-struc?) (load-minus-flux "refl-flux" refl))
(run-sources+
	(stop-when-fields-decayed 50 Hz
		(vector3 0 (+ (* -0.5 sy0) 1))
		1e-2)
	(at-beginning output-epsilon)
	(at-end (output-png Hz "-Zc bluered")))

; save data, etc.
(if no-struc? (save-flux "refl-flux" refl))
(display-fluxes trans refl) ; print out the flux spectrum
