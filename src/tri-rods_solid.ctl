(define-param no-struc? true) ; if TRUE, no structure; do TRUE first (to normalize), then FALSE
(set! eps-averaging? false)

; parameters describing the geometry

(define-param ker 2.43) ; permittivity of keratin (RI ^2)
(define-param mel 4.00) ; permittivity of melanin
(define-param diam 0.110)
(define-param a1 0.140)
(define a2 a1) ; *** a2 and a1 the same - changed this
(define-param theta 60) ; tilt angle at intersection of a1 and a2
(define-param cor 2) ; width of cortex in um
(define-param N 4) ; number of melanin layers

;calculated based on above

(define theta1 (* (/ theta 180) pi)) ; tilt angle at intersection of a1 & a2
(define shift (- (* a2 (cos theta1)) (* 0.5 a1)))

(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5

(define-param k1 0) ; extinction coeff keratin
(define-param k2 0.01) ; extincion coeff melanin
(define mel2 (make medium (epsilon mel) (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
(define ker2 (make medium (epsilon ker) (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))

; the cell dimensions
(define spc (* a2 (sin theta1))) ; spacing between rods (um)
(define sy a1) ; size of cell in y direction (in um)
(define r (/ diam 2)) ; radius of melanosomes
(define sx1 (+ (* spc (- N 1)) (* 2 r))) ; size of hexagonal layer in x direction
(define pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define sub 2) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness

(define sx (+ pad sub cor (* 2 dpml) sx1))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes
(set! k-point (vector3 0 0 0)); periodic boundary conditions
(set! ensure-periodicity true)

; geometry of the cell

(set! geometry
	
	(if no-struc?
	
	(list 
			(make block (center 0 0) (size sx sy infinity)
				(material air)))
	
	(list 
			(make block (center (/ (+ pad dpml) 2) 0) (size (+ cor sx1 sub dpml) sy infinity) 
				(material ker2))
		
			(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
				(radius r) (height infinity) (material mel2))
			(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
				(radius r) (height infinity) (material mel2))

			(make cylinder (center (+ pad dpml cor r spc (* -0.5 sx)) (+ 0 shift)) 
				(radius r) (height infinity) (material mel2))
			(make cylinder (center (+ pad dpml cor r spc (* -0.5 sx)) (+ (- 0 a1) shift)) 
				(radius r) (height infinity) (material mel2))
		
			(make cylinder (center (+ pad dpml cor r (* spc 2) (* -0.5 sx)) (+ (* 0.5 sy) (* 2 shift))) 
				(radius r) (height infinity) (material mel2))
			(make cylinder (center (+ pad dpml cor r (* spc 2) (* -0.5 sx)) (+ (* -0.5 sy) (* 2 shift))) 
				(radius r) (height infinity) (material mel2))
		
			(make cylinder (center (+ pad dpml cor r (* spc 3) (* -0.5 sx)) (+ 0 (* 3 shift))) 
				(radius r) (height infinity) (material mel2))
			(make cylinder (center (+ pad dpml cor r (* spc 3) (* -0.5 sx)) (+ (- 0 a1) (* 3 shift))) 
				(radius r) (height infinity) (material mel2))

			(make cylinder (center (+ pad dpml cor r (* spc 4) (* -0.5 sx)) (+ (* 0.5 sy) (* 4 shift))) 
				(radius r) (height infinity) (material mel2))
			(make cylinder (center (+ pad dpml cor r (* spc 4) (* -0.5 sx)) (+ (* -0.5 sy) (* 4 shift))) 
				(radius r) (height infinity) (material mel2))
)))
	;		(make cylinder (center (+ pad dpml cor r (* spc 5) (* -0.5 sx)) (+ 0 (* 5 shift))) 
	;			(radius r) (height infinity) (material mel2))
	;		(make cylinder (center (+ pad dpml cor r (* spc 5) (* -0.5 sx)) (+ (- 0 a1) (* 5 shift))) 
	;			(radius r) (height infinity) (material mel2))
;)))		
	;		(make cylinder (center (+ pad dpml cor r (* spc 6) (* -0.5 sx)) (+ (* 0.5 sy) (* 6 shift))) 
	;			(radius r) (height infinity) (material mel2))
	;		(make cylinder (center (+ pad dpml cor r (* spc 6) (* -0.5 sx)) (+ (* -0.5 sy) (* 6 shift))) 
	;			(radius r) (height infinity) (material mel2))
;)))

(set! resolution 150)

(define df 2) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 50); number of frequencies at which to compute flux

(set! sources (list
	(make source
		(src (make gaussian-src (frequency fcen) (fwidth df)))
			(component Ez)
			(center (+ dpml 1 (* -0.5 sx)) 0)
			(size 0 sy))))
					
(set! symmetries (list (make mirror-sym (direction Y))))

(define refl ; reflected flux
	(add-flux fcen df nfreq
		(make flux-region
			(center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

(define trans ; transmitted flux
	(add-flux fcen df nfreq
		(make flux-region
			(center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ;changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

(run-until 100
	(at-beginning output-epsilon)
	(at-end (output-png Ez "-Zc bluered")))

(if no-struc? (save-flux "refl-flux" refl))

(display-fluxes trans refl) ; print out the flux spectrum
