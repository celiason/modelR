(define-param no-struc? true) ; if TRUE, no structure; do TRUE first (to normalize), then FALSE

; parameters describing geometry

(define ker 2.43) 			; permittivity of keratin (RI ^2)
(define mel 4.00) 			; permittivity of melanin
(define k 0.01) 			; melanin extinction coeff
(define k2 0.01) 			; keratin extinction coeff
(define e2 (* 2 (sqrt mel) k)) ; epsilon2 of melanin (basically, extinction coeff)
(define e3 (* 2 (sqrt ker) k2)) ; epsilon2 of keratin
(define-param diam 0.225)		; radius of air space in um
(define r1 (/ diam 2))

; the cell dimensions

(define-param a .398)				; lattice constant ***
(define r (/ a 2))					; radius of outer melanin layer
(define spc (* a (/ (sqrt 3) 2)))	; spacing between rod layers
(define sy a)						; size of cell in y direction (in um)

(define-param N 1)	; number of layers

(define sx1 (+ (* spc (- N 1)) (* 2 r))) ; size of hexagonal layer in x direction
(define-param cor 0) ; width of cortex in um
(define pad 6) ; air space above (( FIX!!! and behind)) cortex in a
(define sub 2) ; keratin substrate thickness (in units of a = 1um)
(define dpml 2); perfectly matched layer (PML) thickness

(define sx (+ pad sub cor (* 2 dpml) sx1))

; computational cell

(set! geometry-lattice (make lattice (size sx sy no-size)))

; set absorbing boundaries

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes
(set! k-point (vector3 0 0 0)); periodic boundary conditions
(set! ensure-periodicity true)

(define fcen 2.4) 	; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 2) 		; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 50)	; number of frequencies at which to compute flux

; setting the geometry of the cell

(set! geometry	
  (if no-struc?
	(list 
		(make block (center 0 0) (size sx sy infinity)
		  (material air)))
	(list 
		(make block (center (/ (+ pad dpml) 2) 0) 
			(size (+ cor sx1 sub dpml) sy infinity) 
			 (material (make medium (epsilon ker) 
					   (D-conductivity (/ (* 2 pi fcen e3) ker)))))
		(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
			(radius r) (height infinity)
			 (material (make medium (epsilon mel) 
					   (D-conductivity (/ (* 2 pi fcen e2) mel)))))
		(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
			(radius r) (height infinity)
			 (material (make medium (epsilon mel) 
					   (D-conductivity (/ (* 2 pi fcen e2) mel)))))
)))

(set! resolution 100)

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
			(center (- (* 0.5 sx) dpml) 1) (size 0 (* sy 2)))))

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

(run-until 100
	(at-beginning output-epsilon)
	(at-end (output-png Ez "-Zc bluered")))

(if no-struc? (save-flux "refl-flux" refl))

(display-fluxes trans refl) ; print out the flux spectrum