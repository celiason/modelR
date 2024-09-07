(set! eps-averaging? true)
(use-output-directory "output")
(define-param no-struc? true) ; if TRUE, no structure; do TRUE first (to normalize), then FALSE

; some parameters to describe the geometry
(define-param ker 2.43) 	; permittivity of keratin (RI ^2)
(define-param mel 4.00) 	; permittivity of melanin
(define-param k1 0) 		; keratin extinction coeff
(define-param k2 0.05) 	; melanin extinction coeff
(define e2 (* 2 (sqrt ker) k1)) ; epsilon2 of keratin
(define e3 (* 2 (sqrt mel) k2))	; epsilon2 of melanin (basically, extinction coeff)

(define-param ra 0.5)

(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 100); number of frequencies at which to compute flux

(define ker (make medium (epsilon ker) (D-conductivity (/ (* 2 pi fcen e2) ker))))
(define mel (make medium (epsilon mel) (D-conductivity (/ (* 2 pi fcen e3) mel))))

; the cell dimensions
(define-param a1 .680) ; lattice constant #1
(define-param a2 .556) ; lattice constant #2
(define r (/ a2 2)) ; radius of outer melanin layer
(define gamma (* 2 (acos (/ a1 (* 2 a2)))))
(define sy (* 3 a2 (sin gamma))) ; size of cell in y direction (in um)

(define r1 (* ra r))

(define sx1 (* a2 3)) 	; size of hexagonal layer in x direction
(define-param cor 0) 					; width of cortex in um
(define pad 4) 								; air space above (( FIX!!! and behind)) cortex in a
(define sub 2) 								; keratin substrate thickness (in units of a = 1um)
(define dpml 1)								; perfectly matched layer (PML) thickness

(define sx (+ pad sub cor (* 2 dpml) (* sx1 2)))

; computational cell
(set! geometry-lattice (make lattice (size sx sy no-size)))

; set absorbing boundaries
(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes
(set! k-point (vector3 0 0 0)); periodic boundary conditions
;(set! ensure-periodicity true)

; setting the geometry of the cell
(set! geometry
	
	(if no-struc?
	
	(list 
			(make block (center 0 0) (size sx sy infinity)
				(material ker)))
	
	(list 
		;	(make block (center (/ (+ pad dpml) 2) 0) (size (+ cor sx1 sub dpml) sy infinity) 
		;		(material ker))
		
			(make block (center 0 0) (size sx sy infinity)
				(material ker))
	
			(make cylinder 
				(center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r a2 (* -0.5 sx)) (* 0.5 sy)) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r a2 (* -0.5 sx)) (* -0.5 sy)) 
				(radius r) (height infinity) (material mel))


			(make cylinder 
				(center (+ pad dpml cor r (* a2 2) (* -0.5 sx)) (* 0.5 sy)) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r (* a2 2) (* -0.5 sx)) (* -0.5 sy)) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r (* a2 3) (* -0.5 sx)) (* 0.5 sy)) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r (* a2 3) (* -0.5 sx)) (* -0.5 sy)) 
				(radius r) (height infinity) (material mel))

			(make cylinder 
				(center (+ pad dpml cor r a2 (/ a2 3.0) (* -0.5 sx)) (- (* 0.5 sy) (/ sy 3))) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r a2 (/ a2 1.5) (* -0.5 sx)) (+ (* -0.5 sy) (/ sy 3))) 
				(radius r) (height infinity) (material mel))

			(make cylinder 
				(center (+ pad dpml cor r (* 2 a2) (/ a2 3.0) (* -0.5 sx)) (- (* 0.5 sy) (/ sy 3))) 
				(radius r) (height infinity) (material mel))
			(make cylinder 
				(center (+ pad dpml cor r (* 2 a2) (/ a2 1.5) (* -0.5 sx)) (+ (* -0.5 sy) (/ sy 3))) 
				(radius r) (height infinity) (material mel))
				
				
			(make cylinder 
				(center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r a2 (* -0.5 sx)) (* 0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r a2 (* -0.5 sx)) (* -0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (/ a2 3) (* -0.5 sx)) (- (* 0.5 sy) (/ sy 3))) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (/ a2 1.5) (* -0.5 sx)) (+ (* -0.5 sy) (/ sy 3))) 
				(radius r1) (height infinity) (material air))

			(make cylinder 
				(center (+ pad dpml cor r (* a2 2) (* -0.5 sx)) (* 0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (* a2 2) (* -0.5 sx)) (* -0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (* a2 3) (* -0.5 sx)) (* 0.5 sy)) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (* a2 3) (* -0.5 sx)) (* -0.5 sy)) 
				(radius r1) (height infinity) (material air))

			(make cylinder 
				(center (+ pad dpml cor r a2 (/ a2 3.0) (* -0.5 sx)) (- (* 0.5 sy) (/ sy 3))) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r a2 (/ a2 1.5) (* -0.5 sx)) (+ (* -0.5 sy) (/ sy 3))) 
				(radius r1) (height infinity) (material air))

			(make cylinder 
				(center (+ pad dpml cor r (* 2 a2) (/ a2 3.0) (* -0.5 sx)) (- (* 0.5 sy) (/ sy 3))) 
				(radius r1) (height infinity) (material air))
			(make cylinder 
				(center (+ pad dpml cor r (* 2 a2) (/ a2 1.5) (* -0.5 sx)) (+ (* -0.5 sy) (/ sy 3))) 
				(radius r1) (height infinity) (material air))


)))

(define-param res 50)
(set! resolution res)

(set! sources (list
	(make source
		(src (make gaussian-src (frequency fcen) (fwidth df)))
			(component Ez)
			(center (+ dpml 1 (* -0.5 sx)) 0)
			(size 0 sy))))
					
;(set! symmetries (list (make mirror-sym (direction Y))))

(define refl ; reflected flux
	(add-flux fcen df nfreq
		(make flux-region
			(center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))
					
(define trans ; transmitted flux
	(add-flux fcen df nfreq
		(make flux-region
			(center (- (* 0.5 sx) dpml) 1) (size 0 (* sy 2)))))

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

(run-sources+
	(stop-when-fields-decayed 50 Ez
		(vector3 0 (- (* 0.5 sx) dpml 1))
		1e-2)
	(at-beginning output-epsilon)
	(at-end (output-png Ez "-Zc bluered")))
	

(if no-struc? (save-flux "refl-flux" refl))

(display-fluxes trans refl) ; print out the flux spectrum