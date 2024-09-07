; this script calculates reflectance for a hexagonal lattice tilted 19.1º

(define-param no-struc? true)
(use-output-directory "output")

; input parameters

(define-param r 0.632) ; melanosome radius
(define a (* 2 r)) ; lattice constant
(define-param ra 0.25) ; air radius divided by rod radius
(define-param cor 1) ; width of cortex in um
(define-param nx 4) ; # rod layers

; source attributes				

(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 400) ; number of frequencies at which to compute flux
(define-param TE? true)

; material properties				

(define-param ker 2.43)	; permittivity of keratin
(define-param mel 4.00)	; permittivity of melanin
(define-param k1 0.01) ; extinction coeff keratin
(define-param k2 0.01) ; extincion coeff melanin
(define mel2 
	(make medium (epsilon mel) 
	  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
(define ker2 
	(make medium (epsilon ker) 
	  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))

; geometrical parameters

(define theta1 (atan (/ (sin (/ pi 3)) (+ 2 (cos (/ pi 3))))))
(define theta2 (- (/ pi 2) (/ pi 3) theta1))
(define theta3 (- (/ pi 3) theta2))
(define dx (* a (sin theta1)))
(define dy (* a (cos theta1)))
(define sx1 (* 14 dx))
(define sy (+ (* 2 dy) (* a (sin theta3))))
(define-param pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define-param sub 1) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

; make geometry

(set! geometry-lattice (make lattice (size sx sy no-size)))
(set! geometry	
	(if no-struc?	
		(list 
			(make block (center 0 0) (size sx sy infinity) (material ker2)))
		(append
	 		(list 
	 			(make block (center 0 0) (size sx sy infinity)
	 			  (material ker2)))
			(geometric-objects-duplicates (vector3 (* a (cos theta2)) (* -1 a (sin theta2))) 0 5
				(list
					(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
						(radius r) (height infinity) 
						(material mel2))
					 (make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
						(radius (* ra r)) (height infinity) 
						(material air))))
			(geometric-objects-duplicates (vector3 (* a (cos theta2)) (* -1 a (sin theta2))) 0 5
				(list
					(make cylinder (center (+ pad dpml cor r (* -0.5 sx) dx) (- (* 0.5 sy) dy)) 
						(radius r) (height infinity) 
						(material mel2))
					 (make cylinder (center (+ pad dpml cor r (* -0.5 sx) dx) (- (* 0.5 sy) dy)) 
						(radius (* ra r)) (height infinity)
						(material air))))
			(geometric-objects-duplicates (vector3 (* a (cos theta2)) (* -1 a (sin theta2))) 0 5
				(list
					(make cylinder (center (+ pad dpml cor r (* -0.5 sx) (* dx 2)) (- (* 0.5 sy) (* 2 dy))) 
						(radius r) (height infinity) 
						(material mel2))
					 (make cylinder (center (+ pad dpml cor r (* -0.5 sx) (* dx 2)) (- (* 0.5 sy) (* 2 dy))) 
						(radius (* ra r)) (height infinity) 
						(material air))))
			(list
				(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
					(radius r) (height infinity)
					(material mel2))
				(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
					(radius (* ra r)) (height infinity) 
					(material air))))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 100)
(set! resolution res)
(set! k-point (vector3 0 0 0))
(set! ensure-periodicity true)

; source characteristics

(set! sources 
	(list
	 (make source
	  (src (make gaussian-src (frequency fcen) (fwidth df)))
		(if TE?
	 		(component Hz)
	 		(component Ez))
		(center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy))))

(define refl ; reflected flux
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

(define trans ; transmitted flux
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ;changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

; run sources

(if TE?
	(run-sources+ 
		(stop-when-fields-decayed 50 Hz
	  		(vector3 (+ dpml 2 (* -0.5 sx)) 0) 1e-2)
		 	 (at-beginning output-epsilon)
		  	  (at-end (output-png Hz "-Zc bluered")))
	(run-sources+ 
		(stop-when-fields-decayed 50 Ez
	 		(vector3 (+ dpml 2 (* -0.5 sx)) 0) 1e-2)
		 	 (at-beginning output-epsilon)
	 	  	  (at-end (output-png Ez "-Zc bluered"))))

; output data, load flux, etc.

(if no-struc? (save-flux "refl-flux" refl))
(display-fluxes trans refl) 					; print out the flux spectrum