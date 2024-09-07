; this script calculates reflectance for a hexagonal lattice tilted 19.1º, or other angles

(define-param no-struc? true)
(use-output-directory "output")

; input parameters
(define-param r_out 0.5) ; mel radius / a1
(define-param r_in 0.5) ; air radius / mel radius
(define-param a1 .3) ; lattice constant parallel to surface
(define-param a2 0) ; other lattice constant (60º angle for hex lattice)
(if (= a2 0) (define a2 a1) (define a2 a2)) ; checks if both lattice constants defined
(define-param cor 1) ; width of cortex in um
(define-param nx 4) ; # rod layers
(define r (* r_out a1)) ; air radius divided by rod radius


; source attributes				
(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 200) ; number of frequencies at which to compute flux
(define-param TE? true)

; material properties
(define-param n_ker 1.56) ; refractive index of keratin
(define-param n_mel 2.00)	; refractive index of melanin
(define-param k_ker 0.01) ; extinction coefficient of keratin
(define-param k_mel 0.01) ; extinction coefficient of melanin
(define eps_ker (* n_ker n_ker))
(define eps_mel (* n_mel n_mel))
(define melanin
  (make medium (epsilon eps_mel) 
    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_mel) k_mel)) eps_mel))))
(define keratin
  (make medium (epsilon eps_ker) 
    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_ker) k_ker)) eps_ker))))

; geometrical parameters
(define d1 (sqrt (- (* a2 a2) (* (/ a1 2) (/ a1 2)))))  ; distance b/w bragg planes
(define-param rot 1.5) ; 2.5 for [3-1-20] plane (19.1º if hex); 1.5 for [2-1-10] plane (30º if hex)
(define d2 (* rot a1))  
(define theta1 (atan (/ d1 d2)))  ; angle b/w surface and plane of interest [311?]
(define theta3 (acos (/ (/ a1 2) a2)))
(define theta2 (- (/ pi 2) theta3 theta1))
(define d3 (sqrt (+ (* d1 d1) (* d2 d2))))
(define dx (* a1 (sin theta1)))
(define dy (* a1 (cos theta1)))
(define sx1 (* 14 dx))
(define sy (sqrt (+ (* d1 d1) (* d2 d2))))
(define-param pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define-param sub 1) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

; make geometry
(set! geometry-lattice (make lattice (size sx sy no-size)))
(set! geometry	
	(if no-struc?	
		(list 
			(make block (center 0 0) (size sx sy infinity) (material keratin)))
		(append
	 		(list 
	 			(make block (center 0 0) (size sx sy infinity)
	 			  (material keratin)))
			(geometric-objects-duplicates (vector3 (* a2 (cos theta2)) (* -1 a2 (sin theta2))) 0 5
				(list
					(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
						(radius r) (height infinity) 
						(material melanin))
					 (make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
						(radius (* r_in r)) (height infinity) 
						(material air))))
			(geometric-objects-duplicates (vector3 (* a2 (cos theta2)) (* -1 a2 (sin theta2))) 0 5
				(list
					(make cylinder (center (+ pad dpml cor r (* -0.5 sx) dx) (- (* 0.5 sy) dy)) 
						(radius r) (height infinity) 
						(material melanin))
					 (make cylinder (center (+ pad dpml cor r (* -0.5 sx) dx) (- (* 0.5 sy) dy)) 
						(radius (* r_in r)) (height infinity)
						(material air))))
			(geometric-objects-duplicates (vector3 (* a2 (cos theta2)) (* -1 a2 (sin theta2))) 0 5
				(list
					(make cylinder (center (+ pad dpml cor r (* -0.5 sx) (* dx 2)) (- (* 0.5 sy) (* 2 dy))) 
						(radius r) (height infinity) 
						(material melanin))
					 (make cylinder (center (+ pad dpml cor r (* -0.5 sx) (* dx 2)) (- (* 0.5 sy) (* 2 dy))) 
						(radius (* r_in r)) (height infinity) 
						(material air))))
			(list
				(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
					(radius r) (height infinity)
					(material melanin))
				(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
					(radius (* r_in r)) (height infinity) 
					(material air))))))

; some other cell settings
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

; flux planes
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
