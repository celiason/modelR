; this script calculates reflectance for a hexagonal lattice tilted 19.1º, or other angles

(define-param no-struc? true)
(use-output-directory "output")

; input parameters

(define-param ra1 0.5) ; mel radius / a1
(define-param ra 0.5) ; air radius / mel radius
(define-param a1 .3) ; lattice constant parallel to surface
(define-param a2 .3) ; other lattice constant (60º angle for hex lattice)
(define r (* ra1 a1)) ; air radius divided by rod radius
(define-param cor 0) ; width of cortex in um
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

(define d1 (sqrt (- (* a2 a2) (* (/ a1 2) (/ a1 2)))))  ; distance b/w bragg planes
(define d2 (* 1.5 a1))  ; 2.5 for [3-1-20] plane (19.1º if hex); 1.5 for [2-1-10] plane (30º if hex)
(define theta1 (atan (/ d1 d2)))  ; angle b/w surface and plane of interest [311?]
(define theta3 (acos (/ (/ a1 2) a2)))
(define theta2 (- (/ pi 2) theta3 theta1))
(define d3 (sqrt (+ (* d1 d1) (* d2 d2))))
(define dx (* a1 (sin theta1)))
(define dy (* a1 (cos theta1)))
(define sx1 (* 14 dx))
(define sy (sqrt (+ (* d1 d1) (* d2 d2))))
(define-param pad 3) ; air space above cortex
(define-param sub 0) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))
(define theta (- (/ pi 2) (acos (/ (/ a1 2) a2))))  ; 90º - angle b/w lattice vectors
(define b (/ a1 2)) ; short axis of ellipse
(define r1 (/ a2 2)) ; radius from center of ellipse to edge along a2 direction
(define x (* r1 (cos theta))) ; x distance
(define y (* r1 (sin theta))) ; y distance
(define a (sqrt (/ (* x x) (- 1 (/ (* y y) (* b b)))))) ; long axis of ellipse
(define A (* 2 a))
(define B (* 2 b))

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
			(geometric-objects-duplicates (vector3 (* a2 (cos theta2)) (* -1 a2 (sin theta2))) 0 5
				(list
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
						(material mel2) (size A B B)
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (sin theta1) (cos theta1) 0)
							(e3 0 0 1))
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
						(material air) (size (* ra A) (* ra B) (* ra B))
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))))
			(geometric-objects-duplicates (vector3 (* a2 (cos theta2)) (* -1 a2 (sin theta2))) 0 5
				(list
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx) dx) (- (* 0.5 sy) dy)) 
						(material mel2) (size A B B)
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx) dx) (- (* 0.5 sy) dy)) 
						(material air) (size (* ra A) (* ra B) (* ra B))
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))))
			(geometric-objects-duplicates (vector3 (* a2 (cos theta2)) (* -1 a2 (sin theta2))) 0 5
				(list
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx) (* 2 dx)) (- (* 0.5 sy) (* 2 dy))) 
						(material mel2) (size A B B)
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx) (* 2 dx)) (- (* 0.5 sy) (* 2 dy))) 
						(material air) (size (* ra A) (* ra B) (* ra B))
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))))
			(list
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
						(material mel2) (size A B B)
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))
					(make ellipsoid (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
						(material air) (size (* ra A) (* ra B) (* ra B))
							(e1 (cos theta1) (sin theta1) 0)
							(e2 (* -1 (sin theta1)) (cos theta1) 0)
							(e3 0 0 1))))))

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