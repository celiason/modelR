(define-param no-struc? true)
(define-param no-cortex? false)
(use-output-directory "output")

; parameters describing the geometry

(define-param ra 0.5) ; air radius proportional to melanosome radius
;(define-param r 0.110) ; outer radius of melanosome
;(define-param r2 0.110) ; inner radius of air space
(define-param a1 0.200)	; spacing b/w melanosomes
(define-param a2 0.200) ; changed from (define a2 a1) 03-31-12
(define-param cor 1) ; width of cortex in um
(define-param nx 4) ; # rod layers
(define theta (- (/ pi 2) (acos (/ (/ a1 2) a2))))

(define b (/ a1 2)) ; short axis of ellipse
(define r (/ a2 2)) ; radius from center of ellipse to edge along a2 direction
(define x (* r (cos theta))) ; x distance
(define y (* r (sin theta))) ; y distance
(define a (sqrt (/ (* x x) (- 1 (/ (* y y) (* b b)))))) ; long axis of ellipse
(define A (* 2 a))
(define B (* 2 b))

; source attributes				

(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 400) ; number of frequencies at which to compute flux

; material properties				

(define-param ker 2.43) ; permittivity of keratin
(define-param mel 4.00)	; permittivity of melanin
(define-param k1 0.01) ; extinction coeff keratin
(define-param k2 0.01) ; extincion coeff melanin
(define mel2 
	(make medium (epsilon mel) 
	  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
(define ker2 
	(make medium (epsilon ker) 
	  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))

; cell dimensions

(define spc (sqrt (- (* a2 a2) (* (/ a1 2) (/ a1 2))))) ; spacing between rods (um)
(define sy a1) ; size of cell in y direction (in um)
(define sx1 (+ (* spc (- nx 1)) (* 2 r))) ; size of hexagonal layer in x direction
(define-param pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define-param sub 2) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

(if (odd? nx)
   (list 
   	  (define nx1 (+ (floor (/ nx 2)) 1)) 
	  (define nx2 (floor (/ nx 2))))
   (list 
	  (define nx1 (/ nx 2)) 
	  (define nx2 (/ nx 2))))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(if no-cortex?
  	(set! geometry 
  		(if no-struc?
	  		(list 
	  			(make block (center 0 0) (size sx sy infinity) (material ker2)))
	 		(append
				(list 
					 (make block (center 0 0) (size sx sy infinity)
						(material ker2)))
					(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
						(list 
						 (make ellipsoid 
						 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* 0.5 sy) 0) 
						 		(size A B B) 
									(material mel2))
						 (make ellipsoid 
						 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* 0.5 sy) 0) 
						 		(size (* ra A) (* ra B) (* ra B))
									(material air))))
					(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
						(list
						 (make ellipsoid 
						 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* -0.5 sy) 0) 
						 		(size A B B)
									(material mel2))
						 (make ellipsoid 
						 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* -0.5 sy) 0) 
						 		(size (* ra A) (* ra B) (* ra B))
									(material air))))
					(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx2 1)
						(if (> nx1 0)
							(list
							 (make ellipsoid 
							 	(center (+ pad dpml cor (/ spc 2) spc (* -0.5 sx)) 0 0) 
							 		(size A B B)
										(material mel2))
							 (make ellipsoid 
							 	(center (+ pad dpml cor (/ spc 2) spc (* -0.5 sx)) 0 0) 
							 		(size (* ra A) (* ra B) (* ra B)) 
										(material air)))
							 (make ellipsoid 
							 	(center (+ pad dpml cor (/ spc 2) spc (* -0.5 sx)) 0 0) 
							 		(size A B B) 
										(material ker2)))))))
	(set! geometry 
		(if no-struc?
			(list 
			 (make block (center 0 0) (size sx sy infinity) 
			 	(material air)))
			(append
			 (list 
			  (make block 
			  	(center (/ (+ pad dpml) 2) 0) 
			  	(size (+ cor sx1 sub dpml) sy infinity) 
				(material ker2)))
			  (geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
				(list
				 (make ellipsoid 
				 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* 0.5 sy)) 
				 		(size A B B) 
							(material mel2))
				 (make ellipsoid 
				 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* 0.5 sy)) 
				 		(size (* ra A) (* ra B) (* ra B))
							(material air))))
				(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
					(list
					 (make ellipsoid 
					 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* -0.5 sy)) 
					 		(size A B B) 
								(material mel2))
					 (make ellipsoid 
					 	(center (+ pad dpml cor (/ spc 2) (* -0.5 sx)) (* -0.5 sy)) 
					 	(size (* ra A) (* ra B) (* ra B)) 
						(material air))))
				(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx2 1)
				 (if (> nx1 0)
					 (list
					  (make ellipsoid 
					  	(center (+ pad dpml cor (/ spc 2) spc (* -0.5 sx)) 0) 
					  		(size A B B) 
									(material mel2))
					  (make ellipsoid 
					  	(center (+ pad dpml cor (/ spc 2) spc (* -0.5 sx)) 0) 
					  		(size (* ra A) (* ra B) (* ra B)) 
									(material air))) 
					  (make ellipsoid 
					  	(center (+ pad dpml cor (/ spc 2) spc (* -0.5 sx)) 0) 
					  		(size A B B) 
									(material ker2))))))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 100)
(set! resolution res)

;(define-param theta 0)						; inc angle in degrees
;(define theta_rad (* pi (/ theta 180)))		; angle of incidence (with respect to y-axis)
;(define ky (* 4.4 (sin theta_rad)))			; specify value here (2)

(set! k-point (vector3 0 0 0))
(set! ensure-periodicity true)

; define custom amp function

;(define (my-amp-func p) (exp (* 0+2i pi ky (vector3-y p))))
(set! sources (list
	 		   (make source
	  	 		 (src (make gaussian-src (frequency fcen) (fwidth df)))
	     		 (component Hz) 
	     		 (center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy))))

(define refl ; reflected flux
	(add-flux fcen df nfreq
	 		(make flux-region
	  		  (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

(define trans ; transmitted flux
	(add-flux fcen df nfreq
			(make flux-region
	  		  (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ; changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

(run-sources+ 
 	(stop-when-fields-decayed 50 Hz
        (vector3 (+ dpml 2 (* -0.5 sx)) 0)    
		  1e-3)
	(at-beginning output-epsilon)
	(at-end (output-png Hz "-Zc bluered")))

(if no-struc? (save-flux "refl-flux" refl))
(display-fluxes trans refl) ; print out the flux spectrum