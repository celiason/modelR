(define-param no-struc? true)
(define-param TE? true)
(use-output-directory "output")

;
;
;    a2
; o <--> o  
; ^
; | a1
; v
; o      o 
;
; parameters describing the geometry
(define-param a1 0.140)
(define-param a2 0.140)
(define-param rmel 0.5)
(define-param rair 0.3)
(define rmel (* a2 rmel))
(define rair (* a2 rair))
; (define a1 a)
; (define a2 a)
(define-param cor 0) ; width of cortex in um
(define-param nx 5) ; number of layers

(define-param k0 2.4)  ; ??

; source attributes
(define-param fcen 2.4)	; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define-param df 2) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define-param nfreq 100) ; number of frequencies at which to compute flux

; material properties
(define-param ker 2.43) ; permittivity of keratin
(define-param mel 4.00)	; permittivity of melanin
(define-param k1 0.01) ; extinction coeff keratin
(define-param k2 0.01) ; extincion coeff melanin
(define mel2 (make medium (epsilon mel) 
  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
(define ker2 (make medium (epsilon ker) 
  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))

; cell dimensions
(define sy a2) ; size of cell in y direction (in um)
(define sx1 (+ (* a1 (- nx 1)) (* 2 rmel))) ; size of hexagonal layer in x direction
(define-param pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define-param sub 2) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

; (if (odd? nx) 
; 	(list (define nx1 (+ (floor (/ nx 2)) 1)) (define nx2 (floor (/ nx 2))))
; 	(list (define nx1 (/ nx 2)) (define nx2 (/ nx 2))))

; fix for new version of scheme:
(define nx1 (if (odd? nx) (+ (floor (/ nx 2)) 1) (/ nx 2)))
(define nx2 (if (odd? nx) (floor (/ nx 2)) (/ nx 2)))


(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! geometry	
	(if no-struc?
		(list 
			(make block (center 0 0) (size sx sy infinity) (material air)))
		(append
	 		(list 
	 			(make block (center (/ (+ pad dpml) 2) 0) (size (+ cor sx1 sub dpml) sy infinity)
	 			; (make block (center 0 0) (size sx sy infinity)
	 			  (material ker2)))
			(geometric-object-duplicates (vector3 a1 0) 0 (- nx1 1)
    	 		(make cylinder (center (+ pad dpml cor rmel (* -0.5 sx)) (* 0.5 sy)) 
    	 		  (radius rmel) (height infinity) 
		 		  (material mel2)))
			(geometric-object-duplicates (vector3 a1 0) 0 (- nx1 1)
				(make cylinder (center (+ pad dpml cor rmel (* -0.5 sx)) (* -0.5 sy)) 
				  (radius rmel) (height infinity) 
				  (material mel2)))
			; air spaces between melanosomes
			(geometric-object-duplicates (vector3 a1 0) 0 (- nx1 2)
				(make cylinder (center (+ pad dpml cor rair (/ a1 2) (* -0.5 sx)) 0) 
				  (radius rair) (height infinity)
				  (material air))))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 100)
(set! resolution res)

;(define-param theta 0) ; inc angle in degrees
;(define theta_rad (* pi (/ theta 180)))	; angle of incidence (with respect to y-axis)
;(define k (* 3.4 (sin theta_rad))) ; specify value here (2)

(set! k-point (vector3 0 0 0))
(set! ensure-periodicity true)

; define custom amp function
;(define (my-amp-func p) (exp (* 0+2i pi k (vector3-x p))))

(set! sources 
	(list
	 (make source
	  (src (make gaussian-src (frequency fcen) (fwidth df)))
		(if TE?
	 		(component Hz)
	 		(component Ez))
		(center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy))))

; reflected flux plane
(define refl
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

; transmitted flux plane
(define trans
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ; changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

; run sources
(if TE?
	(run-sources+ 
		(stop-when-fields-decayed 50 Hz
	  		(vector3 (+ dpml 2 (* -0.5 sx)) 0) 1e-3)
		 	 (at-beginning output-epsilon)
		  	  (at-end (output-png Hz "-Zc bluered")))
	(run-sources+ 
		(stop-when-fields-decayed 50 Ez
	 		(vector3 (+ dpml 2 (* -0.5 sx)) 0) 1e-3)
		 	 (at-beginning output-epsilon)
	 	  	  (at-end (output-png Ez "-Zc bluered"))))

(if no-struc? (save-flux "refl-flux" refl))
(display-fluxes trans refl) ; print out the flux spectrum
