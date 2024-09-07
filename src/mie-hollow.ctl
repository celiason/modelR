; MIE scattering code -- hollow rods

(define-param no-struc? true)
(define-param no-cortex? true)
(use-output-directory "output")

; parameters describing the geometry

(define-param r_out 0.110) ; outer radius of melanosome
(define-param r_in 0.110) ; inner radius of air space
(define a1 (* r_out 2)) ; size of cell in y direction

; source attributes				

(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 100) ; number of frequencies at which to compute flux
(define-param TE? true)

; material properties				

(define-param ker 2.43)	; permittivity of keratin
(define-param mel 4.00)	; permittivity of melanin
(define-param k1 0.01) ; extinction coeff keratin
(define-param k2 0.01) ; extincion coeff melanin
(define mel2 (make medium (epsilon mel) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
(define ker2 (make medium (epsilon ker) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))

; set geometry

(define-param pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad (* 2 dpml) a1)) ; size of cell in x direction
(define sy (+ (* 4 dpml) a1)) ; size of cell in y direction
(set! geometry-lattice (make lattice (size sx sy no-size)))
(set! geometry	
	(if no-struc?
		(list 
			(make block (center 0 0) (size sx sy infinity) (material ker2)))
		(list 
			(make block (center 0 0) (size sx sy infinity) (material ker2))
			(make cylinder (center (+ pad dpml r_out (* -0.5 sx)) 0)
			  (radius r_out) (height infinity) 
			   (material mel2))
			(make cylinder (center (+ pad dpml r_out (* -0.5 sx)) 0)
		 	  (radius r_in) (height infinity) 
		  	   (material air)))))

; other stuff

(set! pml-layers (list (make pml (thickness dpml))))
(set! symmetries (list (make mirror-sym (direction Y) (phase -1)))) ; for TE
; (set! symmetries (list (make mirror-sym (direction Y)))) ; for TM
(define-param res 100)
(set! resolution res)

; source characteristics

(define-param theta 0)	 ; inc angle in degrees
(define theta_rad (* pi (/ theta 180)))  ; angle of incidence (with respect to y-axis)
(define ky (* fcen (sin theta_rad)))  ; specify value here (2)
(set! k-point (vector3 0 ky 0))
(set! ensure-periodicity true)
;(set! k-point (vector3 0 0 0))

; define custom amp function

(define (my-amp-func p) 
   (exp (* 0+2i pi ky (vector3-y p))))

;  source

(set! sources 
	(list
	 (make source
	  (src (make gaussian-src (frequency fcen) (fwidth df)))
		(if TE?
	 		(component Hz)
	 		(component Ez))
		(center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy)
		(amp-func my-amp-func))))

; define flux planes

(define refl1 ; reflected flux normal incidence
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))
(define trans ; transmitted flux
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ;changed from -1 on right side
(if (not no-struc?) (load-minus-flux "refl-flux" refl1))

;(run-until 20
;	(at-beginning output-epsilon)
;	 (at-every 0.5 (output-png Ez "-Zc bluered"))
;	 (at-end (output-png Ez "-Zc bluered")))

; excite sources

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

; output

(if no-struc? (save-flux "refl-flux" refl1))
(display-fluxes trans refl1) ; print out the flux spectrum