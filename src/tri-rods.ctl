; Simulation of reflection from a triangular lattice of rods (either solid or hollow)
; User can specify polarization with 'TE?' argument and incident angle with 'theta'
; Written by Chad M. Eliason, Jan 24, 2013

(define-param no-struc? true)
(define-param no-cortex? true)
(define-param TE? true)
(use-output-directory "output")

; parameters describing the geometry
(define-param r_out 0.5) ; normalized melanin radius (r/a)
(define-param r_in 0) ; radius of inner air space (radius divided by r_out)
(define-param a1 .1) ; lattice constant 1 (parallel to surface)
(define-param a2 0) ; lattice constant 2 (plane at 60º to surface), if both are the same this should be set to 0
(if (= a2 0) (define a2 a1) (define a2 a2)) ; checks if both lattice constants defined
(define-param cor 1) ; width of cortex in um
(define-param nx 4) ; number of rod layers
(define r (* (min a1 a2) r_out))
(define r2 (* r r_in)) ; inner radius of air space

; source attributes				
(define fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4) ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 400) ; number of frequencies at which to compute flux

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
        (make block (center 0 0) (size sx sy infinity) (material keratin)))
      (append
        (list 
          (make block (center 0 0) (size sx sy infinity)
            (material keratin)))
          (geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
            (list 
              (make cylinder 
                (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
                (radius r) (height infinity) 
                (material melanin))
              (make cylinder 
                (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
                (radius r2) (height infinity) 
                (material air))))
          (geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
            (list
              (make cylinder 
                (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
                (radius r) (height infinity) 
                (material melanin))
              (make cylinder 
                (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
                (radius r2) (height infinity) 
                (material air))))
          (geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx2 1)
            (if (> nx1 0)
              (list
                (make cylinder 
                  (center (+ pad dpml cor r spc (* -0.5 sx)) 0) 
                  (radius r) (height infinity) 
                  (material melanin))
                (make cylinder 
                  (center (+ pad dpml cor r spc (* -0.5 sx)) 0) 
                  (radius r2) (height infinity) 
                  (material air)))
                (make cylinder 
                  (center (+ pad dpml cor r spc (* -0.5 sx)) 0) 
                  (radius r) (height infinity) 
                  (material keratin)))))))
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
				(material keratin)))
			  (geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
				(list
				 (make cylinder 
				 	(center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
				 	(radius r) (height infinity) 
					(material melanin))
				 (make cylinder 
				 	(center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy)) 
				 	(radius r2) (height infinity) 
					(material air))))
				(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx1 1)
					(list
					 (make cylinder 
					 	(center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
					 	(radius r) (height infinity) 
						(material melanin))
					 (make cylinder 
					 	(center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy)) 
					 	(radius r2) (height infinity) 
						(material air))))
				(geometric-objects-duplicates (vector3 (* spc 2) 0) 0 (- nx2 1)
				 (if (> nx1 0)
					 (list
					  (make cylinder 
					  	(center (+ pad dpml cor r spc (* -0.5 sx)) 0) 
					  	(radius r) (height infinity) 
						(material melanin))
					  (make cylinder 
					  	(center (+ pad dpml cor r spc (* -0.5 sx)) 0) 
					  	(radius r2) (height infinity) 
						(material air))) 
					  (make cylinder 
					  	(center (+ pad dpml cor r spc (* -0.5 sx)) 0) 
					  	(radius r) (height infinity) 
						(material keratin))))))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 10)
(set! resolution res)

(define-param theta 0)	 ; inc angle in degrees
(define theta_rad (* pi (/ theta 180)))  ; angle of incidence (with respect to y-axis)
(define ky (* fcen (sin theta_rad)))  ; specify value here (2)

(set! k-point (vector3 0 ky 0))
(set! ensure-periodicity true)

;(set! k-point (vector3 0 0 0))

; define custom amp function
(define (my-amp-func p) 
   (exp (* 0+2i pi ky (vector3-y p))))

(set! sources 
	(list
	 (make source
	  (src (make gaussian-src (frequency fcen) (fwidth df)))
		(if TE?
	 		(component Hz)
	 		(component Ez))
		(center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy)
		(amp-func my-amp-func))))

(define refl ; reflected flux
	(add-flux fcen df nfreq
	 		(make flux-region
	  		  (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2))
          (direction X))))

(define trans ; transmitted flux
	(add-flux fcen df nfreq
			(make flux-region
	  		  (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))
          (direction X)))) ; changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

; uncomment this to test the script
; (run-until 10
;   (at-beginning output-epsilon)
;     (at-every 0.6 (output-png Ez "-Zc bluered")))

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
(display-fluxes trans refl)  ; print out the flux spectrum
