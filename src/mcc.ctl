; User can specify polarization with 'TE?' argument and incident angle with 'theta'

(define-param no-struc? true)
; (define-param no-cortex? false)
(define-param TE? true)
(use-output-directory "output")

; geometry parameters
(define-param r_out 0.5) ; normalized melanin radius (r/a)
(define-param r_in 0)    ; radius of inner air space (radius divided by r_out)
(define-param a1 .1)     ; lattice constant parallel to surface
(define-param cor 0)     ; width of cortex (um)
(define r (* a1 r_out))  ; particle radius
(define r2 (* r r_in))   ; inner radius of air space

; light source
(define-param fcen 2.4) ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define-param df 4)     ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 200)      ; number of frequencies at which to compute flux

; material properties				
(define-param n_ker 1.56) ; refractive index of keratin
(define-param n_mel 2)    ; refractive index of melanin
(define-param n_mat 1)    ; refractive index of air (material around mels)
(define-param n_out 1)    ; refractive index of outside material (air usually)
(define-param n_in 1)     ; refractive index of inside material (air usually)
(define-param n_sub 1)    ; refractive index of the substrate beneath melanosomes

(define-param k_ker 0)    ; extinction coefficient of keratin
(define-param k_mel 0.1)  ; extinction coefficient of melanin
(define-param k_mat 0)    ; extinction coefficient of matrix
(define-param k_out 0)    ; extinction coeff outside material
(define-param k_in 0)     ; extinction coeff inside material
(define-param k_sub 0)

(define eps_ker (* n_ker n_ker))
(define eps_mel (* n_mel n_mel))
(define eps_mat (* n_mat n_mat))
(define eps_out (* n_out n_out))
(define eps_in (* n_in n_in))
(define eps_sub (* n_sub n_sub))

(define melanin (make medium 
                    (epsilon eps_mel) 
                    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_mel) k_mel)) eps_mel))))
(define keratin (make medium
                    (epsilon eps_ker) 
                    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_ker) k_ker)) eps_ker))))
(define matrix (make medium
                    (epsilon eps_mat) 
                    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_mat) k_mat)) eps_mat))))
(define outside (make medium
                    (epsilon eps_out) 
                    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_out) k_out)) eps_out))))
(define inside (make medium
                    (epsilon eps_in) 
                    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_in) k_in)) eps_in))))
(define substrate (make medium
                    (epsilon eps_sub) 
                    (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps_sub) k_sub)) eps_sub))))

; whether to have material "bleed" in-between melanosomes
(define-param bleed_top? false)
(define-param bleed_bot? false)
(if bleed_top?
  (define bleed_top r)
  (define bleed_top 0))
(if bleed_bot?
  (define bleed_bot r)
  (define bleed_bot 0))


; cell dimensions
(define sy a1)         ; size of cell in y direction (in um)
(define sx1 a1)        ; size of hexagonal layer in x direction
(define-param pad 4)   ; air space above (( FIX!!! and behind)) cortex in a
(define-param sub 2)   ; keratin substrate thickness (in units of a = 1um)
(define dpml 1)        ; perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! geometry
  (if no-struc?
      (list
        (make block (center 0 0)
                    (size sx sy infinity)
                    (material outside)))
		 (list
        (make block (center 0 0)
                    (size sx sy infinity)
                    (material outside))
        (make block (center (+ (* -0.5 sx) dpml pad cor r) 0)
                    (size (* r 2) sy)
                    (material matrix))
        (make block (center (+ dpml pad (/ (+ cor bleed_top) 2) (* -0.5 sx)) 0)
                    (size (+ cor bleed_top) sy)
			              (material keratin))
        (make block (center (- (* 0.5 sx) (/ (+ dpml sub bleed_bot) 2)) 0)
                    (size (+ sub dpml bleed_bot) sy)
                    (material substrate))
    		(make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* 0.5 sy))
    		               (radius r)
                       (height infinity)
    		               (material melanin))
        (make cylinder (center (+ pad dpml cor r (* -0.5 sx)) (* -0.5 sy))
                       (radius r)
                       (height infinity)
                       (material melanin)))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 10)
(set! resolution res)

(define-param theta 0)	 ; inc angle in degrees
(define theta_rad (* pi (/ theta 180)))  ; angle of incidence (with respect to y-axis)
(define ky (* fcen (sin theta_rad)))  ; specify value here (2)

(set! k-point (vector3 0 ky 0))
(set! ensure-periodicity true)
;(set! k-point (vector3 0 0 0))

; symmetry (not working right now)
; (if (= 0 theta)
;   (if TE? (set! symmetries (list (make mirror-sym (direction Y) (phase -1))))
;           (set! symmetries (list (make mirror-sym (direction Y) (phase 1))))))


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
;     (at-every 0.6 (output-png Hz "-Zc bluered")))

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
