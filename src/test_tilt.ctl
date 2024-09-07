; not sure what this script does

(define-param no-struc? true)
(use-output-directory "output")

; parameters describing the geometry

(define-param r 0.110)
(define-param ra 0.5)
(define a (* r 2))
(define d (/ (* a (sqrt 3)) 2))

; source attributes				

(define fcen 2.4)			; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 2)				; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 400) 			; number of frequencies at which to compute flux

; material properties				

(define-param ker 2.43)		; permittivity of keratin
(define-param mel 4.00)		; permittivity of melanin
(define-param k1 0) 		; extinction coeff keratin
(define-param k2 0.01) 		; extincion coeff melanin
(define mel2 
	(make medium (epsilon mel) 
	  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
(define ker2 
	(make medium (epsilon ker) 
	  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))

; cell dimensions

(define sy1 (* 6 d)) ; size of cell in y direction (in um)
(define sx1 (* 6 a)) ; size of hexagonal layer in x direction
(define-param pad 4) ; air space above (( FIX!!! and behind)) cortex in a
(define-param sub 2) ; keratin substrate thickness (in units of a = 1um)
(define dpml 1); perfectly matched layer (PML) thickness
(define sx (+ pad sub (* 2 dpml) sx1))
(define sy (+ (* 2 dpml) sy1))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! geometry	

(if no-struc?	

(list 
		(make block (center 0 0) (size sx sy infinity) (material ker2)))

(append
	(list 
		(make block (center 0 0) (size sx sy infinity)
	  	(material ker2))

(make cylinder (center (+ pad dpml (* -0.5 sx)) (* 2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* -0.5 sx)) (* 2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml a (* -0.5 sx)) (* 2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml a (* -0.5 sx)) (* 2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* 2 a) (* -0.5 sx)) (* 2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* 2 a) (* -0.5 sx)) (* 2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* 3 a) (* -0.5 sx)) (* 2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* 3 a) (* -0.5 sx)) (* 2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
	
(make cylinder (center (+ pad dpml (/ a 2) (* -0.5 sx)) d)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (/ a 2) (* -0.5 sx)) d)
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (+ a (/ a 2)) (* -0.5 sx)) d)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (+ a (/ a 2)) (* -0.5 sx)) d)
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (+ (* 2 a) (/ a 2)) (* -0.5 sx)) d)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (+ (* 2 a) (/ a 2)) (* -0.5 sx)) d)
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (+ (* 3 a) (/ a 2)) (* -0.5 sx)) d)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (+ (* 3 a) (/ a 2)) (* -0.5 sx)) d)
	(radius (* r ra)) (height infinity) 
	(material air))

(make cylinder (center (+ pad dpml a (* -0.5 sx)) 0)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml a (* -0.5 sx)) 0)
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* 2 a) (* -0.5 sx)) 0)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* 2 a) (* -0.5 sx)) 0)
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* 3 a) (* -0.5 sx)) 0)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* 3 a) (* -0.5 sx)) 0)
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* 4 a) (* -0.5 sx)) 0)
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* 4 a) (* -0.5 sx)) 0)
	(radius (* r ra)) (height infinity) 
	(material air))

(make cylinder (center (+ pad dpml (/ (* 3 a) 2) (* -0.5 sx)) (* -1 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (/ (* 3 a) 2) (* -0.5 sx)) (* -1 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (+ a (/ (* 3 a) 2)) (* -0.5 sx)) (* -1 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (+ a (/ (* 3 a) 2)) (* -0.5 sx)) (* -1 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (+ (* 2 a) (/ (* 3 a) 2)) (* -0.5 sx)) (* -1 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (+ (* 2 a) (/ (* 3 a) 2)) (* -0.5 sx)) (* -1 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (+ (* 3 a) (/ (* 3 a) 2)) (* -0.5 sx)) (* -1 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (+ (* 3 a) (/ (* 3 a) 2)) (* -0.5 sx)) (* -1 d))
	(radius (* r ra)) (height infinity) 
	(material air))
	
(make cylinder (center (+ pad dpml (* a 2) (* -0.5 sx)) (* -2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* a 2) (* -0.5 sx)) (* -2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* a 3) (* -0.5 sx)) (* -2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* a 3) (* -0.5 sx)) (* -2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* a 4) (* -0.5 sx)) (* -2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* a 4) (* -0.5 sx)) (* -2 d))
	(radius (* r ra)) (height infinity) 
	(material air))
(make cylinder (center (+ pad dpml (* a 5) (* -0.5 sx)) (* -2 d))
	(radius r) (height infinity) 
	(material mel2))
(make cylinder (center (+ pad dpml (* a 5) (* -0.5 sx)) (* -2 d))
	(radius (* r ra)) (height infinity) 
	(material air))

))))

(set! pml-layers (list (make pml (thickness dpml))))

(define-param res 100)
(set! resolution res)

;(set! k-point (vector3 0 0 0))
;(set! ensure-periodicity true)


; define source
(set! sources 
	(list
	 (make source
	  (src (make gaussian-src (frequency fcen) (fwidth df)))
	   (component Hz) (center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy)
			)))

(define refl ; reflected flux
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

(define trans ; transmitted flux
	(add-flux fcen df nfreq
	 (make flux-region
	  (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ;changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

(run-sources+ 
 (stop-when-fields-decayed 50 Hz
        (vector3 (+ dpml 2 (* -0.5 sx)) 0)    
		  1e-2)
		  (at-beginning output-epsilon)
		  (at-end (output-png Hz "-Zc bluered")))

(if no-struc? (save-flux "refl-flux" refl))
(display-fluxes trans refl) 					; print out the flux spectrum
