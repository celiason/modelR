(define-param no-struc? true)
(define-param no-cortex? false)
(use-output-directory "output")

; geometry parameters
(define-param d 1)    ; normalized diameter of melanin sheet (d/a)
(define-param ra 0.5) ; proportion of air inside mel
(define-param a 0.2)  ; spacing between melanin sheets
(define d (* d a))
(define d1 (* d ra))  ; diameter of air sheet inside melanosomes
(define-param cor 1)  ; width of cortex in um
(define-param nx 4)   ; # rod layers

; source attributes
(define fcen 2.4)  ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define df 4)      ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define nfreq 100) ; number of frequencies at which to compute flux

; material properties				
(define-param ker 1.56)        ; permittivity of keratin
(define-param mel 4)           ; permittivity of melanin
(define-param sub 1)           ; substrate refractive index
(define-param ker (* ker ker)) ; permittivity of keratin
(define-param mel (* mel mel)) ; permittivity of melanin
(define-param sub (* sub sub)) ; subtrait refractive index
(define-param k1 0)         ; extinction coeff keratin
(define-param k2 0)         ; extincion coeff melanin
; (define-param k3 0.01)       ; extincion coeff melanin

(define ker2 (make medium 
                (epsilon ker) 
                (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt ker) k1)) ker))))
(define mel2 (make medium 
                (epsilon mel) 
                (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt mel) k2)) mel))))
; (define sub (make medium 
;                 (epsilon sub) 
;                 (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt sub) k3)) sub))))

; cell dimensions

(define sy a) 				      ; size of cell in y direction
(define sx1 (+ (* a nx) d)) ; size of layer in x direction
(define-param pad 4) 		    ; air space above cortex
(define-param sub 2) 		    ; keratin substrate thickness
(define dpml 1)				      ; perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

(if (odd? nx) 
  (list (define nx1 (+ (floor (/ nx 2)) 1)) (define nx2 (floor (/ nx 2))))
  (list (define nx1 (/ nx 2)) (define nx2 (/ nx 2))))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(if no-cortex?
  (set! geometry	
    (if no-struc?
      (list 
        (make block (center 0 0) (size sx sy infinity) (material ker2)))
      (append
        (list (make block (center 0 0) (size sx sy infinity)
                          (material ker2)))
        (geometric-object-duplicates (vector3 a 0) 0 (- nx 1)
          (make block (center (+ pad dpml cor (/ d 2) (* -0.5 sx)) 0) 
                      (size d sy infinity)
                      (material mel2)))
        (geometric-object-duplicates (vector3 a 0) 0 (- nx 1)
          (make block (center (+ pad dpml cor (/ d 2) (* -0.5 sx)) 0) 
                      (size d1 sy infinity)
                      (material air))))))

  (set! geometry	
    (if no-struc?
      (list 
        (make block (center 0 0) (size sx sy infinity) (material air)))
      (append
        (list (make block (center (/ (+ pad dpml) 2) 0) 
                          (size (+ cor sx1 sub dpml) sy infinity) 
                          (material ker2)))
        (geometric-object-duplicates (vector3 (* a 2) 0) 0 (- nx1 1)
          (make block (center (+ pad dpml cor (/ d 2) (* -0.5 sx)) 0)
                      (size d sy infinity)
                      (material mel2)))
        (geometric-object-duplicates (vector3 (* a 2) 0) 0 (- nx2 1)
          (if (> nx1 0)
            (make block (center (+ pad dpml cor (/ d 2) a (* -0.5 sx)) 0)
                        (size d sy infinity)
                        (material mel2))
            (make block (center (+ pad dpml cor (/ d 2) a (* -0.5 sx)) 0)
                        (size d sy infinity)
                        (material ker2))))))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 100)
(set! resolution res)

(define-param theta 0)						; incident angle in degrees
(define theta_rad (* pi (/ theta 180)))		; angle of incidence with respect to y-axis
(define k (* 2 (sin theta_rad)))			; incident angle in radians

(set! k-point (vector3 0 k 0))
(set! ensure-periodicity true)

; define custom amp function

(define (my-amp-func p) (exp (* 0+2i pi k (vector3-x p))))

; excite sources

(set! sources 
  (list
    (make source
      (src (make gaussian-src (frequency fcen) (fwidth df)))
           (component Ez) (center (+ dpml 1 (* -0.5 sx)) 0) (size 0 sy)
           (amp-func my-amp-func))))

(define refl 							; reflected flux
  (add-flux fcen df nfreq
    (make flux-region
      (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

(define trans 							; transmitted flux
  (add-flux fcen df nfreq
    (make flux-region
      (center (- (* 0.5 sx) dpml) 0) (size 0 (* sy 2))))) ;changed from -1 on right side

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

(run-sources+ 
  (stop-when-fields-decayed 50 Ez
    (vector3 (+ dpml 2 (* -0.5 sx)) 0) 1e-2)
  (at-beginning output-epsilon)
  (at-end (output-png Ez "-Zc bluered")))

(if no-struc? (save-flux "refl-flux" refl))
(display-fluxes trans refl) ; print out the flux spectrum
