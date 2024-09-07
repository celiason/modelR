(define-param no-struc? true)

(use-output-directory "output")

; geometry parameters
(define-param d1 .1)   ; thickness of material 1
(define-param d2 .1)   ; thickness of material 2
(define-param cor 0)  ; width of cortex in um
(define-param nx 4)   ; # number of periods layers
(define a (+ d1 d2))  ; lattice constant

; source attributes
(define-param fcen 2.4)  ; pulse center frequency 2.4 (for 300 - 700 nm) or 5
(define-param df 4)      ; pulse width (in frequency) 2 (for 300 - 700 nm) or 8
(define-param nfreq 100) ; number of frequencies at which to compute flux

; material properties				
(define-param n1 1.56)        ; layer 1 RI
(define-param n2 2)           ; layer 2 RI
; (define-param absorbing false)           ; whether to use absorbing pigment
(define-param n3 1)           ; substrate RI
(define e1 (* n1 n1)) ; permittivity of keratin
(define e2 (* n2 n2)) ; permittivity of melanin
(define e3 (* n3 n3)) ; subtrait refractive index
(define-param k1 0)         ; layer 1 extinction coef
(define-param k2 0.1)         ; layer 2 extinction coef
(define-param k3 0)        ; substrate extinction coef

(define mat1 (make medium 
                (epsilon e1) 
                (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt e1) k1)) e1))))

(define mat2 (make medium 
                (epsilon e2) 
                (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt e2) k2)) e2))))

(define mat3 (make medium 
                (epsilon e3) 
                (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt e3) k3)) e3))))

; need option to toggle pigment absorption
(define mat4 
  (make dielectric (epsilon 2.59)
            (E-susceptibilities 
             (make lorentzian-susceptibility
               (frequency 3.78) (gamma .32) (sigma .32))
             (make lorentzian-susceptibility
               (frequency 3.51) (gamma 3) (sigma .36))
            (make lorentzian-susceptibility
               (frequency 2.23) (gamma 2.4) (sigma -.2))
             )))

; cell dimensions
(define sy a) 				      ; size of cell in y direction
(define sx1 (* a nx))       ; size of layer in x direction
(define-param pad 4) 		    ; air space above cortex
(define-param sub 2) 		    ; substrate thickness
(define dpml 1)				      ; perfectly matched layer (PML) thickness
(define sx (+ pad sub cor (* 2 dpml) sx1))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! geometry	
  (if no-struc?
    (list
      (make block (center 0 0) (size sx sy infinity) (material air)))
    (append
      (list
      (make block (center 0 0) (size sx sy infinity) (material air))
      (make block (center (+ dpml pad (/ cor 2) (* -0.5 sx)) 0)
                        (size cor sy infinity)
                        (material mat1)))
      (geometric-object-duplicates (vector3 a 0) 0 (- nx 1)
            (make block (center (+ dpml pad cor (/ d1 2) (* -0.5 sx)) 0)
                        (size d1 sy infinity)
                        (material mat1)))
        (geometric-object-duplicates (vector3 a 0) 0 (- nx 1)
            (make block (center (+ dpml pad cor d1 (/ d2 2) (* -0.5 sx)) 0)
                        (size d2 sy infinity)
                        (if absorbing
                          (material mat4)
                          (material mat2)
                          ))))))

(set! pml-layers (list (make pml (direction X) (thickness dpml)))) ; only on x extremes

(define-param res 40)
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

; reflected flux plane
(define refl
  (add-flux fcen df nfreq
    (make flux-region
      (center (+ dpml 2 (* -0.5 sx)) 0) (size 0 (* sy 2)))))

; transmitted flux plane
(define trans
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
