; This code calculates the extinction, scattering and absorption cross-sections of 2d cylinder.
; Code written by Bala krishna Juluri, http://juluribk.com, Feb 21, 2012
; Code modified by Chad Eliason on Mar 27, 2012

(define-param fcen 2.4) ; pulse center frequency
(define-param df 4) ; turn-on bandwidth
(define-param resl 40)
(define-param rad 0.5); radius of particle
(define-param ra 0) ; proportional amount of air within melanin rod
(define-param dpml 1) ; thickness of PML layers
(define-param sx (* 10 rad) ) ; the size of the comp cell in X, not including PML
(define sx0 (+ sx (* 2 dpml))) ; cell size in Y direction, including PML
(define sy sx); Aspect ratio of simulation domain
(define sy0 (+ sy (* 2 dpml))) ; cell size in Y direction, including PML
(define-param sample 10)
(define-param nfreq 20)
(define-param structure? false)
(define-param scattering? true)
(define-param substrate? true)
(define-param TE? true)
(set! eps-averaging? false)
(set! force-complex-fields? false)
; what do these params do? i turned them off
; (set! Courant 0.25)
(set! output-single-precision? true)

(define-param low 1)
(define-param high 4)
(define-param sub 2.43)
(define-param k_low 0)
(define-param k_high 0.1)
(define-param k_sub 0)

(define high (make medium (epsilon high) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt high) k_high)) high))))
(define background (make medium (epsilon low) 
	(D-conductivity (/ (* 2 pi fcen (* 2 (sqrt low) k_low)) low))))
(define sub (make medium (epsilon sub) 
  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt sub) k_sub)) sub))))


; Define size of lattice
(set! geometry-lattice (make lattice (size sx0 sy0 no-size)))
(set! geometry  
   (if structure?	
	   (if substrate?
          (list
             (make block (center 0 0) (size sx sy infinity) 
                  (material background))
              (make block (center 0 (* rad 3)) (size sx (* 4 rad) infinity) 
                  (material sub))
              (make cylinder (center 0 0) (radius rad) (height infinity) 
                	(material high))
            	(make cylinder (center 0 0) (radius (* ra rad)) (height infinity) 
                	(material air)))
          (list
             (make block (center 0 0) (size sx sy infinity) 
                  (material background))
              (make cylinder (center 0 0) (radius rad) (height infinity) 
                  (material high))
              (make cylinder (center 0 0) (radius (* ra rad)) (height infinity) 
                  (material air))))
          (list
                (make block (center 0 0) (size sx sy infinity)
             	    (material background)))))

; (set! symmetries (list (make mirror-sym (direction X) (phase -1)))) ; Our structure and source has mirror symmetry with odd phase. This will reduce the simulation time by half

;PML layers in X and Y-direction
(set! pml-layers (list (make pml (thickness dpml))))


;Source definition
(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (if TE? 
                  (component Hz)
                  (component Ez))
               (center 0 (* rad 4) 0)             
               		(size sx0 0))))

(set! resolution resl)

; Define a monitor half the size
(define monitor (volume (center 0 0 0) (size sx sy 0)))


(if (not structure?)
(define incident ; transmitted flux                                                
	(add-flux fcen df nfreq
		(make flux-region
		(center 0 (* rad -4)  0) (size (* rad 4) 0)))))

(define left                                       
	(add-flux fcen df nfreq
        (make flux-region
        (center (* rad -2)  0 0) (size 0 (* rad 4) ))))
(define right                                                
	(add-flux fcen df nfreq
        (make flux-region
        (center  (* rad 2)  0 0) (size 0 (* rad 4) ))))
(define top                                                
	(add-flux fcen df nfreq
        (make flux-region
        (center 0 (* rad 2) 0) (size  (* rad 4)  0))))
(define bottom                                                
	(add-flux fcen df nfreq
        (make flux-region
        (center 0 (* rad -2) 0) (size (* rad 4) 0))))


(if TE?
  (if scattering? 
  	(list 
  		(if structure? 
  			(list 
  			  (load-minus-flux "left_flux" left) 
  			  (load-minus-flux "right_flux" right) 
  			  (load-minus-flux "top_flux" top) 
  			  (load-minus-flux "bottom_flux" bottom)))
  			(run-sources+	(stop-when-fields-decayed 1 Hz (vector3 0 (* rad -4) 0) 1e-6))
  		(if (not structure?) 
  			(list 
  			  (save-flux "left_flux" left)
  			  (save-flux "right_flux" right)
  			  (save-flux "top_flux" top)
  			  (save-flux "bottom_flux" bottom)))     
  		(if (not structure?) (display-fluxes incident))
  		(if structure? (display-fluxes left right top bottom)))

  	(list 
  		(run-sources+	(stop-when-fields-decayed 1 Hz (vector3 0 (* rad -4) 0) 1e-6)
    		(in-volume monitor (at-beginning output-epsilon))
  		  (in-volume monitor (to-appended "Hz" (at-every (/ 1 (+ fcen (* 0.5 df)) sample) output-hfield-z))))
  		(if structure? (display-fluxes left right top bottom))))
  (if scattering? 
    (list 
      (if structure? 
        (list 
          (load-minus-flux "left_flux" left) 
          (load-minus-flux "right_flux" right) 
          (load-minus-flux "top_flux" top) 
          (load-minus-flux "bottom_flux" bottom)))
        (run-sources+ (stop-when-fields-decayed 1 Ez (vector3 0 (* rad -4) 0) 1e-6))
      (if (not structure?) 
        (list 
          (save-flux "left_flux" left)
          (save-flux "right_flux" right)
          (save-flux "top_flux" top)
          (save-flux "bottom_flux" bottom)))     
      (if (not structure?) (display-fluxes incident))
      (if structure? (display-fluxes left right top bottom)))

    (list 
      (run-sources+ (stop-when-fields-decayed 1 Ez (vector3 0 (* rad -4) 0) 1e-6)
        (in-volume monitor (at-beginning output-epsilon))
        (in-volume monitor (to-appended "Ez" (at-every (/ 1 (+ fcen (* 0.5 df)) sample) output-efield-z))))
      (if structure? (display-fluxes left right top bottom))))
)

; (run-until 500 
;   (at-beginning output-epsilon)
;     (to-appended "hz" 
;       (in-volume (volume (center 0 0) (size 1 1))
;       output-hfield-z)))

