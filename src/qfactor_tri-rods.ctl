; line defect slab 
; website: http://npl.kongju.ac.kr/wiki/images/3/35/Line-defect-slab.txt)

(define-param h 0.571) ; the thickness of the slab
(define-param eps 12.25) ; the dielectric constant of the slab
(define-param loweps 1.0) ; the dielectric constant of the clad
(define-param r 0.312) ; the radius of the holes
(define-param supercell-x 7) ; the (odd) number of horizontal supercell periods
(define-param supercell-y 7) ; the (odd) number of vertical supercell periods
(define-param supercell-h 4) ; height of the supercell

; triangular lattice with supercell:
(set! geometry-lattice (make lattice (size supercell-x supercell-y no-size)
      (basis1 (/ (sqrt 3) 2) 0.5)
      (basis2 (/ (sqrt 3) 2) -0.5)))

(set! geometry
      (list
        (make cylinder (center 0 0) (radius r) (height infinity) 
          (material (make dielectric (epsilon eps))))))
        ; (make cylinder (material air)
        ;   (center 0) (radius r) (height supercell-h))))

(set! geometry
  (append
 ; duplicate the bulk crystal rods over the supercell:
    (geometric-objects-lattice-duplicates geometry 1 1 4)
 ; add a rod of slab, to erase a row of air rods and form a waveguide:
 (list
  (make cylinder (center 0) (radius r) (height h)
    (material (make dielectric (epsilon eps)))))))

; (define Gamma (vector3 0 0 0 ))
; (define K' (lattice->reciprocal (vector3 0.5 0 0))) ; edge of Brillouin zone.
; (set! k-points (interpolate 4 (list Gamma K')))

(set-param! resolution 30)

(define-param fcen .15) ; pulse center frequency
(define-param df .1) ; pulse width (in frequency)
(set! sources (list
        (make source
          (src (make gaussian-src (frequency fcen) (fwidth df)))
          (component Ez) (center 0 0))))

(run-sources+ 300
    (at-beginning output-epsilon)
    (after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))
