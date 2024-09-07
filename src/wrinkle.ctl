; Wrinkle nanostructures simulated using sinusoidal gratins

(define-param no-struc? true)

(use-output-directory "output")

(define-param L .2) ;Period of the sinusoidal boundary
(define-param q (/ (* 2 pi) L)) ;Wave vector of the sinusoidal boundary
(define-param A .1) ;Amplitude of the " "
(define-param N 5) ;Number of single period repetions to right and left of central period

(define-param dx (* (+ (* N 2) 1) L)) ;Box basic dimensions
(define-param dpml 2)
(define-param dy (+ dx (* 2 dpml))) ;Box dimensions
(define-param frec 1)
(define-param mm 10) ;number of slices I need to create boundary
(define-param ht (/ A mm)) ;rectangles' height
; (define Dy (/ A m)) ;rectangles' height

(define-param theta (/ pi 4))

(define kx (* frec (sin theta)))
(define (my-amp-func p)(exp (* 0+2i pi kx (vector3-x p))))

(set! k-point (vector3 kx 0 0)) ;periodic boundary condition along direction without mirror (kx cons.)

(set! geometry-lattice (make lattice (size dx dy no-size)))

; (define-param block_list (list)) ;creation of list of superimposed rectangles, THIS DOESN'T WORK

(define (doit x x-max dx)
    (if (<= x x-max)
        (begin
            (set! geometry
                (append ;geometry
                    (list
                        ; (make block (center 0 0 0)
                        ;     (size dx dy no-size)
                        ;     (material air))
                            ; (material (make dielectric (epsilon 30))))
                        ; (make block (center 0 0 0)
                        ;     (size (/ dx 2) (/ dy 2) no-size)
                        ;     (material (make dielectric (epsilon 30))))
                        (geometric-object-duplicates (vector3 L 0) 0 N
                            (make block
                                (center 0 (* dy (- x (/ 1 2))))
                                (size (- L (* (/ 2 q) (asin (* x (/ ht A))))) ht)
                                (material (make dielectric (epsilon 30)))))
                        ; (geometric-object-duplicates (vector3 (* -1 L) 0) 0 N
                        ;     (make block
                        ;         (center 0 (* dy (- x (/ 1 2))))
                        ;         (size (- L (* (/ 2 q) (asin (* x (/ ht A))))) ht)
                        ;         (material (make dielectric (epsilon 30)))))
                    )
                )
            )
        (doit (+ x dx) x-max dx))))

    ; ((  i 1;
    ;   (+ i 1)
    ; ))
    ; ((< m))

;     (set! block_list
;         (append block_list
;              (make block
;                  (center 0 (* Dy (- i (/ 1 2))))
;                  (size (- L (* (/ 2 q) asin ((* i (/ Dy A))) )) Dy)
;                  (material (make dielectric (epsilon 30)))
;              )
;         )
;     )
; )

(if (not no-struc?)
    (doit 1 mm 1))

; (set! geometry             ; I don't know if this works...
;       (append
;            (list
;                (make block (center 0 (/ s 2) ) (size s (/ s 4)
; infinity)(material (make dielectric (epsilon 30))))

;                (geometric-object-duplicates (vector3 L 0) 0  N  block_list )

;                (geometric-object-duplicates (vector3 -L 0) 0  N  block_list
; )
;            )
;       )
; ) ; creation of geometric domain

(set! sources (list
             (make source
              (src (make continuous-src (frequency frec)))
              (component Ez)
              (center 0 (/ dx -2))
              (size dx 0)
              (amp-func my-amp-func))))

(set! pml-layers (list (make pml (thickness dpml) (direction Y))))

(set! resolution 40)

(run-until 100
             (at-beginning output-epsilon)
             ;(at-end output-efield-z)
             ;(at-end (output-png Ez "-Zc dkbluered"))
             )
