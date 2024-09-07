; ############   run 1    #################

; morphological parameters
(define-param ra .5)    ; radius of melanosomes
(define-param a .200)   ; spacing between melanosomes
(define-param res 5)
(define r (* a ra))
(define pad 2.5)
(define dpml 1)
(define mag (/ a .06273424)) ; expands points so spacing between = a, value is mean nnd of simulated points

; source properties
(define fcen 3.4)
(define df 4)
(define nfreq 100)

; geometry
(define sx mag) ; width of original image
(define sy0 mag) ; height of original image
(define sy (+ pad sy0 (* 2 dpml)))
(define sz mag)

(set! geometry-lattice (make lattice (size sx sy sz)))

; material properties
(define-param n1 1.56) ; refractive index of matrix
(define-param n2 1) ; RI of rods
(define eps1 (* n1 n1)) ; permittivity of material 1 (low index)
(define eps2 (* n2 n2)) ; permittivity of material 2 (high index)
(define-param k1 0.03) ; extinction coeff mat 1
(define-param k2 0)  ; extinction coeff mat 2
(define mat1 (make medium (epsilon eps1) (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps1) k1)) eps1))))
(define mat2 (make medium (epsilon eps2) (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps2) k2)) eps2))))
(set! default-material mat1)

; define geometry
(set! geometry (list 
  (make block (center 0 0 0) (size sx sy sz) (material mat1))))

; boundary conditions, etc.
(set! pml-layers (list (make pml (direction Y) (thickness dpml))))
(set! ensure-periodicity true)
(set! k-point (vector3 0 0 0))
(set! resolution res)

; excite sources
(set! sources (list
   (make source
    (src (make gaussian-src (frequency fcen) (fwidth df)))
      (component Hz)
    (center 0 (- (* 0.5 sy) dpml 1) 0)
    (size sx 0 sz))))

; reflected flux plane
(define refl
  (add-flux fcen df nfreq
   (make flux-region
    (center 0 (- (* 0.5 sy) dpml 2) 0) (size (* sx 2) 0 (* 2 sz)))))

; transmitted flux plane
(define trans
  (add-flux fcen df nfreq
   (make flux-region
   (center 0 (+ (* -0.5 sy) dpml 1) 0) (size (* sx 2) 0 (* 2 sz)))))

(run-sources+
  (stop-when-fields-decayed 50 Hz
    (vector3 0 (- (* 0.5 sy) dpml 1) 0) 1e-2))

(save-flux "refl-flux" refl)

(display-fluxes trans refl) ; print out the flux spectrum


(reset-meep)
; (restart-fields)

; ########  run 2  #####################

(set! geometry-lattice (make lattice (size sx sy sz)))

; material properties
(define-param n1 1.56) ; refractive index of matrix
(define-param n2 1) ; RI of rods
(define eps1 (* n1 n1)) ; permittivity of material 1 (low index)
(define eps2 (* n2 n2)) ; permittivity of material 2 (high index)
(define-param k1 0.03) ; extinction coeff mat 1
(define-param k2 0)  ; extinction coeff mat 2
(define mat1 (make medium (epsilon eps1) (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps1) k1)) eps1))))
(define mat2 (make medium (epsilon eps2) (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps2) k2)) eps2))))
(set! default-material mat1)

; simulated xyz points for quasi-ordered spheres
(define xx #(0.5904994 0.5018023 0.1021994 0.2843065 0.1329397 0.1691881 0.1441221 0.3904615 0.112207 0.5538542 0.2484779 0.5155655 0.03547686 0.1623268 0.1457533 0.3077682 0.4791726 0.1982473 0.995703 0.7915184 0.921145 0.7570738 0.6516778 0.09590671 0.2986273 0.9425245 0.1663836 0.5327356 0.01847456 0.7861075 0.3682979 0.696567 0.2305814 0.8082924 0.3568684 0.9765326 0.3945234 0.2876561 0.5375429 0.8411503 0.6614553 0.4429531 0.9526596 0.6354951 0.1773679 0.3399255 0.7248083 0.6525165 0.6479251 0.8660474 0.9335623 0.4997934 0.4791287 0.8111467 0.3680828 0.3430075 0.8804687 0.3248468 0.6856731 0.5175809 0.8436074 0.9676664 0.8767241 0.03400352 0.07400258 0.003808549 0.59748 0.4614059 0.8322891 0.8214395 0.2526761 0.7301279 0.3346351 0.6128408 0.1561519 0.3170207 0.7756495 0.324539 0.4618381 0.6204036 0.6943196 0.6708936 0.3236995 0.4473595 0.8125735 0.804995 0.9818182 0.01702803 0.5024566 0.954905 0.8299605 0.674311 0.6278789 0.6329678 0.1097851 0.1546071 0.02723719 0.4444337 0.4834717 0.9757844))
(define yy #(0.5653934 0.05328814 0.9915978 0.846291 0.4320679 0.1290471 0.8053674 0.3538076 0.5886792 0.8654629 0.1387124 0.1647353 0.4279343 0.6069222 0.5351509 0.4342975 0.9607038 0.2378415 0.6141566 0.7832297 0.1707797 0.1888213 0.9827568 0.2252998 0.7166738 0.315893 0.7939768 0.2609946 0.1146217 0.2783078 0.1293871 0.7916664 0.3111401 0.382813 0.7275132 0.7706389 0.9336634 0.3431245 0.7567722 0.08037528 0.5784394 0.251254 0.6871551 0.3670848 0.01201259 0.1437349 0.3793814 0.9632638 0.4672984 0.1945135 0.605057 0.5493827 0.4510478 0.9014183 0.2481268 0.6405913 0.2959529 0.5498113 0.6751625 0.6591553 0.7091058 0.8108417 0.959631 0.1954787 0.3454989 0.9835646 0.1594574 0.8593634 0.9913529 0.4945412 0.4698695 0.05579262 0.9230615 0.2663257 0.9080967 0.9216562 0.6924878 0.5210643 0.340179 0.3710916 0.1724729 0.08041821 0.04754693 0.7446779 0.8662109 0.4961864 0.9034526 0.09592327 0.0623995 0.4075688 0.5949423 0.5768669 0.7727188 0.8700452 0.3284481 0.663932 0.8500991 0.6433635 0.4588218 0.5254537))
(define zz #(0.1744701 0.3467322 0.2410735 0.0782678 0.1284397 0.8144199 0.6767487 0.07981191 0.7835597 0.6864213 0.3938877 0.0283344 0.6340173 0.557042 0.3364913 0.5014574 0.1338081 0.0181354 0.1774071 0.4504937 0.8389027 0.002755516 0.9793372 0.2290861 0.8403748 0.0629055 0.4186815 0.2393412 0.419709 0.6787979 0.6374625 0.854438 0.6828915 0.8852846 0.6154606 0.8341184 0.5147534 0.2881183 0.4666608 0.5695504 0.9521033 0.453557 0.6179279 0.5151588 0.5755757 0.179858 0.3002358 0.4994857 0.7246803 0.246404 0.4012719 0.5647221 0.8837178 0.1411947 0.8595633 0.3998711 0.4665525 0.1863893 0.6272885 0.787357 0.01465932 0.2985665 0.9258108 0.6371388 0.8512542 0.7286154 0.6127324 0.9253099 0.3546214 0.572373 0.9327952 0.773045 0.739208 0.8364324 0.8878433 0.2958591 0.2369988 0.7173142 0.6675137 0.05015419 0.4017972 0.1895243 0.9686999 0.2200138 0.6666242 0.1159128 0.5113908 0.051538 0.8197521 0.2762749 0.7926495 0.4121722 0.07616179 0.2872248 0.4370954 0.01507575 0.0725321 0.009980036 0.3512282 0.9608616))
(define xx (vector-map (lambda (n) (- n .5)) xx))
(define yy (vector-map (lambda (n) (- n .5)) yy))
(define zz (vector-map (lambda (n) (- n .5)) zz))

; loop to place rods at specified locations
(define (doit x x-max dx)
  (if (<= x x-max)
    (begin
      (set! geometry (append geometry
        (list 
          (make sphere
            (center (* (vector-ref xx x) mag) 
                    (- (* (vector-ref yy x) mag) (/ pad 2))
                    (* (vector-ref zz x) mag))
            (radius r)
            (material mat2)))))
    (doit (+ x dx) x-max dx))))

; define geometry
(doit 0 99 1) ; do loop from a to b in steps of dx

; boundary conditions, etc.
(set! pml-layers (list (make pml (direction Y) (thickness dpml))))
(set! ensure-periodicity true)
(set! k-point (vector3 0 0 0))
(set! resolution res)

; excite sources
(set! sources (list
   (make source
    (src (make gaussian-src (frequency fcen) (fwidth df)))
      (component Hz)
    (center 0 (- (* 0.5 sy) dpml 1) 0)
    (size sx 0 sz))))

; reflected flux plane
(define refl
  (add-flux fcen df nfreq
   (make flux-region
    (center 0 (- (* 0.5 sy) dpml 2) 0) (size (* sx 2) 0 (* 2 sz)))))

; transmitted flux plane
(define trans
  (add-flux fcen df nfreq
   (make flux-region
   (center 0 (+ (* -0.5 sy) dpml 1) 0) (size (* sx 2) 0 (* 2 sz)))))

(load-minus-flux "refl-flux" refl)

(run-sources+
  (stop-when-fields-decayed 50 Hz
    (vector3 0 (- (* 0.5 sy) dpml 1) 0) 1e-2))

(display-fluxes trans refl) ; print out the flux spectrum

(new-line)
