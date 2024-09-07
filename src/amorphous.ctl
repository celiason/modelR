(define-param no-struc? true)
(use-output-directory "output")
(set! eps-averaging? true)
(define-param periodic? true)
(define-param TE? true)

; morphological parameters
(define-param a .2)    ; space between melanosome centers
(define-param ra .5)   ; normalized mel size
(define-param r2 0)    ; air core diameter (relative to particle diameter)
(define-param pad 4)   ; space b/w rods and boundaries
(define-param sub 2)   ; substrate thickness
(define-param dpml 1)  ; width of absorbing boundary layers
(define-param d0 .082) ; unscaled particle diameter (from molec. dynam. simulations)
(define-param cor 0)   ; cortex thickness
(define mag (/ a d0))  ; expands points so spacing between = a
(define r (* ra a))    ; new rod radius
(define r2 (* r r2))   ; new air core radius
(define-param rot 0)   ; angle to rotate point set (in degrees)
(define rot (* (/ rot 180) pi))

; whether to have cortex "bleed" in-between melanosomes or not
(define-param bleed? false)
(if bleed? (define bleed r) (define bleed 0))

; account for small reflectance plane (matching experimental setup)
(define-param Omega 5) ; solid angle in degrees
(define Omega (* (/ Omega 180) pi))
(define workdist (- pad 2)) ; distance between reflectance plane and sample
(define ap_width (* workdist (tan (/ Omega 2)))) ; width of reflectance plane aperture

; default point positions
(define-param xx #(-0.412853359525165 -0.366995486054022 0.411049243181321 -0.244200141906658 0.433143294560511 0.261309745964321 -0.292650400290792 -0.212819671326066 -0.199319013107479 0.127744422197873 -0.131199221328975 0.075294982974932 0.323802926026089 -0.435403974366241 0.0573820422332009 0.416930492399945 0.330126874880211 0.418705248845652 -0.0942878461518095 -0.411830964538009 -0.21903018469239 -0.0386941422319938 0.28694190887313 -0.394909979332182 0.324709966493964 -0.370932839788204 -0.136976453147808 -0.312146085707475 -0.388372484651402 -0.00314944332657113 0.235445889918162 0.401424485990284 0.168496096462599 -0.0542737120878133 0.491681279099692 0.328600873142637 -0.0209617328168623 -0.305980552570411 -0.386283942984451 -0.153981278638005 0.241911824201397 -0.311007760730917 0.161630519585191 0.255899235257119 -0.200472065303031 -0.49412212204123 -0.150484655164677 -0.284620659488642 0.276006679456424 0.229857537852791 0.432932811699216 0.49412212204123 -0.4535386642612 0.155801057592991 0.138027454550668 -0.293719699079692 0.190964949315911 0.157812228537668 -0.194454206029756 -0.0569627844723583 0.410104805751487 -0.136997887596278 -0.0382522372296367 0.0170333050249484 -0.0725478537459235 0.247747550957811 0.0110315474847043 0.0748762318470312 -0.328270362781992 -0.209010881804389 0.0208860216612536 0.122576144190031 0.342241806710174 -0.118279724401205 0.166438336386827 -0.347026988794311 0.0313321717369533 -0.268273585669235 -0.103754548194182 0.445486821687189 -0.299210550148633 0.363991949282897 -0.492271171221776 0.339338258612909 0.251390276768566 0.0770376334829062 -0.142304732326042 -0.470295725669663 0.208524681237129 0.434716932234442 -0.487032286292815 0.126436791212908 -0.0110092110110543 0.0649900380983204 -0.395648009923454 0.0051521277990822 0.323463747860705 -0.459520867584729 0.0598584580388848 -0.102906870948684))
(define-param yy #(-0.15419750654617 0.0791308965119263 -0.173859149605544 -0.480378962403349 0.0421078183793341 -0.351905930734978 0.119369922172458 0.042528754108364 -0.270191739575757 0.493892672297404 0.385465337487702 -0.401425439712756 -0.254025116262751 -0.441257198241478 -0.0239637444320669 -0.358068414921286 -0.410135319668326 -0.272847088169908 -0.255011879236452 -0.253118355111625 0.193673933272057 -0.00260718622332956 0.144696018551216 0.00201467129244659 0.312964512960894 -0.34724914208533 0.193622716251347 -0.199102314327231 0.16648572624477 -0.093605848738278 0.448388787946598 0.160994565655783 -0.316562640507243 -0.175885883254212 0.172359723657781 -0.493892672297404 0.322712305065729 0.328962233743724 0.263283484862319 -0.419845922927862 0.356641555717988 0.212630605567424 0.395936406411801 -0.449399534808367 -0.184247872691742 -0.13017190659486 0.47356717626398 0.432200559045722 0.0612790082313586 0.204682050877211 -0.0864828395019221 0.488434140808453 0.358063300450721 -0.427095851013611 -0.23127829763531 -0.388472011706626 0.111663506937917 -0.0942421223815474 0.330575697092201 0.111609006187661 0.484455422615139 -0.121797048718388 0.241641238554152 0.440418018257563 -0.447629173809161 -0.218729810290379 -0.472019065278 -0.175941688047672 -0.484434859503982 -0.0421075646687726 0.0707729878439456 0.194652734994256 0.228245996325073 0.0205868707940657 0.265375300831334 -0.0690865561210724 0.171020836332001 -0.109255798230587 -0.342291026212677 0.246600324477957 -0.286634258339929 0.414463565648986 -0.0320693762724592 -0.0549260387377591 -0.0763052560607504 -0.316517703668002 0.102454347321762 0.259951549151108 0.000779138183330197 0.335441559913048 -0.300382901895126 0.0247619927333676 -0.374133509310235 0.258017032268099 0.444763273538309 -0.247418831402322 -0.167344671329699 0.0792795449141843 0.359406927055753 0.307667982259349))

; scale particle coordinates by multiplying by mag
; (define xx (vector-map (lambda (n) (* mag n)) xx))
; (define yy (vector-map (lambda (n) (* mag n)) yy))

; source properties
(define-param fcen 2.4)
(define-param df 4)
(define-param nfreq 200)

; geometry
(define-param sx0 1) ; width of original image
(define-param sy0 1) ; height of original image

(if periodic?  ; if not periodic, adds space for absorbing layers on all sides
  (define sx (* sx0 mag))
  (define sx (+ (* sx0 mag) (* 2 dpml))))

(define sy (+ (* sy0 mag) pad cor sub (* 2 dpml)))

(set! geometry-lattice (make lattice (size sx sy no-size)))


; material properties
(define-param n_matrix 1)           ; refractive index of matrix
(define-param n_rod 2)              ; RI of rods
(define-param n_cortex 1.56)        ; RI of cortex
(define-param n_out 1)              ; RI of outside material (air usually)
(define eps1 (* n_matrix n_matrix)) ; permittivity of material 1
(define eps2 (* n_rod n_rod))       ; permittivity of material 2
(define eps3 (* n_cortex n_cortex)) ; permittivity of material 3
(define eps4 (* n_out n_out))       ; permittivity of material 3
(define-param k_matrix 0)           ; extinction coeff mat 1
(define-param k_rod 0)              ; extinction coeff mat 2
(define-param k_cortex 0)           ; extinction coeff mat 2
(define-param k_out 0)              ; extinction coeff mat 2
(define matrix (make medium 
                (epsilon eps1) 
                  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps1) k_matrix)) eps1))))
(define rod (make medium 
                (epsilon eps2) 
                  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps2) k_rod)) eps2))))
(define cortex (make medium 
                (epsilon eps3) 
                  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps3) k_cortex)) eps3))))
(define outside (make medium 
                (epsilon eps4) 
                  (D-conductivity (/ (* 2 pi fcen (* 2 (sqrt eps4) k_out)) eps4))))
(set! default-material matrix)


; loop to place rods at specified locations
(define (doit x x-max dx)
	(if (<= x x-max)
 	  (begin
  		(set! geometry
        (append geometry
  			   (list
  			     (make cylinder (center 
                                (+ (* (vector-ref xx x) (cos rot) mag) (* (vector-ref yy x) mag (sin rot)))
                                (+ (- (* 0.5 sy) dpml pad cor (/ (* sy0 mag) 2)) (+ (* (* (vector-ref xx x) -1) (sin rot) mag) (* (vector-ref yy x) mag (cos rot)))))
  				 	                (radius r)
                            (height infinity) (material rod))
             (make cylinder (center 
                                (+ (* (vector-ref xx x) (cos rot) mag) (* (vector-ref yy x) mag (sin rot)))
                                (+ (- (* 0.5 sy) dpml pad cor (/ (* sy0 mag) 2)) 
                                   (+ (* (* (vector-ref xx x) -1) (sin rot) mag) (* (vector-ref yy x) mag (cos rot)))))
                            (radius r2)
                            (height infinity) (material air)))))
  	     (doit (+ x dx) x-max dx))))

; define geometry
(if no-struc?
	(set! geometry 
   (list 
		(make block (center 0 0 0) 
                (size sx sy no-size)
                (material outside))))
	(set! geometry
    (list 
      (make block (center 0 (- (* 0.5 sy) (/ (+ dpml pad) 2)))
                  (size sx (+ pad cor))
                  (material outside))
      (make block (center 0 (+ (- (* 0.5 sy) dpml pad (/ cor 2) (/ bleed 2)) r))
                  (size sx (+ cor bleed) infinity) 
                  (material cortex)))))

(if (not no-struc?)
  (doit 0 (- (vector-length xx) 1) 1))

; boundary conditions, etc.
(if periodic? 
  (set! pml-layers (list (make pml (direction Y) (thickness dpml)))) ; only on Y extremes
  (set! pml-layers (list (make pml (thickness dpml))))) ; all around cell

(define-param theta 0)   ; inc angle in degrees
(define theta_rad (* pi (/ theta 180)))  ; angle of incidence (with respect to y-axis)
(define k (* fcen (sin theta_rad)))
(set! k-point (vector3 k 0 0))

; (set! ensure-periodicity true)

; are these parameters necessary? {
; (set! ensure-periodicity true)
; (set! k-point (vector3 0 0 0))
; }

(define-param res 50)

(set! resolution res)

; define custom amp function
(define (my-amp-func p) 
   (exp (* 0+2i pi k (vector3-x p))))

; excite source
(set! sources 
  (list
   (make source
    (src (make gaussian-src (frequency fcen) (fwidth df)))
    (if TE?
      (component Hz)
      (component Ez))
    (center 0 (- (* -0.5 sy) dpml 1)) (size sx 0)
    (amp-func my-amp-func))))

; reflected flux
(define refl
	(add-flux fcen df nfreq
	 (make flux-region
	  (center 0 (- (* 0.5 sy) dpml 2)) (size (* sx 2) 0))))
    ; (center 0 (- (* 0.5 sy) dpml 2)) (size ap_width 0))))

; transmitted flux
(define trans
	(add-flux fcen df nfreq
	 (make flux-region
	  (center 0 (+ (* -0.5 sy) dpml 1)) (size (* sx 2) 0))))
    ; (center 0 (+ (* -0.5 sy) dpml 1)) (size ap_width 0))))

(if (not no-struc?) (load-minus-flux "refl-flux" refl))

; run sources
(if TE?
  (run-sources+ 
    (stop-when-fields-decayed 50 Hz
      (vector3 0 (- (* 0.5 sy) dpml 1)) 1e-2)
       (at-beginning output-epsilon)
          (at-end (output-png Hz "-Zc bluered")))
  (run-sources+ 
    (stop-when-fields-decayed 50 Ez
      (vector3 0 (- (* 0.5 sy) dpml 1)) 1e-2)
       (at-beginning output-epsilon)
          (at-end (output-png Ez "-Zc bluered"))))

; run source to get fields in x plane
; (if no-struc?
;   (run-until 200 
;     (at-beginning output-epsilon)
;       (to-appended "hz_nostruc"
;         (in-volume (volume (center 0 2) (size (* sx 2) 0))
;         output-hfield))))

; (if (not no-struc?)
;   (run-until 200 
;     (at-beginning output-epsilon)
;       (to-appended "hz_struc"
;         (in-volume (volume (center 0 -2) (size (* sx 2) 0))
;         output-hfield))))

; uncomment this to test the script
; (run-until 10
;   (at-beginning output-epsilon)
;     (at-every 0.6 (output-png Hz "-Zc bluered")))

; output data

(if no-struc? (save-flux "refl-flux" refl))

(display-fluxes trans refl) ; print out the flux spectrum

; TEMP - compute local density of states
; (define-param Th 500)
; (run-sources+ Th (after-sources (harminv Ez (vector3 0) fcen df)))
; (define f (harminv-freq-re (car harminv-results)))
; (define Q (harminv-Q (car harminv-results)))
; (define Vmode (* 0.25 sx sy))
; (print "ldos0:, " (/ Q Vmode (* 2 pi f pi 0.5)))
