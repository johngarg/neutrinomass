(import [adderall.dsl [*]])
(import [hydiomatic.core [*]])
(import [adderall.bitnum [+o -o]])

(require [adderall.dsl [*]])
(require [hydiomatic.macros [*]])

;; define macro for dealing with combinations
(defrules [operators]
  ;; O1
  [`(* (L ~?a ~?i) (L ~?b ~?j) (H ~?k) (H ~?l) (eps ~?a ~?b)
       (* (eps ~?i ~?k) (eps ~?j ~?l)))
   'O1]
  [`(* (L ~?a ~?i) (L ~?b ~?j) (H ~?k) (H ~?l) (eps ~?a ~?b)
       (* (eps ~?j ~?l) (eps ~?i ~?k)))
   'O1]
  [`(* (L ~?a ~?i) (L ~?b ~?j) (H ~?l) (H ~?k) (eps ~?a ~?b)
       (* (eps ~?i ~?k) (eps ~?j ~?l)))
   'O1]
  [`(* (L ~?a ~?i) (L ~?b ~?j) (H ~?l) (H ~?k) (eps ~?a ~?b)
       (* (eps ~?j ~?l) (eps ~?i ~?k)))
   'O1]

  ;; O3
  [`(* (L ~?a ~?i) (L ~?b ~?j) (Q ~?c ~?k) (db ~?d) (H ~?l)
       (* (eps ~?a ~?c) (eps ~?b ~?d))
       (* (eps ~?i ~?j) (eps ~?k ~?l)))
   'O3a]
  [`(* (L ~?a ~?i) (L ~?b ~?j) (Q ~?c ~?k) (db ~?d) (H ~?l)
       (* (eps ~?a ~?c) (eps ~?b ~?d))
       (* (eps ~?i ~?k) (eps ~?j ~?l)))
   'O3b]
  ;; [`(* (eps i j) (eps k l)) 'O3a]
  ;; [`(* (eps i l) (eps k j)) 'O3b]
  ;; [`(* (eps k j) (eps i l)) 'O3b]
  ;; [`(* (eps i l) (eps j k)) '(- O3b)]
  ;; [`(* (eps j k) (eps i l)) '(- O3b)]
  ;; [`(* (eps i j) (eps k l)) 'O3a]
  ;; [`(* (eps i l) (eps k j)) 'O3b]
  ;; [`(* (eps k j) (eps i l)) 'O3b]
  ;; [`(* (eps i l) (eps j k)) '(- O3b)]
  ;; [`(* (eps j k) (eps i l)) '(- O3b)]
  )

;; TODO Idea, unify what you have immediately with the patterns you are matching
;; for, if a set of epsilons match, only work with the other set

(defrules [epsilon-rules]
  ;; Higgs symmetry
  (prep
    (appendo ?xs (cons `(H ~?i) `(H ~?j) ?rest) expr)
    (== ?xs (cons '* ?r))
    (appendo ?ys
             (cons `(* ~(cons 'eps1 ?left-indices)
                       ~(cons 'eps1 ?right-indices))
                   ?zs)
             ?rest)
    (conde
      [(== ?left-indices `(~?i ~?j)) (== out 0)]
      [(== ?left-indices `(~?j ~?i)) (== out 0)]
      [(== ?right-indices `(~?i ~?j)) (== out 0)]
      [(== ?right-indices `(~?j ~?i)) (== out 0)]))

  ;; Schouten identity
  [`(* (eps ~?a ~?b) (eps ~?c ~?d)) `(+ (* (eps1 ~?a ~?c) (eps1 ~?b ~?d))
                                        (* (eps1 ~?a ~?d) (eps1 ~?c ~?b)))]

  ;; [(cons ?X `(* (eps ~?a ~?b) (eps ~?c ~?d)) ?XS) `(~?X (+ (* (eps ~?a ~?c) (eps ~?b ~?d)
  ;;                                                          (* (eps ~?a ~?d) (eps ~?b ~?c)))) ~?XS)]
  ;; [`(~?y (* (eps ~?a ~?b) (eps ~?c ~?d))) ?y]

  ;; Higgs antisymmetry


;; (+ (* (L a i) (L b j) (H k) (H l) (eps a b) (* (eps1 j k) (eps1 i l)))
;;    (* (L a i) (L b j) (H k) (H l) (eps a b) (* (eps1 j l) (eps1 k i))))
 )

(defrules [arithmetic]
  ;; (+ 0 x), (+ x 0) => x
  [`(+ 0 ~?x) `(+ ~?x 0)]
  [`(+ ~?x 0) ?x]

  ;; (* 1 x), (* x 1) => x
  [`(* (n 0) ~?x) ?x]
  [`(* ~?x (n 0)) ?x]

  ;; (* (+ 1 2) 5) => (+ (* 1 5) (* 2 5))
  [`(* (+ ~?x ~?y) ~?z) `(+ (* ~?x ~?z) (* ~?y ~?z))]
  [`(* ~?z (+ ~?x ~?y)) `(+ (* ~?z ~?x) (* ~?z ~?y))]

  ;; (+ a a) => (* 2 a)
  ;; [`(+ ~?x ~?x) `(* (n (n 0)) ~?x)]

  ;; (+ (* 2 a) a) => (* 3 a)
  (prep
    (== expr `(+ (* (n ~?n) ~?a) ~?a))
    (== out `(* (n (n ~?n)) ~?a)))

  ;; (+ a (* 2 a)) => (* 3 a)
  (prep
    (== expr `(+ ~?a (* (n ~?n) ~?a)))
    (== out `(* (n (n ~?n)) ~?a)))

  ;; (* ... (+ a b) ...) => (+ (* ... a ...) (* ... b ...))
  (prep
    (appendo (cons '* ?xs) (cons `(+ ~?a ~?b) ?rest) expr) ;; this is significant and useful
    (appendo (cons '* ?xs) (cons ?a ?rest) ?left)
    (appendo (cons '* ?xs) (cons ?b ?rest) ?right)
    (== out `(+ ~?left ~?right)))

  ;; (+ x (+ ...)) => (+ x ...)
  ;; [`(+ ~?x ~(cons '+ ?xs)) (cons '+ ?x ?xs)]

  ;; (* x (* ...)) => (* x ...)
  ;; [`(* ~?x ~(cons '* ?xs)) (cons '* ?x ?xs)]
  ;; [`(* ~(cons '* ?xs) ~?x) (cons '* ?x ?xs)]
)
