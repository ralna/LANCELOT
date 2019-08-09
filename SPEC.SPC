BEGIN
* maximizer-sought
  check-derivatives
  ignore-derivative-bugs
* use-scalings
* print-scalings
* finite-difference-gradients
  exact-second-derivatives-used
* bfgs-approximate-second-derivatives-used
* sr1-approximate-second-derivatives-used
  bandsolver-preconditioned-cg-solver-used   5
* diagonal-preconditioned-cg-solver-used
* gill-murray-ponceleon-saunders-preconditioned-cg-solver-used
* schnabel-eskow-preconditioned-cg-solver-used
* munksgaards-preconditioned-cg-solver-used
  exact-cauchy-point-required
* inexact-cauchy-point-required
* solve-bqp-accurately
* two-norm-trust-region
* gradient-tolerance    1.0D-5
* constraint-tolerance  1.0D-5
  trust-region-radius   1.0D+0
  maximum-number-of-iterations   1000
  print-level                    0
  start-printing-at-iteration    0
  stop-printing-at-iteration     1000
END
