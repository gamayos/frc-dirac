# Supplementary Checks

Run:

```bash
python3 src/finite_checks.py
```

Expected outputs:

- exact image counts and loss factors for several power maps on `F_13^x`,
- verification that the Cayley propagator in the Schr\"odinger example satisfies `U^13 = I`,
- exact admissibility and order classification for all trace-zero parameters `alpha=k*c` in the `p=13` Schr\"odinger example,
- verification of the Clifford relations, the boost conjugation formulas,
  the transported gamma formulas, and an exact sample covariance check in
  the Dirac example.

The script uses only exact arithmetic over `F_13[c]/(c^2-2)` and does not rely on floating-point computation.

Current exact output:

```text
Using p=13, nu=2, i_t=5, c^2=2
Power map counts on F_p^x
  epsilon= 1: image_size=12, formula=12, loss_factor= 1
  epsilon= 2: image_size= 6, formula= 6, loss_factor= 2
  epsilon= 3: image_size= 4, formula= 4, loss_factor= 3
  epsilon= 4: image_size= 3, formula= 3, loss_factor= 4
  epsilon= 6: image_size= 2, formula= 2, loss_factor= 6
  epsilon=12: image_size= 1, formula= 1, loss_factor=12
Schr\"odinger check
  det-denominator-invertible = True
  U^13 = I                 = True
  admissible alpha=k*c and exact orders
    k= 0: admissible=True, order=1
    k= 1: admissible=True, order=13
    k= 2: admissible=True, order=13
    k= 3: admissible=True, order=13
    k= 4: admissible=True, order=13
    k= 5: admissible=True, order=13
    k= 6: admissible=True, order=13
    k= 7: admissible=True, order=13
    k= 8: admissible=True, order=13
    k= 9: admissible=True, order=13
    k=10: admissible=True, order=13
    k=11: admissible=True, order=13
    k=12: admissible=True, order=13
Dirac check
  Clifford relations       = True
  boost formula gamma^0    = True
  boost formula gamma^1    = True
  transported gamma^0      = True
  transported gamma^1      = True
  covariance sample check  = True
  A coefficient            = Fp2(10,0)
  B coefficient            = Fp2(2,0)
```
