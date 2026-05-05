# Composite Numerical Quadrature

> A priori and adaptive composite Trapezoidal rules versus the composite two-point Gauss–Legendre rule, implemented in C++.

**Programming Assignment 4 — Foundations of Computational Mathematics 2 (Spring 2026)**
Godfred Antwi Koduah · Florida State University

---

## Overview

This project implements and empirically compares three composite numerical quadrature methods on five test integrals with known closed-form values. All three methods share a common function-evaluation counter so the cost comparison reflects real work done.

| # | Method | Convergence | Per-step cost |
|---|---|---|---|
| 1 | Composite Trapezoidal Rule with **a priori** chosen `H_m` | `O(H²)` | `n + 1` evals |
| 2 | Composite Trapezoidal Rule with **adaptive step-halving** (α = 1/2) | `O(H²)` | `2^L + 1` evals after `L` levels |
| 3 | Composite **Two-point Gauss–Legendre** Rule | `O(H⁴)` | `2n` evals |

The adaptive scheme uses the recursive doubling identity

```
T_{2n} = ½·T_n + H_new · Σ f(midpoints)
```

so each refinement adds only `n` new function evaluations, and stops when the Richardson estimate `|T_{2n} − T_n| / 3` falls below the requested tolerance.

---

## Test Integrals

| ID | Integrand | Interval | Notable feature |
|---|---|---|---|
| (1) | `eˣ` | `[0, 3]` | smooth, monotone |
| (2) | `e^{sin(2x)} · cos(2x)` | `[0, π/3]` | smooth, oscillatory product |
| (3) | `tanh(x)` | `[−2, 1]` | localized 2nd derivative |
| (4) | `x · cos(2πx)` | `[0, 3.5]` | highly oscillatory (~7 periods) |
| (5) | `x + 1/x` | `[0.1, 2.5]` | near-singular at left endpoint |

Closed-form values are hard-coded in `build_tests()` for true-error reporting.

---

## Build & Run

Requires a C++14 compiler (or newer).

```bash
g++ -O2 -std=c++14 main.cpp -o quadrature
./quadrature
```

The program creates a `data/` directory and writes all results as CSV.

---

## Output Files

| File | Contents |
|---|---|
| `data/integrals_info.csv` | Per-integral metadata: `[a, b]`, exact value, numerically estimated `max\|f''\|` and `max\|f⁽⁴⁾\|` |
| `data/task_ii_trap_fixed.csv` | Task (i)+(ii): a priori `H_m` and fixed-`H` Trap results across `ε_rel ∈ {1e−2, …, 1e−9}` |
| `data/task_iii_trap_adaptive.csv` | Task (iii): adaptive step-halving Trap results |
| `data/task_iv_gauss.csv` | Task (iv): composite 2-pt Gauss–Legendre swept over `n = 1, 2, 4, …, 32768` |
| `data/trap_fixed_sweep.csv` | Fixed-`H` Trap swept up to `n = 2¹⁸` for the cross-method cost plot |
| `data/test_polynomial_exactness.csv` | Verifies Trap is exact for deg ≤ 1 and Gauss for deg ≤ 3 |
| `data/test_doubling_identity.csv` | Verifies `T_{2n}` from the recursion matches direct evaluation to machine precision |
| `data/test_empirical_rates.csv` | Observed convergence rates `p̂ = log₂(E(n)/E(2n))` |

---

## Key Results

**Convergence rates match theory.** Empirical `p̂ → 2` (Trap) and `p̂ → 4` (Gauss) for every integrand once `n` is large enough, agreeing with theory to better than 0.5%.

**A priori bound is tight in slope but conservative in constant.** The bound/error ratio is roughly 3 for smooth integrands (1)–(2), grows to ≈ 50 for the near-singular (5), and reaches ≈ 240 for the oscillatory (4) where adjacent subintervals exhibit error cancellation that the pointwise `max\|f''\|` cannot see.

**Adaptive step-halving ≈ a priori sizing in cost.** The two Trap variants are nearly indistinguishable in accuracy-per-evaluation. The adaptive scheme's real advantage is **robustness**: it needs no derivative bound, only a tolerance.

**Gauss–Legendre dominates for smooth integrands.** Empirical speedup grows as `ε^{−1/4}`:

| Target `\|err\|` | `N_Trap / N_Gauss` (integral 1) |
|---|---|
| `10⁻²` | ~8× |
| `10⁻⁴` | 32× |
| `10⁻⁶` | 64× |
| `10⁻⁸` | 256× |

**Round-off saturation.** Gauss reaches the round-off floor `~10⁻¹⁴ · \|I\|` with `n ≤ 1024`; the Trap variants cannot practically reach `10⁻¹⁰` within the tested budget of `n = 2¹⁸`.

**Bottom line:** *Basic-rule quality dominates adaptive cleverness.* A degree-3 basic rule (Gauss 2-pt) beats a degree-1 basic rule (Trap) by a factor of `H²` per step, and no amount of adaptive refinement of a low-order rule recovers that gap. Adaptive refinement is an orthogonal improvement best layered on top of a high-order basic rule.

---

## File Structure

```
.
├── main.cpp                 # All three quadrature codes + tests + driver
├── README.md                # This file
├── report.pdf               # Full write-up with figures and tables
└── data/                    # Generated at runtime
    ├── integrals_info.csv
    ├── task_ii_trap_fixed.csv
    ├── task_iii_trap_adaptive.csv
    ├── task_iv_gauss.csv
    ├── trap_fixed_sweep.csv
    ├── test_polynomial_exactness.csv
    ├── test_doubling_identity.csv
    └── test_empirical_rates.csv
```

---

## References

1. K. E. Atkinson, *An Introduction to Numerical Analysis*, 2nd ed., Wiley, 1989.
2. P. J. Davis and P. Rabinowitz, *Methods of Numerical Integration*, 2nd ed., Academic Press, 1984.
3. W. H. Press et al., *Numerical Recipes*, 3rd ed., Cambridge University Press, 2007.
4. L. N. Trefethen, *Approximation Theory and Approximation Practice*, SIAM, 2013.
5. K. A. Gallivan, Lecture notes for FCM2, FSU, Spring 2026.
