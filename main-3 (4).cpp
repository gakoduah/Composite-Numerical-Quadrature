// =============================================================================
// Programming Assignment 4
// Foundations of Computational Mathematics 2 - Spring 2026
// Godfred Antwi Koduah
//
// Composite Numerical Quadrature:
//   (1) Composite Trapezoidal Rule with a priori chosen H_m
//   (2) Composite Trapezoidal Rule with adaptive step halving (alpha = 1/2)
//   (3) Composite Two-point Gauss-Legendre Rule
//
// =============================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <algorithm>
#include <sys/stat.h>

using Func = std::function<double(double)>;

// ---------------------------------------------------------------------------
// Result structure
// ---------------------------------------------------------------------------
struct QuadResult {
    double value;      // computed integral value
    double true_err;   // absolute true error
    double rel_err;    // relative true error
    long   feval;      // number of function evaluations
    double H;          // subinterval width
    long   n;          // number of subintervals
};

// ---------------------------------------------------------------------------
// Composite Trapezoidal Rule with a priori choice of n subintervals.
// Uses n+1 function evaluations: f(a), f(a+H), ..., f(b).
// ---------------------------------------------------------------------------
QuadResult comp_trap_fixed(const Func& f, double a, double b, long n) {
    double H = (b - a) / static_cast<double>(n);
    double s = 0.5 * (f(a) + f(b));
    for (long i = 1; i < n; ++i) {
        s += f(a + i * H);
    }
    QuadResult r;
    r.value = s * H;
    r.feval = n + 1;
    r.H = H;
    r.n = n;
    return r;
}

// ---------------------------------------------------------------------------
// Composite Trapezoidal Rule with adaptive step halving (alpha = 1/2).
//
// Starts with n = 1 (T_0 uses 2 evaluations: f(a), f(b)).
// Each refinement level halves H (doubles n), and reuses previously computed
// points: the new trapezoidal value is
//
//     T_{2n} = 0.5 * T_n + H_new * sum_{i=0}^{n-1} f(a + (i+0.5)*H_old)
//
// Richardson-based error estimate (Trapezoidal error is O(H^2)):
//     |I - T_{2n}| ~= |T_{2n} - T_n| / 3
//
// Termination when estimated error is below the requested absolute tolerance.
// ---------------------------------------------------------------------------
QuadResult comp_trap_adaptive(const Func& f, double a, double b,
                              double abstol, int max_levels = 28,
                              double* err_est_out = nullptr)
{
    double fa = f(a), fb = f(b);
    long feval = 2;
    long n = 1;
    double H = b - a;
    double T_old = 0.5 * H * (fa + fb);   // T_0
    double T_new = T_old;
    double err_est = std::numeric_limits<double>::infinity();

    for (int level = 1; level <= max_levels; ++level) {
        double H_old = H;
        H = H * 0.5;
        double sum_mid = 0.0;
        // add midpoints of the n current subintervals
        for (long i = 0; i < n; ++i) {
            double x = a + (i + 0.5) * H_old;
            sum_mid += f(x);
            ++feval;
        }
        T_new = 0.5 * T_old + H * sum_mid;
        n *= 2;

        err_est = std::abs(T_new - T_old) / 3.0;
        // Require at least 2 refinements before accepting the estimate.
        if (level >= 2 && err_est < abstol) {
            break;
        }
        T_old = T_new;
    }
    if (err_est_out) *err_est_out = err_est;
    QuadResult r;
    r.value = T_new;
    r.feval = feval;
    r.H = H;
    r.n = n;
    return r;
}

// ---------------------------------------------------------------------------
// Composite Two-Point Gauss-Legendre Rule.
// Basic rule on [-1, 1]:  nodes = +- 1/sqrt(3), weights = 1, 1.
// Mapped to subinterval [x_i, x_{i+1}] of width H:
//     nodes  = midpoint +- (H/2) * (1/sqrt(3))
//     weight = H/2  for each node
// Uses 2n function evaluations, degree of exactness 3 (exact for cubics).
// ---------------------------------------------------------------------------
QuadResult comp_gauss2(const Func& f, double a, double b, long n) {
    const double inv_sqrt3 = 1.0 / std::sqrt(3.0);
    double H = (b - a) / static_cast<double>(n);
    double half_H = 0.5 * H;
    double s = 0.0;
    for (long i = 0; i < n; ++i) {
        double xm = a + (i + 0.5) * H;
        double x1 = xm - half_H * inv_sqrt3;
        double x2 = xm + half_H * inv_sqrt3;
        s += f(x1) + f(x2);
    }
    QuadResult r;
    r.value = s * half_H;
    r.feval = 2 * n;
    r.H = H;
    r.n = n;
    return r;
}

// ---------------------------------------------------------------------------
// Numerical estimate of max|f^(k)| over [a,b] using a finite-difference
// approximation of the k-th derivative on a fine uniform grid.
// Used to obtain a priori error constants without symbolic differentiation.
// ---------------------------------------------------------------------------
// Max |f''| via a 3-point central difference on a fine grid of interior points.
double estimate_max_second_deriv(const Func& f, double a, double b, int N = 4001) {
    double step = std::min(1e-4, (b - a) / 2000.0);
    double x_start = a + step;
    double x_end   = b - step;
    double h = (x_end - x_start) / (N - 1);
    double m = 0.0;
    for (int i = 0; i < N; ++i) {
        double x = x_start + i * h;
        double d2 = (f(x + step) - 2.0 * f(x) + f(x - step)) / (step * step);
        m = std::max(m, std::abs(d2));
    }
    // Also include the endpoints via one-sided 3-point formula
    double d2_a = (f(a + 2*step) - 2*f(a + step) + f(a)) / (step * step);
    double d2_b = (f(b) - 2*f(b - step) + f(b - 2*step)) / (step * step);
    m = std::max(m, std::max(std::abs(d2_a), std::abs(d2_b)));
    return m;
}

// Max |f^(4)| via a 5-point central difference on a fine grid of interior points.
// step is chosen to balance truncation and round-off error.
double estimate_max_fourth_deriv(const Func& f, double a, double b, int N = 4001) {
    double step = std::min(5e-3, (b - a) / 500.0);
    double x_start = a + 2.0 * step;
    double x_end   = b - 2.0 * step;
    if (x_end <= x_start) {
        step = (b - a) / 20.0;
        x_start = a + 2.0 * step;
        x_end   = b - 2.0 * step;
    }
    double h = (x_end - x_start) / (N - 1);
    double s2 = step * step;
    double s4 = s2 * s2;
    double m = 0.0;
    for (int i = 0; i < N; ++i) {
        double x = x_start + i * h;
        double d4 = (f(x + 2*step) - 4.0*f(x + step) + 6.0*f(x)
                    - 4.0*f(x - step) + f(x - 2*step)) / s4;
        m = std::max(m, std::abs(d4));
    }
    return m;
}

// ---------------------------------------------------------------------------
// A priori H_m for Composite Trapezoidal Rule to achieve |I - T_n| <= eps_abs.
//   |I - T_n| <= (b - a) * H^2 / 12 * max|f''|
//   H <= sqrt( 12 * eps_abs / ((b-a) * max|f''|) )
// For relative error eps_rel: eps_abs = eps_rel * |I_exact|.
// ---------------------------------------------------------------------------
double apriori_H_trap(double a, double b, double abs_tol, double max_f2) {
    if (max_f2 <= 0.0) return (b - a);
    return std::sqrt(12.0 * abs_tol / ((b - a) * max_f2));
}

// A priori error bound (upper estimate) for Trapezoidal with given H.
double apriori_bound_trap(double a, double b, double H, double max_f2) {
    return (b - a) * H * H * max_f2 / 12.0;
}

// A priori error bound for 2-point Gauss-Legendre with given H.
//   |I - G_n| <= (b-a) * H^4 / 4320 * max|f^(4)|
double apriori_bound_gauss(double a, double b, double H, double max_f4) {
    return (b - a) * std::pow(H, 4) * max_f4 / 4320.0;
}

// ---------------------------------------------------------------------------
// Test integrals
// ---------------------------------------------------------------------------
struct TestIntegral {
    std::string name;
    std::string latex;
    Func        f;
    double      a, b;
    double      exact;
};

std::vector<TestIntegral> build_tests() {
    std::vector<TestIntegral> T;

    T.push_back({
        "I1_exp",
        "int_0^3 e^x dx = e^3 - 1",
        [](double x) { return std::exp(x); },
        0.0, 3.0,
        std::exp(3.0) - 1.0
    });

    T.push_back({
        "I2_exp_sin",
        "int_0^{pi/3} e^{sin(2x)} cos(2x) dx = 0.5*(-1 + e^{sqrt(3)/2})",
        [](double x) { return std::exp(std::sin(2.0 * x)) * std::cos(2.0 * x); },
        0.0, M_PI / 3.0,
        0.5 * (-1.0 + std::exp(std::sqrt(3.0) / 2.0))
    });

    T.push_back({
        "I3_tanh",
        "int_{-2}^1 tanh(x) dx = ln(cosh(1)/cosh(2))",
        [](double x) { return std::tanh(x); },
        -2.0, 1.0,
        std::log(std::cosh(1.0) / std::cosh(2.0))
    });

    T.push_back({
        "I4_xcos",
        "int_0^{3.5} x*cos(2*pi*x) dx = -1/(2*pi^2)",
        [](double x) { return x * std::cos(2.0 * M_PI * x); },
        0.0, 3.5,
        -1.0 / (2.0 * M_PI * M_PI)
    });

    T.push_back({
        "I5_xinvx",
        "int_{0.1}^{2.5} (x + 1/x) dx = (2.5^2 - 0.1^2)/2 + ln(25)",
        [](double x) { return x + 1.0 / x; },
        0.1, 2.5,
        (2.5 * 2.5 - 0.1 * 0.1) / 2.0 + std::log(2.5 / 0.1)
    });

    return T;
}

// ---------------------------------------------------------------------------
// Helper: make sure output directory exists
// ---------------------------------------------------------------------------
void ensure_dir(const std::string& d) {
    mkdir(d.c_str(), 0775);
}

// ---------------------------------------------------------------------------
//                            MAIN  DRIVER
// ---------------------------------------------------------------------------
int main() {
    ensure_dir("data");

    auto tests = build_tests();

    // Target relative errors (ranging from loose to tight)
    std::vector<double> eps_rel_list = {
        1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9
    };

    // -------------------------------------------------------------
    // Write summary information about each integral
    // -------------------------------------------------------------
    {
        std::ofstream fout("data/integrals_info.csv");
        fout << "name,a,b,exact,max_f2,max_f4\n";
        fout << std::setprecision(16);
        for (auto& T : tests) {
            double M2 = estimate_max_second_deriv(T.f, T.a, T.b);
            double M4 = estimate_max_fourth_deriv(T.f, T.a, T.b);
            fout << T.name << "," << T.a << "," << T.b << "," << T.exact
                 << "," << M2 << "," << M4 << "\n";
        }
    }

    // =====================================================================
    // TASK (i), (ii):  A priori H_m and Fixed-H Composite Trapezoidal Rule
    // =====================================================================
    std::cout << "\n===== TASK (i), (ii): A priori H_m + Fixed-H Composite "
                 "Trapezoidal =====\n\n";

    {
        std::ofstream fout("data/task_ii_trap_fixed.csv");
        fout << "integral,eps_rel,eps_abs,max_f2,H_predicted,n_used,H_used,"
                "approx,true_err,apriori_bound,feval\n";
        fout << std::setprecision(16);

        for (auto& T : tests) {
            double M2 = estimate_max_second_deriv(T.f, T.a, T.b);
            std::cout << "Integral " << T.name
                      << "  exact = " << std::setprecision(12) << T.exact
                      << "  max|f''| ~= " << M2 << "\n";
            std::cout << std::setw(10) << "eps_rel"
                      << std::setw(14) << "H_pred"
                      << std::setw(10) << "n"
                      << std::setw(16) << "true_err"
                      << std::setw(16) << "bound"
                      << std::setw(10) << "feval" << "\n";

            for (double eps : eps_rel_list) {
                double abs_tol = eps * std::abs(T.exact);
                double H_pred = apriori_H_trap(T.a, T.b, abs_tol, M2);
                long n = std::max<long>(1, (long)std::ceil((T.b - T.a) / H_pred));
                QuadResult r = comp_trap_fixed(T.f, T.a, T.b, n);
                double true_err = std::abs(r.value - T.exact);
                double bound = apriori_bound_trap(T.a, T.b, r.H, M2);

                fout << T.name << "," << eps << "," << abs_tol << "," << M2
                     << "," << H_pred << "," << n << "," << r.H
                     << "," << r.value << "," << true_err << "," << bound
                     << "," << r.feval << "\n";

                std::cout << std::setw(10) << eps
                          << std::setw(14) << H_pred
                          << std::setw(10) << n
                          << std::setw(16) << true_err
                          << std::setw(16) << bound
                          << std::setw(10) << r.feval << "\n";
            }
            std::cout << "\n";
        }
    }

    // =====================================================================
    // TASK (iii):  Adaptive Step-Halving Composite Trapezoidal
    // =====================================================================
    std::cout << "\n===== TASK (iii): Adaptive Step-Halving Trapezoidal =====\n\n";

    {
        std::ofstream fout("data/task_iii_trap_adaptive.csv");
        fout << "integral,eps_rel,eps_abs,n_final,H_final,"
                "approx,true_err,err_est,feval\n";
        fout << std::setprecision(16);

        for (auto& T : tests) {
            std::cout << "Integral " << T.name << "\n";
            std::cout << std::setw(10) << "eps_rel"
                      << std::setw(10) << "n"
                      << std::setw(16) << "true_err"
                      << std::setw(16) << "err_est"
                      << std::setw(10) << "feval" << "\n";

            for (double eps : eps_rel_list) {
                double abs_tol = eps * std::abs(T.exact);
                double err_est = 0.0;
                QuadResult r = comp_trap_adaptive(T.f, T.a, T.b, abs_tol, 28,
                                                  &err_est);
                double true_err = std::abs(r.value - T.exact);

                fout << T.name << "," << eps << "," << abs_tol
                     << "," << r.n << "," << r.H
                     << "," << r.value << "," << true_err
                     << "," << err_est << "," << r.feval << "\n";

                std::cout << std::setw(10) << eps
                          << std::setw(10) << r.n
                          << std::setw(16) << true_err
                          << std::setw(16) << err_est
                          << std::setw(10) << r.feval << "\n";
            }
            std::cout << "\n";
        }
    }

    // =====================================================================
    // TASK (iv):  Composite Two-Point Gauss-Legendre, sweep over H_m
    // =====================================================================
    std::cout << "\n===== TASK (iv): Composite 2-point Gauss-Legendre =====\n\n";

    {
        std::ofstream fout("data/task_iv_gauss.csv");
        fout << "integral,n,H,approx,true_err,apriori_bound,feval\n";
        fout << std::setprecision(16);

        // Sweep n = 1, 2, 4, 8, ..., up to 2^15 = 32768
        std::vector<long> n_list;
        for (int k = 0; k <= 15; ++k) n_list.push_back(1L << k);

        for (auto& T : tests) {
            double M4 = estimate_max_fourth_deriv(T.f, T.a, T.b);
            std::cout << "Integral " << T.name
                      << "  max|f^(4)| ~= " << M4 << "\n";
            std::cout << std::setw(8)  << "n"
                      << std::setw(14) << "H"
                      << std::setw(16) << "true_err"
                      << std::setw(16) << "bound"
                      << std::setw(10) << "feval" << "\n";

            for (long n : n_list) {
                QuadResult r = comp_gauss2(T.f, T.a, T.b, n);
                double true_err = std::abs(r.value - T.exact);
                double bound = apriori_bound_gauss(T.a, T.b, r.H, M4);
                fout << T.name << "," << n << "," << r.H
                     << "," << r.value << "," << true_err
                     << "," << bound << "," << r.feval << "\n";

                std::cout << std::setw(8)  << n
                          << std::setw(14) << r.H
                          << std::setw(16) << true_err
                          << std::setw(16) << bound
                          << std::setw(10) << r.feval << "\n";
            }
            std::cout << "\n";
        }
    }

    // =====================================================================
    // Extended: also sweep a range of n for fixed-H Trap, for convergence
    // study (useful for plotting error vs feval comparing all 3 methods)
    // =====================================================================
    {
        std::ofstream fout("data/trap_fixed_sweep.csv");
        fout << "integral,n,H,approx,true_err,feval\n";
        fout << std::setprecision(16);
        std::vector<long> n_list;
        for (int k = 0; k <= 18; ++k) n_list.push_back(1L << k);
        for (auto& T : tests) {
            for (long n : n_list) {
                QuadResult r = comp_trap_fixed(T.f, T.a, T.b, n);
                double true_err = std::abs(r.value - T.exact);
                fout << T.name << "," << n << "," << r.H
                     << "," << r.value << "," << true_err
                     << "," << r.feval << "\n";
            }
        }
    }

    // =====================================================================
    // EXTRA TEST A: Polynomial exactness
    //   Trap is exact for polynomials of degree <= 1
    //   2-pt Gauss is exact for polynomials of degree <= 3
    // =====================================================================
    std::cout << "\n===== EXTRA TEST A: Polynomial Exactness =====\n\n";
    {
        std::ofstream fout("data/test_polynomial_exactness.csv");
        fout << "degree,method,n,approx,exact,abs_err\n";
        fout << std::setprecision(17);

        struct PolyTest {
            int d;
            Func f;
            double exact;
        };
        // All on [0, 1]
        std::vector<PolyTest> ptests = {
            {0, [](double x){(void)x; return 1.0;}, 1.0},
            {1, [](double x){return x;}, 0.5},
            {2, [](double x){return x*x;}, 1.0/3.0},
            {3, [](double x){return x*x*x;}, 0.25},
            {4, [](double x){return x*x*x*x;}, 0.2},
        };
        std::vector<long> n_try = {1, 2, 4, 8, 16};
        for (auto& P : ptests) {
            for (long n : n_try) {
                auto rT = comp_trap_fixed(P.f, 0.0, 1.0, n);
                auto rG = comp_gauss2(P.f, 0.0, 1.0, n);
                fout << P.d << ",trap," << n << "," << rT.value << ","
                     << P.exact << "," << std::abs(rT.value - P.exact) << "\n";
                fout << P.d << ",gauss2," << n << "," << rG.value << ","
                     << P.exact << "," << std::abs(rG.value - P.exact) << "\n";
            }
        }
        std::cout << "  -> data/test_polynomial_exactness.csv\n";
    }

    // =====================================================================
    // EXTRA TEST B: Doubling identity for the adaptive rule
    //   Compare T_{2n} computed from scratch vs T_{2n} obtained by
    //   T_{2n} = 0.5 * T_n + H_new * sum_{i=0}^{n-1} f(a + (i+0.5)*H_old)
    // They should agree to machine precision.
    // =====================================================================
    std::cout << "\n===== EXTRA TEST B: Doubling Identity =====\n\n";
    {
        std::ofstream fout("data/test_doubling_identity.csv");
        fout << "integral,n,T2n_direct,T2n_recursion,diff\n";
        fout << std::setprecision(17);
        for (auto& T : tests) {
            for (int k = 0; k <= 8; ++k) {
                long n = (1L << k);
                auto rA = comp_trap_fixed(T.f, T.a, T.b, n);     // T_n
                auto rB = comp_trap_fixed(T.f, T.a, T.b, 2*n);   // T_{2n} direct
                // recursion: T_{2n} = 0.5*T_n + H_new*sum_mid
                double H_old = rA.H;
                double H_new = 0.5 * H_old;
                double sum_mid = 0.0;
                for (long i = 0; i < n; ++i) {
                    sum_mid += T.f(T.a + (i + 0.5) * H_old);
                }
                double T_recursion = 0.5 * rA.value + H_new * sum_mid;
                double diff = std::abs(T_recursion - rB.value);
                fout << T.name << "," << n << "," << rB.value << ","
                     << T_recursion << "," << diff << "\n";
            }
        }
        std::cout << "  -> data/test_doubling_identity.csv\n";
    }

    // =====================================================================
    // EXTRA TEST C: Empirical convergence rates
    //   For each integral and each method, compute
    //   p_hat = log2(E(n)/E(2n))  and verify p_hat -> 2 (Trap) or 4 (Gauss)
    // =====================================================================
    std::cout << "\n===== EXTRA TEST C: Empirical Convergence Rates =====\n\n";
    {
        std::ofstream fout("data/test_empirical_rates.csv");
        fout << "integral,method,n_small,n_big,err_small,err_big,p_observed\n";
        fout << std::setprecision(12);
        auto compute_rate = [](double e1, double e2) {
            if (e1 <= 0 || e2 <= 0) return std::numeric_limits<double>::quiet_NaN();
            return std::log2(e1 / e2);
        };
        // use doubling n sequence 16,32,64,128,256,512,1024
        std::vector<long> seq = {16, 32, 64, 128, 256, 512, 1024};
        for (auto& T : tests) {
            for (size_t i = 0; i + 1 < seq.size(); ++i) {
                long n1 = seq[i], n2 = seq[i+1];
                auto t1 = comp_trap_fixed(T.f, T.a, T.b, n1);
                auto t2 = comp_trap_fixed(T.f, T.a, T.b, n2);
                auto g1 = comp_gauss2(T.f, T.a, T.b, n1);
                auto g2 = comp_gauss2(T.f, T.a, T.b, n2);
                double eT1 = std::abs(t1.value - T.exact);
                double eT2 = std::abs(t2.value - T.exact);
                double eG1 = std::abs(g1.value - T.exact);
                double eG2 = std::abs(g2.value - T.exact);
                fout << T.name << ",trap," << n1 << "," << n2
                     << "," << eT1 << "," << eT2 << ","
                     << compute_rate(eT1, eT2) << "\n";
                fout << T.name << ",gauss2," << n1 << "," << n2
                     << "," << eG1 << "," << eG2 << ","
                     << compute_rate(eG1, eG2) << "\n";
            }
        }
        std::cout << "  -> data/test_empirical_rates.csv\n";
    }

    std::cout << "\nAll data written to ./data/*.csv\n";
    std::cout << "Files:\n";
    std::cout << "  data/integrals_info.csv\n";
    std::cout << "  data/task_ii_trap_fixed.csv\n";
    std::cout << "  data/task_iii_trap_adaptive.csv\n";
    std::cout << "  data/task_iv_gauss.csv\n";
    std::cout << "  data/trap_fixed_sweep.csv\n";
    std::cout << "  data/test_polynomial_exactness.csv\n";
    std::cout << "  data/test_doubling_identity.csv\n";
    std::cout << "  data/test_empirical_rates.csv\n";

    return 0;
}
