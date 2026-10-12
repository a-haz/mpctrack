#include "bpt/qp_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace bpt {
namespace {

bool cholesky_factor(const Matrix& A, Matrix& L, double min_pivot) {
    const int n = A.rows();
    L = Matrix(n, n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = A(i, j);
            for (int k = 0; k < j; ++k) {
                sum -= L(i, k) * L(j, k);
            }
            if (i == j) {
                if (sum <= min_pivot) {
                    return false;
                }
                L(i, i) = std::sqrt(sum);
            } else {
                L(i, j) = sum / L(j, j);
            }
        }
    }
    return true;
}

std::vector<double> chol_solve(const Matrix& L, const std::vector<double>& b) {
    const int n = L.rows();
    std::vector<double> y(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int k = 0; k < i; ++k) {
            sum -= L(i, k) * y[static_cast<std::size_t>(k)];
        }
        y[static_cast<std::size_t>(i)] = sum / L(i, i);
    }
    std::vector<double> x(static_cast<std::size_t>(n));
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[static_cast<std::size_t>(i)];
        for (int k = i + 1; k < n; ++k) {
            sum -= L(k, i) * x[static_cast<std::size_t>(k)];
        }
        x[static_cast<std::size_t>(i)] = sum / L(i, i);
    }
    return x;
}

bool factor_spd(Matrix A, Matrix& L) {
    const int n = A.rows();
    double jitter = 0.0;
    for (int attempt = 0; attempt < 6; ++attempt) {
        if (jitter > 0.0) {
            for (int i = 0; i < n; ++i) {
                A(i, i) += jitter;
            }
        }
        if (cholesky_factor(A, L, 1e-18)) {
            return true;
        }
        jitter = (jitter == 0.0) ? 1e-12 : jitter * 10.0;
    }
    return false;
}

bool lu_solve(Matrix A, std::vector<double> b, std::vector<double>& x) {
    const int n = A.rows();
    if (n == 0) {
        x.clear();
        return true;
    }
    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double best = std::abs(A(k, k));
        for (int i = k + 1; i < n; ++i) {
            const double value = std::abs(A(i, k));
            if (value > best) {
                best = value;
                pivot = i;
            }
        }
        if (best < 1e-14) {
            return false;
        }
        if (pivot != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A(k, j), A(pivot, j));
            }
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(pivot)]);
        }
        for (int i = k + 1; i < n; ++i) {
            A(i, k) /= A(k, k);
            const double lik = A(i, k);
            for (int j = k + 1; j < n; ++j) {
                A(i, j) -= lik * A(k, j);
            }
        }
    }
    std::vector<double> y(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int j = 0; j < i; ++j) {
            sum -= A(i, j) * y[static_cast<std::size_t>(j)];
        }
        y[static_cast<std::size_t>(i)] = sum;
    }
    x.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j) {
            sum -= A(i, j) * x[static_cast<std::size_t>(j)];
        }
        x[static_cast<std::size_t>(i)] = sum / A(i, i);
    }
    return true;
}

double step_length(const std::vector<double>& v, const std::vector<double>& dv) {
    double alpha = 1.0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (dv[i] < -1e-16) {
            alpha = std::min(alpha, -v[i] / dv[i]);
        }
    }
    return std::max(0.0, alpha);
}

QpSolution solve_unconstrained(const QpProblem& problem) {
    QpSolution solution;
    Matrix L;
    if (!factor_spd(problem.P, L)) {
        return solution;
    }
    std::vector<double> rhs(problem.q.size());
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        rhs[i] = -problem.q[i];
    }
    solution.x = chol_solve(L, rhs);
    solution.success = true;
    solution.iterations = 0;
    return solution;
}

QpSolution solve_equalities(const QpProblem& problem) {
    QpSolution solution;
    const int n = problem.P.rows();
    const int ne = problem.A.rows();
    Matrix K(n + ne, n + ne, 0.0);
    std::vector<double> rhs(static_cast<std::size_t>(n + ne), 0.0);
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            K(r, c) = problem.P(r, c);
        }
        rhs[static_cast<std::size_t>(r)] = -problem.q[static_cast<std::size_t>(r)];
    }
    for (int r = 0; r < ne; ++r) {
        for (int c = 0; c < n; ++c) {
            K(n + r, c) = problem.A(r, c);
            K(c, n + r) = problem.A(r, c);
        }
        rhs[static_cast<std::size_t>(n + r)] = problem.b[static_cast<std::size_t>(r)];
    }
    std::vector<double> full;
    if (!lu_solve(K, rhs, full)) {
        return solution;
    }
    solution.x.assign(full.begin(), full.begin() + n);
    solution.success = true;
    return solution;
}

QpSolution solve_inequalities(const QpProblem& problem, const QpSettings& settings,
                              const std::vector<double>* warm_start) {
    QpSolution solution;
    const int n = problem.P.rows();
    const int m = problem.G.rows();
    Matrix Lp;
    if (!factor_spd(problem.P, Lp)) {
        return solution;
    }

    // W = P^{-1} G^T, GPG = G P^{-1} G^T. Both are constant across IPM iterations.
    Matrix W(n, m, 0.0);
    for (int j = 0; j < m; ++j) {
        std::vector<double> column(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            column[static_cast<std::size_t>(i)] = problem.G(j, i);
        }
        const std::vector<double> solved = chol_solve(Lp, column);
        for (int i = 0; i < n; ++i) {
            W(i, j) = solved[static_cast<std::size_t>(i)];
        }
    }
    const Matrix GPG = problem.G * W;

    std::vector<double> x(static_cast<std::size_t>(n), 0.0);
    if (warm_start != nullptr && static_cast<int>(warm_start->size()) == n) {
        x = *warm_start;
    }
    std::vector<double> z(static_cast<std::size_t>(m), 1.0);
    std::vector<double> s(static_cast<std::size_t>(m), 1.0);
    if (warm_start != nullptr) {
        const std::vector<double> slack = problem.h;
        const std::vector<double> Gx = problem.G.multiply(x);
        for (int i = 0; i < m; ++i) {
            s[static_cast<std::size_t>(i)] =
                std::max(1e-4, slack[static_cast<std::size_t>(i)] - Gx[static_cast<std::size_t>(i)]);
        }
    }

    const double tol = settings.tolerance;
    const double tau = settings.fraction_to_boundary;
    std::vector<double> best_x = x;
    double best_score = std::numeric_limits<double>::infinity();

    for (int iteration = 0; iteration < settings.max_iterations; ++iteration) {
        const std::vector<double> Px = problem.P.multiply(x);
        const std::vector<double> Gtz = problem.G.multiply_transpose(z);
        const std::vector<double> Gx = problem.G.multiply(x);
        std::vector<double> rx(static_cast<std::size_t>(n));
        std::vector<double> rz(static_cast<std::size_t>(m));
        double mu = 0.0;
        for (int i = 0; i < n; ++i) {
            rx[static_cast<std::size_t>(i)] =
                Px[static_cast<std::size_t>(i)] + problem.q[static_cast<std::size_t>(i)] +
                Gtz[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < m; ++i) {
            rz[static_cast<std::size_t>(i)] =
                Gx[static_cast<std::size_t>(i)] + s[static_cast<std::size_t>(i)] -
                problem.h[static_cast<std::size_t>(i)];
            mu += s[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
        }
        mu /= static_cast<double>(m);
        const double dual_inf = inf_norm(rx);
        const double primal_inf = inf_norm(rz);
        const double score = std::max(dual_inf, std::max(primal_inf, mu));
        if (score < best_score) {
            best_score = score;
            best_x = x;
            solution.dual_residual = dual_inf;
            solution.primal_residual = primal_inf;
            solution.complementarity = mu;
            solution.iterations = iteration;
        }
        if (dual_inf < tol && primal_inf < tol && mu < tol) {
            solution.x = x;
            solution.success = true;
            return solution;
        }

        std::vector<double> D(static_cast<std::size_t>(m));
        for (int i = 0; i < m; ++i) {
            D[static_cast<std::size_t>(i)] =
                s[static_cast<std::size_t>(i)] / std::max(z[static_cast<std::size_t>(i)], 1e-16);
        }
        Matrix Schur = GPG;
        for (int i = 0; i < m; ++i) {
            Schur(i, i) += D[static_cast<std::size_t>(i)];
        }
        Schur.reflect_upper();
        Matrix Ls;
        if (!factor_spd(Schur, Ls)) {
            break;
        }

        auto newton = [&](const std::vector<double>& rt, std::vector<double>& dx, std::vector<double>& dz,
                          std::vector<double>& ds) {
            std::vector<double> rhs1(static_cast<std::size_t>(n));
            std::vector<double> rhs2(static_cast<std::size_t>(m));
            for (int i = 0; i < n; ++i) {
                rhs1[static_cast<std::size_t>(i)] = -rx[static_cast<std::size_t>(i)];
            }
            for (int i = 0; i < m; ++i) {
                const double zi = std::max(z[static_cast<std::size_t>(i)], 1e-16);
                rhs2[static_cast<std::size_t>(i)] =
                    -rz[static_cast<std::size_t>(i)] + rt[static_cast<std::size_t>(i)] / zi;
            }
            const std::vector<double> p_rhs1 = chol_solve(Lp, rhs1);
            std::vector<double> schur_rhs = problem.G.multiply(p_rhs1);
            for (int i = 0; i < m; ++i) {
                schur_rhs[static_cast<std::size_t>(i)] -= rhs2[static_cast<std::size_t>(i)];
            }
            dz = chol_solve(Ls, schur_rhs);
            const std::vector<double> Gtdz = problem.G.multiply_transpose(dz);
            std::vector<double> dx_rhs(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                dx_rhs[static_cast<std::size_t>(i)] =
                    rhs1[static_cast<std::size_t>(i)] - Gtdz[static_cast<std::size_t>(i)];
            }
            dx = chol_solve(Lp, dx_rhs);
            ds.assign(static_cast<std::size_t>(m), 0.0);
            for (int i = 0; i < m; ++i) {
                const double zi = std::max(z[static_cast<std::size_t>(i)], 1e-16);
                ds[static_cast<std::size_t>(i)] =
                    (-rt[static_cast<std::size_t>(i)] - s[static_cast<std::size_t>(i)] * dz[static_cast<std::size_t>(i)]) /
                    zi;
            }
        };

        std::vector<double> rt(static_cast<std::size_t>(m));
        for (int i = 0; i < m; ++i) {
            rt[static_cast<std::size_t>(i)] = s[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
        }
        std::vector<double> dx;
        std::vector<double> dz;
        std::vector<double> ds;
        newton(rt, dx, dz, ds);

        const double alpha_p_aff = step_length(s, ds);
        const double alpha_d_aff = step_length(z, dz);
        double mu_aff = 0.0;
        for (int i = 0; i < m; ++i) {
            const double s_aff = s[static_cast<std::size_t>(i)] + alpha_p_aff * ds[static_cast<std::size_t>(i)];
            const double z_aff = z[static_cast<std::size_t>(i)] + alpha_d_aff * dz[static_cast<std::size_t>(i)];
            mu_aff += s_aff * z_aff;
        }
        mu_aff /= static_cast<double>(m);
        double sigma = 0.0;
        if (mu > 0.0) {
            const double ratio = mu_aff / mu;
            sigma = std::min(1.0, ratio * ratio * ratio);
            sigma = std::max(0.0, sigma);
        }
        for (int i = 0; i < m; ++i) {
            rt[static_cast<std::size_t>(i)] =
                s[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)] +
                ds[static_cast<std::size_t>(i)] * dz[static_cast<std::size_t>(i)] - sigma * mu;
        }
        newton(rt, dx, dz, ds);

        const double alpha_p = std::min(1.0, tau * step_length(s, ds));
        const double alpha_d = std::min(1.0, tau * step_length(z, dz));
        if (alpha_p == 0.0 && alpha_d == 0.0) {
            break;
        }
        for (int i = 0; i < n; ++i) {
            x[static_cast<std::size_t>(i)] += alpha_p * dx[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < m; ++i) {
            s[static_cast<std::size_t>(i)] =
                std::max(1e-16, s[static_cast<std::size_t>(i)] + alpha_p * ds[static_cast<std::size_t>(i)]);
            z[static_cast<std::size_t>(i)] =
                std::max(1e-16, z[static_cast<std::size_t>(i)] + alpha_d * dz[static_cast<std::size_t>(i)]);
        }
    }

    solution.x = best_x;
    solution.success = best_score < std::max(1e-5, settings.tolerance * 100.0);
    return solution;
}

}  // namespace

QpSolution solve_qp(const QpProblem& problem, const QpSettings& settings,
                    const std::vector<double>* warm_start) {
    const int n = problem.P.rows();
    if (n <= 0 || static_cast<int>(problem.q.size()) != n || problem.P.cols() != n) {
        return {};
    }
    const int m = problem.G.rows();
    const int ne = problem.A.rows();
    if (m > 0 && (problem.G.cols() != n || static_cast<int>(problem.h.size()) != m)) {
        return {};
    }
    if (ne > 0 && (problem.A.cols() != n || static_cast<int>(problem.b.size()) != ne)) {
        return {};
    }
    if (m == 0 && ne == 0) {
        return solve_unconstrained(problem);
    }
    if (m == 0) {
        return solve_equalities(problem);
    }
    if (ne == 0) {
        return solve_inequalities(problem, settings, warm_start);
    }

    // General equalities plus inequalities: full KKT Newton system. Used by small
    // problems; the path tracker condenses dynamics and keeps only inequalities.
    QpSolution solution;
    std::vector<double> x(static_cast<std::size_t>(n), 0.0);
    std::vector<double> y(static_cast<std::size_t>(ne), 0.0);
    std::vector<double> z(static_cast<std::size_t>(m), 1.0);
    std::vector<double> s(static_cast<std::size_t>(m), 1.0);
    const double tol = settings.tolerance;
    const double tau = settings.fraction_to_boundary;

    for (int iteration = 0; iteration < settings.max_iterations; ++iteration) {
        const std::vector<double> Px = problem.P.multiply(x);
        const std::vector<double> Gtz = problem.G.multiply_transpose(z);
        const std::vector<double> Aty = problem.A.multiply_transpose(y);
        const std::vector<double> Gx = problem.G.multiply(x);
        const std::vector<double> Ax = problem.A.multiply(x);
        std::vector<double> rx(static_cast<std::size_t>(n));
        std::vector<double> re(static_cast<std::size_t>(ne));
        std::vector<double> rz(static_cast<std::size_t>(m));
        double mu = 0.0;
        for (int i = 0; i < n; ++i) {
            rx[static_cast<std::size_t>(i)] = Px[static_cast<std::size_t>(i)] +
                                              problem.q[static_cast<std::size_t>(i)] +
                                              Gtz[static_cast<std::size_t>(i)] +
                                              Aty[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < ne; ++i) {
            re[static_cast<std::size_t>(i)] = Ax[static_cast<std::size_t>(i)] - problem.b[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < m; ++i) {
            rz[static_cast<std::size_t>(i)] = Gx[static_cast<std::size_t>(i)] + s[static_cast<std::size_t>(i)] -
                                              problem.h[static_cast<std::size_t>(i)];
            mu += s[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
        }
        mu /= static_cast<double>(m);
        const double dual_inf = inf_norm(rx);
        const double eq_inf = inf_norm(re);
        const double primal_inf = std::max(inf_norm(rz), eq_inf);
        solution.iterations = iteration;
        solution.dual_residual = dual_inf;
        solution.primal_residual = primal_inf;
        solution.complementarity = mu;
        if (dual_inf < tol && primal_inf < tol && mu < tol) {
            solution.x = x;
            solution.success = true;
            return solution;
        }

        std::vector<double> D(static_cast<std::size_t>(m));
        for (int i = 0; i < m; ++i) {
            D[static_cast<std::size_t>(i)] =
                s[static_cast<std::size_t>(i)] / std::max(z[static_cast<std::size_t>(i)], 1e-16);
        }

        auto newton = [&](const std::vector<double>& rt, std::vector<double>& dx, std::vector<double>& dy,
                          std::vector<double>& dz, std::vector<double>& ds) {
            const int N = n + ne + m;
            Matrix K(N, N, 0.0);
            std::vector<double> rhs(static_cast<std::size_t>(N));
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    K(r, c) = problem.P(r, c);
                }
                rhs[static_cast<std::size_t>(r)] = -rx[static_cast<std::size_t>(r)];
            }
            for (int r = 0; r < ne; ++r) {
                for (int c = 0; c < n; ++c) {
                    K(n + r, c) = problem.A(r, c);
                    K(c, n + r) = problem.A(r, c);
                }
                rhs[static_cast<std::size_t>(n + r)] = -re[static_cast<std::size_t>(r)];
            }
            for (int r = 0; r < m; ++r) {
                for (int c = 0; c < n; ++c) {
                    K(n + ne + r, c) = problem.G(r, c);
                    K(c, n + ne + r) = problem.G(r, c);
                }
                K(n + ne + r, n + ne + r) = -D[static_cast<std::size_t>(r)];
                const double zi = std::max(z[static_cast<std::size_t>(r)], 1e-16);
                rhs[static_cast<std::size_t>(n + ne + r)] =
                    -rz[static_cast<std::size_t>(r)] + rt[static_cast<std::size_t>(r)] / zi;
            }
            std::vector<double> sol;
            if (!lu_solve(K, rhs, sol)) {
                for (int i = 0; i < N; ++i) {
                    K(i, i) += (i < n) ? 1e-10 : -1e-10;
                }
                if (!lu_solve(K, rhs, sol)) {
                    dx.assign(static_cast<std::size_t>(n), 0.0);
                    dy.assign(static_cast<std::size_t>(ne), 0.0);
                    dz.assign(static_cast<std::size_t>(m), 0.0);
                    ds.assign(static_cast<std::size_t>(m), 0.0);
                    return;
                }
            }
            dx.assign(sol.begin(), sol.begin() + n);
            dy.assign(sol.begin() + n, sol.begin() + n + ne);
            dz.assign(sol.begin() + n + ne, sol.end());
            ds.assign(static_cast<std::size_t>(m), 0.0);
            for (int i = 0; i < m; ++i) {
                const double zi = std::max(z[static_cast<std::size_t>(i)], 1e-16);
                ds[static_cast<std::size_t>(i)] =
                    (-rt[static_cast<std::size_t>(i)] -
                     s[static_cast<std::size_t>(i)] * dz[static_cast<std::size_t>(i)]) /
                    zi;
            }
        };

        std::vector<double> rt(static_cast<std::size_t>(m));
        for (int i = 0; i < m; ++i) {
            rt[static_cast<std::size_t>(i)] = s[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
        }
        std::vector<double> dx, dy, dz, ds;
        newton(rt, dx, dy, dz, ds);
        const double ap_aff = step_length(s, ds);
        const double ad_aff = step_length(z, dz);
        double mu_aff = 0.0;
        for (int i = 0; i < m; ++i) {
            mu_aff += (s[static_cast<std::size_t>(i)] + ap_aff * ds[static_cast<std::size_t>(i)]) *
                      (z[static_cast<std::size_t>(i)] + ad_aff * dz[static_cast<std::size_t>(i)]);
        }
        mu_aff /= static_cast<double>(m);
        double sigma = 0.0;
        if (mu > 0.0) {
            const double ratio = mu_aff / mu;
            sigma = std::clamp(ratio * ratio * ratio, 0.0, 1.0);
        }
        for (int i = 0; i < m; ++i) {
            rt[static_cast<std::size_t>(i)] = s[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)] +
                                              ds[static_cast<std::size_t>(i)] * dz[static_cast<std::size_t>(i)] -
                                              sigma * mu;
        }
        newton(rt, dx, dy, dz, ds);
        const double ap = std::min(1.0, tau * step_length(s, ds));
        const double ad = std::min(1.0, tau * step_length(z, dz));
        for (int i = 0; i < n; ++i) {
            x[static_cast<std::size_t>(i)] += ap * dx[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < ne; ++i) {
            y[static_cast<std::size_t>(i)] += ad * dy[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < m; ++i) {
            s[static_cast<std::size_t>(i)] =
                std::max(1e-16, s[static_cast<std::size_t>(i)] + ap * ds[static_cast<std::size_t>(i)]);
            z[static_cast<std::size_t>(i)] =
                std::max(1e-16, z[static_cast<std::size_t>(i)] + ad * dz[static_cast<std::size_t>(i)]);
        }
    }
    solution.x = x;
    return solution;
}

}  // namespace bpt
