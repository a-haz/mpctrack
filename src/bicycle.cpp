#include "bpt/bicycle.hpp"

#include <array>

namespace bpt {
namespace {

using Vec6 = std::array<double, 6>;
using Vec2 = std::array<double, 2>;
using Mat6 = std::array<std::array<double, 6>, 6>;
using Mat62 = std::array<std::array<double, 2>, 6>;

Mat6 zero6() {
    Mat6 A{};
    return A;
}

Mat62 zero62() {
    Mat62 B{};
    return B;
}

Mat6 eye6() {
    Mat6 A = zero6();
    for (int i = 0; i < 6; ++i) {
        A[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1.0;
    }
    return A;
}

Mat6 add6(const Mat6& A, const Mat6& B) {
    Mat6 C{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                B[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    return C;
}

Mat6 scale6(const Mat6& A, double s) {
    Mat6 C{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                s * A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    return C;
}

Mat6 mul6(const Mat6& A, const Mat6& B) {
    Mat6 C{};
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 6; ++k) {
            const double aik = A[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            if (aik == 0.0) {
                continue;
            }
            for (int j = 0; j < 6; ++j) {
                C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
                    aik * B[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
            }
        }
    }
    return C;
}

Mat62 mul62(const Mat6& A, const Mat62& B) {
    Mat62 C{};
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 6; ++k) {
            const double aik = A[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            if (aik == 0.0) {
                continue;
            }
            for (int j = 0; j < 2; ++j) {
                C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
                    aik * B[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
            }
        }
    }
    return C;
}

Mat62 add62(const Mat62& A, const Mat62& B) {
    Mat62 C{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 2; ++j) {
            C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                B[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    return C;
}

Mat62 scale62(const Mat62& A, double s) {
    Mat62 C{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 2; ++j) {
            C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                s * A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    return C;
}

Vec6 axpy(const Vec6& x, const Vec6& k, double scale) {
    Vec6 y{};
    for (int i = 0; i < 6; ++i) {
        y[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)] + scale * k[static_cast<std::size_t>(i)];
    }
    return y;
}

void dynamics(const Vec6& x, const Vec2& u, const VehicleParams& params, Vec6& dx, Mat6* A, Mat62* B) {
    const double cpsi = std::cos(x[2]);
    const double spsi = std::sin(x[2]);
    const double v = x[3];
    const double delta = x[4];
    const double accel = x[5];
    const double cdelta = std::cos(delta);
    const double tdelta = std::tan(delta);
    const double L = params.wheelbase;
    const double cdelta2 = cdelta * cdelta;

    dx[0] = v * cpsi;
    dx[1] = v * spsi;
    dx[2] = v * tdelta / L;
    dx[3] = accel;
    dx[4] = (u[0] - delta) / params.tau_delta;
    dx[5] = (u[1] - accel) / params.tau_accel;

    if (A != nullptr) {
        *A = zero6();
        (*A)[0][2] = -v * spsi;
        (*A)[0][3] = cpsi;
        (*A)[1][2] = v * cpsi;
        (*A)[1][3] = spsi;
        (*A)[2][3] = tdelta / L;
        (*A)[2][4] = v / (L * cdelta2);
        (*A)[3][5] = 1.0;
        (*A)[4][4] = -1.0 / params.tau_delta;
        (*A)[5][5] = -1.0 / params.tau_accel;
    }
    if (B != nullptr) {
        *B = zero62();
        (*B)[4][0] = 1.0 / params.tau_delta;
        (*B)[5][1] = 1.0 / params.tau_accel;
    }
}

Matrix to_matrix6(const Mat6& A) {
    Matrix M(6, 6, 0.0);
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 6; ++c) {
            M(r, c) = A[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
        }
    }
    return M;
}

Matrix to_matrix62(const Mat62& B) {
    Matrix M(6, 2, 0.0);
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 2; ++c) {
            M(r, c) = B[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
        }
    }
    return M;
}

}  // namespace

std::array<double, 6> continuous_derivative(const VehicleState& state, const Input& input,
                                                    const VehicleParams& params) {
    Vec6 dx{};
    dynamics(state.to_array(), input.to_array(), params, dx, nullptr, nullptr);
    return dx;
}

void continuous_jacobians(const VehicleState& state, const Input& input, const VehicleParams& params, Matrix& A,
                          Matrix& B) {
    Vec6 dx{};
    Mat6 Ac{};
    Mat62 Bc{};
    dynamics(state.to_array(), input.to_array(), params, dx, &Ac, &Bc);
    A = to_matrix6(Ac);
    B = to_matrix62(Bc);
}

VehicleState integrate_rk4(const VehicleState& state, const Input& input, const VehicleParams& params, double dt) {
    const Vec6 x = state.to_array();
    const Vec2 u = input.to_array();
    Vec6 k1{};
    Vec6 k2{};
    Vec6 k3{};
    Vec6 k4{};
    dynamics(x, u, params, k1, nullptr, nullptr);
    dynamics(axpy(x, k1, 0.5 * dt), u, params, k2, nullptr, nullptr);
    dynamics(axpy(x, k2, 0.5 * dt), u, params, k3, nullptr, nullptr);
    dynamics(axpy(x, k3, dt), u, params, k4, nullptr, nullptr);
    Vec6 xn{};
    for (int i = 0; i < 6; ++i) {
        xn[static_cast<std::size_t>(i)] =
            x[static_cast<std::size_t>(i)] +
            (dt / 6.0) * (k1[static_cast<std::size_t>(i)] + 2.0 * k2[static_cast<std::size_t>(i)] +
                          2.0 * k3[static_cast<std::size_t>(i)] + k4[static_cast<std::size_t>(i)]);
    }
    return VehicleState::from_array(xn);
}

void discrete_jacobians(const VehicleState& state, const Input& input, const VehicleParams& params, double dt,
                        Matrix& A, Matrix& B) {
    const Vec6 x = state.to_array();
    const Vec2 u = input.to_array();

    Vec6 k1{};
    Vec6 k2{};
    Vec6 k3{};
    Vec6 k4{};
    Mat6 A1{};
    Mat6 A2{};
    Mat6 A3{};
    Mat6 A4{};
    Mat62 B1{};
    Mat62 B2{};
    Mat62 B3{};
    Mat62 B4{};

    dynamics(x, u, params, k1, &A1, &B1);

    dynamics(axpy(x, k1, 0.5 * dt), u, params, k2, &A2, &B2);
    const Mat6 dx2 = add6(eye6(), scale6(A1, 0.5 * dt));
    const Mat6 Dk2 = mul6(A2, dx2);
    const Mat62 Du2 = add62(mul62(A2, scale62(B1, 0.5 * dt)), B2);

    dynamics(axpy(x, k2, 0.5 * dt), u, params, k3, &A3, &B3);
    const Mat6 dx3 = add6(eye6(), scale6(Dk2, 0.5 * dt));
    const Mat6 Dk3 = mul6(A3, dx3);
    const Mat62 Du3 = add62(mul62(A3, scale62(Du2, 0.5 * dt)), B3);

    dynamics(axpy(x, k3, dt), u, params, k4, &A4, &B4);
    const Mat6 dx4 = add6(eye6(), scale6(Dk3, dt));
    const Mat6 Dk4 = mul6(A4, dx4);
    const Mat62 Du4 = add62(mul62(A4, scale62(Du3, dt)), B4);

    Mat6 sumA = add6(A1, scale6(Dk2, 2.0));
    sumA = add6(sumA, scale6(Dk3, 2.0));
    sumA = add6(sumA, Dk4);
    const Mat6 Ad = add6(eye6(), scale6(sumA, dt / 6.0));

    Mat62 sumB = add62(B1, scale62(Du2, 2.0));
    sumB = add62(sumB, scale62(Du3, 2.0));
    sumB = add62(sumB, Du4);
    const Mat62 Bd = scale62(sumB, dt / 6.0);

    A = to_matrix6(Ad);
    B = to_matrix62(Bd);
}

}  // namespace bpt
