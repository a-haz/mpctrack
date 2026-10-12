#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace bpt {

class Matrix {
public:
    Matrix() = default;
    Matrix(int rows, int cols, double fill = 0.0)
        : rows_(rows), cols_(cols), data_(static_cast<std::size_t>(rows * cols), fill) {}

    static Matrix identity(int n) {
        Matrix out(n, n, 0.0);
        for (int i = 0; i < n; ++i) {
            out(i, i) = 1.0;
        }
        return out;
    }

    int rows() const { return rows_; }
    int cols() const { return cols_; }
    bool empty() const { return rows_ == 0 || cols_ == 0; }

    double& operator()(int r, int c) {
        return data_[static_cast<std::size_t>(r * cols_ + c)];
    }
    double operator()(int r, int c) const {
        return data_[static_cast<std::size_t>(r * cols_ + c)];
    }

    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }

    Matrix transpose() const {
        Matrix out(cols_, rows_);
        for (int r = 0; r < rows_; ++r) {
            for (int c = 0; c < cols_; ++c) {
                out(c, r) = (*this)(r, c);
            }
        }
        return out;
    }

    Matrix operator*(const Matrix& rhs) const {
        assert(cols_ == rhs.rows_);
        Matrix out(rows_, rhs.cols_, 0.0);
        for (int r = 0; r < rows_; ++r) {
            for (int k = 0; k < cols_; ++k) {
                const double aik = (*this)(r, k);
                if (aik == 0.0) {
                    continue;
                }
                for (int c = 0; c < rhs.cols_; ++c) {
                    out(r, c) += aik * rhs(k, c);
                }
            }
        }
        return out;
    }

    std::vector<double> multiply(const std::vector<double>& x) const {
        assert(static_cast<int>(x.size()) == cols_);
        std::vector<double> y(static_cast<std::size_t>(rows_), 0.0);
        for (int r = 0; r < rows_; ++r) {
            double sum = 0.0;
            for (int c = 0; c < cols_; ++c) {
                sum += (*this)(r, c) * x[static_cast<std::size_t>(c)];
            }
            y[static_cast<std::size_t>(r)] = sum;
        }
        return y;
    }

    std::vector<double> multiply_transpose(const std::vector<double>& x) const {
        assert(static_cast<int>(x.size()) == rows_);
        std::vector<double> y(static_cast<std::size_t>(cols_), 0.0);
        for (int r = 0; r < rows_; ++r) {
            const double xr = x[static_cast<std::size_t>(r)];
            if (xr == 0.0) {
                continue;
            }
            for (int c = 0; c < cols_; ++c) {
                y[static_cast<std::size_t>(c)] += (*this)(r, c) * xr;
            }
        }
        return y;
    }

    void add_scaled_column(int col, const std::vector<double>& vec, double scale) {
        assert(static_cast<int>(vec.size()) == rows_);
        for (int r = 0; r < rows_; ++r) {
            (*this)(r, col) += scale * vec[static_cast<std::size_t>(r)];
        }
    }

    void reflect_upper() {
        assert(rows_ == cols_);
        for (int r = 0; r < rows_; ++r) {
            for (int c = r + 1; c < cols_; ++c) {
                (*this)(c, r) = (*this)(r, c);
            }
        }
    }

private:
    int rows_{0};
    int cols_{0};
    std::vector<double> data_;
};

inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

inline double inf_norm(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) {
        m = std::max(m, std::abs(x));
    }
    return m;
}

}  // namespace bpt
