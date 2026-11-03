#include "bpt/matrix.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

namespace {

// A = [1 2 3]
//     [4 5 6]
bpt::Matrix make_a() {
    bpt::Matrix a(2, 3);
    a(0, 0) = 1.0;
    a(0, 1) = 2.0;
    a(0, 2) = 3.0;
    a(1, 0) = 4.0;
    a(1, 1) = 5.0;
    a(1, 2) = 6.0;
    return a;
}

}  // namespace

TEST(Matrix, ConstructorSetsSizeAndFill) {
    const bpt::Matrix m(2, 3, 1.5);
    EXPECT_EQ(m.rows(), 2);
    EXPECT_EQ(m.cols(), 3);
    EXPECT_FALSE(m.empty());
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_DOUBLE_EQ(m(r, c), 1.5);
        }
    }

    EXPECT_TRUE(bpt::Matrix().empty());
    EXPECT_TRUE(bpt::Matrix(0, 4).empty());
}

TEST(Matrix, StorageIsRowMajor) {
    const bpt::Matrix a = make_a();
    // Row 0 is stored first, so data() reads 1 2 3 4 5 6.
    const double expected[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(a.data()[i], expected[i]);
    }
}

TEST(Matrix, Identity) {
    const bpt::Matrix eye = bpt::Matrix::identity(3);
    ASSERT_EQ(eye.rows(), 3);
    ASSERT_EQ(eye.cols(), 3);
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_DOUBLE_EQ(eye(r, c), r == c ? 1.0 : 0.0);
        }
    }
}

TEST(Matrix, Transpose) {
    // A^T = [1 4]
    //       [2 5]
    //       [3 6]
    const bpt::Matrix at = make_a().transpose();
    ASSERT_EQ(at.rows(), 3);
    ASSERT_EQ(at.cols(), 2);
    EXPECT_DOUBLE_EQ(at(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(at(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(at(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(at(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(at(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(at(2, 1), 6.0);
}

TEST(Matrix, MatrixProduct) {
    // B = [ 7  8]
    //     [ 9 10]
    //     [11 12]
    bpt::Matrix b(3, 2);
    b(0, 0) = 7.0;
    b(0, 1) = 8.0;
    b(1, 0) = 9.0;
    b(1, 1) = 10.0;
    b(2, 0) = 11.0;
    b(2, 1) = 12.0;

    // A * B = [1*7 + 2*9 + 3*11   1*8 + 2*10 + 3*12]   [ 58  64]
    //         [4*7 + 5*9 + 6*11   4*8 + 5*10 + 6*12] = [139 154]
    const bpt::Matrix ab = make_a() * b;
    ASSERT_EQ(ab.rows(), 2);
    ASSERT_EQ(ab.cols(), 2);
    EXPECT_DOUBLE_EQ(ab(0, 0), 58.0);
    EXPECT_DOUBLE_EQ(ab(0, 1), 64.0);
    EXPECT_DOUBLE_EQ(ab(1, 0), 139.0);
    EXPECT_DOUBLE_EQ(ab(1, 1), 154.0);
}

TEST(Matrix, ProductWithZeroEntries) {
    // operator* skips zero entries of the left matrix, so check a case with
    // zeros on the diagonal.
    // [0 2] * [1 2] = [0*1 + 2*3   0*2 + 2*4] = [6 8]
    // [3 0]   [3 4]   [3*1 + 0*3   3*2 + 0*4]   [3 6]
    bpt::Matrix left(2, 2);
    left(0, 1) = 2.0;
    left(1, 0) = 3.0;
    bpt::Matrix right(2, 2);
    right(0, 0) = 1.0;
    right(0, 1) = 2.0;
    right(1, 0) = 3.0;
    right(1, 1) = 4.0;

    const bpt::Matrix out = left * right;
    EXPECT_DOUBLE_EQ(out(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(out(0, 1), 8.0);
    EXPECT_DOUBLE_EQ(out(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(out(1, 1), 6.0);
}

TEST(Matrix, IdentityTimesMatrixIsUnchanged) {
    const bpt::Matrix a = make_a();
    const bpt::Matrix out = bpt::Matrix::identity(2) * a;
    ASSERT_EQ(out.rows(), 2);
    ASSERT_EQ(out.cols(), 3);
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_DOUBLE_EQ(out(r, c), a(r, c));
        }
    }
}

TEST(Matrix, MultiplyVector) {
    // A * [1 0 -1] = [1 - 3, 4 - 6] = [-2, -2]
    const std::vector<double> y = make_a().multiply({1.0, 0.0, -1.0});
    ASSERT_EQ(y.size(), 2u);
    EXPECT_DOUBLE_EQ(y[0], -2.0);
    EXPECT_DOUBLE_EQ(y[1], -2.0);
}

TEST(Matrix, MultiplyTransposeVector) {
    // A^T * [1 2] = [1 + 8, 2 + 10, 3 + 12] = [9, 12, 15]
    const bpt::Matrix a = make_a();
    const std::vector<double> y = a.multiply_transpose({1.0, 2.0});
    ASSERT_EQ(y.size(), 3u);
    EXPECT_DOUBLE_EQ(y[0], 9.0);
    EXPECT_DOUBLE_EQ(y[1], 12.0);
    EXPECT_DOUBLE_EQ(y[2], 15.0);

    // Same answer as building the transpose first.
    const std::vector<double> check = a.transpose().multiply({1.0, 2.0});
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(y[i], check[i]);
    }
}

TEST(Matrix, AddScaledColumn) {
    // I + 0.5 * [3 4]^T in column 1:
    // [1 0 + 1.5]   [1 1.5]
    // [0 1 + 2.0] = [0 3.0]
    bpt::Matrix m = bpt::Matrix::identity(2);
    m.add_scaled_column(1, {3.0, 4.0}, 0.5);
    EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m(0, 1), 1.5);
    EXPECT_DOUBLE_EQ(m(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(m(1, 1), 3.0);
}

TEST(Matrix, ReflectUpperCopiesUpperTriangle) {
    // Upper triangle and diagonal are set, lower triangle holds -1 that
    // should be overwritten.
    // [ 1  2  3]      [1 2 3]
    // [-1  4  5]  ->  [2 4 5]
    // [-1 -1  6]      [3 5 6]
    bpt::Matrix m(3, 3, -1.0);
    m(0, 0) = 1.0;
    m(0, 1) = 2.0;
    m(0, 2) = 3.0;
    m(1, 1) = 4.0;
    m(1, 2) = 5.0;
    m(2, 2) = 6.0;
    m.reflect_upper();

    const double expected[3][3] = {{1.0, 2.0, 3.0}, {2.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_DOUBLE_EQ(m(r, c), expected[r][c]);
        }
    }
}

TEST(Matrix, DotAndInfNorm) {
    // 1*4 + 2*(-5) + 3*6 = 4 - 10 + 18 = 12
    EXPECT_DOUBLE_EQ(bpt::dot({1.0, 2.0, 3.0}, {4.0, -5.0, 6.0}), 12.0);
    EXPECT_DOUBLE_EQ(bpt::dot({}, {}), 0.0);

    // Largest absolute value, sign ignored.
    EXPECT_DOUBLE_EQ(bpt::inf_norm({1.0, -7.5, 3.0}), 7.5);
    EXPECT_DOUBLE_EQ(bpt::inf_norm({}), 0.0);
}
