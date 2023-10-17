#include "matrix.h"

#include <algorithm>
#include <cassert>

Matrix::Matrix() {

}

Matrix::Matrix(const std::vector<std::vector<Rational>> &value_) {
    Set(value_);
}

Matrix::Matrix(size_t m, size_t n)  {
    value.resize(m, std::vector<Rational>(n));
}

void Matrix::Set(const std::vector<std::vector<Rational>> &value_) {
    value = value_;
}

size_t Matrix::M() const {
    return value.size();
}

size_t Matrix::N() const {
    return value[0].size();
}

Matrix& operator+=(Matrix& a, const Matrix& b) {
    for (size_t i = 0; i < a.value.size(); ++i) {
        for (size_t j = 0; j < a.value[i].size(); ++j) {
            a.value[i][j] += b.value[i][j];
        }
    }
    return a;
}

Matrix operator+(const Matrix& a, const Matrix& b) {
    Matrix res = a;
    res += b;
    return res;
}

Matrix& operator*=(Matrix& a, const Matrix& b) {
    assert(a.N() == b.M());

    const size_t K = a.N();

    Matrix res(a.M(), b.N());

    for (size_t i = 0; i < a.M(); i++) {
        for (size_t j = 0; j < b.N(); j++) {
            for (size_t k = 0; k < K; k++) {
                res.value[i][j] += a.value[i][k] * b.value[k][j];
            }
        }
    }

    return a = res;
}

Matrix operator*(const Matrix& a, const Matrix& b) {
    Matrix res = a;
    res *= b;
    return res;
}

const std::vector <std::vector <Rational>>& Matrix::read() const {
    return value;
}

std::istream& operator>>(std::istream& cin, Matrix& a) {
    for (size_t i = 0; i < a.M(); i++) {
        for (size_t j = 0; j < a.N(); j++) {
            cin >> a.value[i][j];
        }
    }
    return cin;
}

std::ostream& operator<<(std::ostream& cout, const Matrix& a) {
    cout << "\\[\n"
            "\\left(\n"
            "\\begin{matrix}\n";
    for (size_t i = 0; i < a.M(); i++) {
        for (size_t j = 0; j < a.N(); j++) {
            cout << a.read()[i][j];
            if (j + 1 == a.N()) {
                cout << "\\\\\n";
            } else {
                cout << "& ";
            }
        }
    }
    cout << "\\end{matrix}\n"
            "\\right)\n"
            "\\]";
    cout << std::endl;
    return cout;
}

bool CmpRows(const std::vector <Rational>& row1, const std::vector <Rational>& row2) {
    assert(row1.size() == row2.size());
    const size_t N = row1.size();
    for (size_t i = 0; i < N; i++) {
        if (row1[i] == row2[i]) {
            continue;
        }

        if (row1[i].GetNumerator() == 0) {
            return false;
        }
        if (row2[i].GetNumerator() == 0) {
            return true;
        }

        int val1 = std::abs(row1[i].GetNumerator()) + row1[i].GetDenominator();
        int val2 = std::abs(row2[i].GetNumerator()) + row2[i].GetDenominator();

        if (val1 != val2) {
            return val1 < val2;
        }
        if (row1[i].GetNumerator() > 0 && row2[i].GetNumerator() < 0) {
            return true;
        }
        return false;
    }
    return false;
}

Matrix Matrix::GaussIteration() const{
    const auto& a = *this;
    Matrix res = a;
    std::sort(res.value.begin(), res.value.end(), CmpRows);

    bool changed = res != a;
    if (changed) {
        return res;
    }

    const size_t M = res.M();
    const size_t N = res.N();

    for (size_t row = 0; row < M; ++row) {
        size_t found_col = N;
        for (size_t col = 0; col < N; ++col) {
            if (res.value[row][col] != 0) {
                found_col = col;
                break;
            }
        }
        if (found_col == N) {
            return res;
        }

        auto val = res.value[row][found_col];
        if (val != 1) {
            for (size_t col = 0; col < N; ++col) {
                res.value[row][col] /= val;
            }
            return res;
        }

        size_t found_row = M;
        bool fixed = false;
        do {
            found_row = M;
            for (size_t search_row = 0; search_row < M; search_row++) {
                if (search_row != row && res.value[search_row][found_col] != 0) {
                    found_row = search_row;
                    break;
                }
            }

            if (found_row == M) {
                break;
            }

            fixed = true;
            const Rational factor = res.value[found_row][found_col];
            for (size_t col = 0; col < N; ++col) {
                res.value[found_row][col] -= factor * res.value[row][col];
            }
        } while (true);

        if (fixed) {
            return res;
        }
    }
    return res;
}

bool operator==(const Matrix& a, const Matrix& b) {
    const size_t M = a.M();
    const size_t N = a.N();
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            if (a.read()[i][j] != b.read()[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool operator!=(const Matrix& a, const Matrix& b) {
    return !(a == b);
}