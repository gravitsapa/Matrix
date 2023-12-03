#include "matrix.h"

#include <algorithm>
#include <cassert>
#include <vector>

#define DOLLAR "$"

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

Matrix::row::row(size_t index, Matrix* data) : index(index), data(data) {
}

Matrix::const_row::const_row(size_t index, const Matrix* data) : index(index), data(data) {
}

Rational &Matrix::row::operator[](size_t column) {
    return data->value[index][column];
}

const Rational &Matrix::const_row::operator[](size_t column) {
    return data->value[index][column];
}

Matrix::const_row Matrix::operator[](size_t index) const {
    return const_row(index, this);
}

Matrix::row Matrix::operator[](size_t index) {
    return row(index, this);
}

Matrix Matrix::operator-() const{
    Matrix res = *this;

    for (size_t i = 0; i < M(); i++) {
        for (size_t j = 0; j < N(); j++) {
            res[i][j] = -res[i][j];
        }
    }

    return res;
}

Matrix operator-(const Matrix& a, const Matrix& b) {
    Matrix res = a + (-b);
    return res;
}

Matrix& Matrix::operator+=(const Matrix& b) {
    for (size_t i = 0; i < value.size(); ++i) {
        for (size_t j = 0; j < value[i].size(); ++j) {
            value[i][j] += b.value[i][j];
        }
    }
    return *this;
}

Matrix operator+(const Matrix& a, const Matrix& b) {
    Matrix res = a;
    res += b;
    return res;
}

Matrix& Matrix::operator*=(const Matrix& b) {
    assert(N() == b.M());

    const size_t K = N();

    Matrix res(M(), b.N());

    for (size_t i = 0; i < M(); i++) {
        for (size_t j = 0; j < b.N(); j++) {
            for (size_t k = 0; k < K; k++) {
                res.value[i][j] += value[i][k] * b.value[k][j];
            }
        }
    }

    return *this = res;
}

Matrix operator*(const Matrix& a, const Matrix& b) {
    Matrix res = a;
    res *= b;
    return res;
}

std::istream& operator>>(std::istream& cin, Matrix& a) {
    for (size_t i = 0; i < a.M(); i++) {
        for (size_t j = 0; j < a.N(); j++) {
            cin >> a[i][j];
        }
    }
    return cin;
}

std::ostream& operator<<(std::ostream& cout, const Matrix& a) {
    cout << "\\left(\n"
            "\\begin{matrix}\n";
    for (size_t i = 0; i < a.M(); i++) {
        for (size_t j = 0; j < a.N(); j++) {
            cout << a[i][j];
            if (j + 1 == a.N()) {
                cout << "\\\\\n";
            } else {
                cout << "& ";
            }
        }
    }
    cout << "\\end{matrix}\n"
            "\\right)";
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
            if (a[i][j] != b[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool operator!=(const Matrix& a, const Matrix& b) {
    return !(a == b);
}

Matrix& Matrix::operator|=(const Matrix& b) {
    const size_t N_1 = N();
    const size_t N_2 = b.N();

    for (size_t i = 0; i < M(); i++) {
        value[i].resize(N_1 + N_2);
        for (size_t j = N_1; j < N_1 + N_2; j++) {
            value[i][j] = b.value[i][j - N_1];
        }
    }

    return *this;
}

Matrix operator|(const Matrix& a, const Matrix& b) {
    Matrix res = a;
    res |= b;
    return res;
}

Matrix unit(size_t n) {
    Matrix res(n, n);
    for (size_t i = 0; i < n; ++i) {
        res[i][i] = 1;
    }
    return res;
}

bool Matrix::is_square() const {
    return N() == M();
}

Matrix GaussSolution(Matrix a, bool print) {
    Matrix changed = a;

    do {
        std::swap(changed, a);
        if (print) {
            std::cout << "\\[\n" << a << "\\]\n" << std::endl;
        }
        changed = a.GaussIteration();
    } while (changed != a);

    return a;
}

Matrix FindInverse(Matrix a, bool print) {
    auto copy = a;
    if (print) {
        std::cout << "Находим обратную методом Гаусса:\n";
    }

    assert(a.is_square());
    a |= unit(a.N());

    auto got = GaussSolution(a, print);

    auto res = Matrix(a.M(), a.M());

    for (size_t i = 0; i < a.M(); i++) {
        for (size_t j = 0; j < a.M(); j++) {
            res[i][j] = got[i][j + a.M()];
        }
    }

    if (print) {
        std::cout << "Получаем обратную: \n" << DOLLAR << res << DOLLAR << " к " << DOLLAR << copy << DOLLAR << "\n\n";
    }

    return res;
}

Matrix Multiply(Matrix a, Matrix b) {
    std::cout << "Перемножим две матрицы: " << DOLLAR << a << DOLLAR << " и " << DOLLAR << b << " = " << DOLLAR << "\n\n";
    Matrix res = a * b;

    std::cout << DOLLAR <<  res << DOLLAR << "\n\n";

    return res;
}

Matrix Subtract(Matrix a, Matrix b) {
    std::cout << "Вычтем две матрицы: " << DOLLAR << a << DOLLAR << " и " << DOLLAR << b << " = ";
    Matrix res = a - b;

    std::cout << res << DOLLAR << "\n\n";

    return res;
}

Matrix Solve(Matrix a) {
    std::cout << "Ответ: " << DOLLAR << a << DOLLAR << "\n\n";
    return a;
}

Rational det(Matrix a) {
    if (a.M() == 0) return 1;

    std::vector <int> sigma(a.N());
    for(int i = 0; i < a.N(); ++i) sigma[i] = i;

    Rational ans = 0;
    do {
        int inversions = 0;

        for (int i = 0; i < a.N(); ++i) {
            for (int j = i + 1; j < a.N(); ++j) {
                inversions += sigma[i] > sigma[j];
            }
        }

        Rational mult = 1;
        if (inversions % 2) mult = -1;

        for (int i = 0; i < a.N(); ++i) {
            mult *= a[i][sigma[i]];
        }

        ans += mult;
    } while (next_permutation(sigma.begin(), sigma.end()));

    return ans;
}

void all_subsets(int i, int n, int k, std::vector <std::vector <int>>& ans, std::vector <int>& now) {
    if (i == n) {
        assert(now.size() == k);
        ans.push_back(now);
        return;
    }

    if (now.size() < k) {
        now.push_back(i);
        all_subsets(i + 1, n, k, ans, now);
        now.pop_back();
    }
    if (n - i - 1 >= k - now.size()) {
        all_subsets(i + 1, n, k, ans, now);
    }
}

std::vector <std::vector <int>> all_subsets(int n, int k) {
    std::vector <std::vector <int>> ans;
    std::vector <int> now;
    all_subsets(0, n, k, ans, now);
    return ans;
}

std::vector <Rational> characteristic(Matrix a) {
    std::vector <Rational> ans(a.N() + 1);

    for (int i = 0; i <= a.N(); ++i) {
        int k = a.N() - i;

        auto subsets = all_subsets(a.N(), k);

        Rational sum = 0;
        for (auto subset : subsets) {
            Matrix chosen(k, k);
            for (int x = 0; x < k; x++) {
                for (int y = 0; y < k; ++y) {
                    chosen[x][y] = a[subset[x]][subset[y]];
                }
            }

            sum += det(chosen);
        }

        ans[i] = (((a.N() - i) % 2 == 0) ? sum : -sum);
    }

    return ans;
}

void PrintChar(std::vector <Rational> ch) {
    std::cout << DOLLAR;
    for (int i = (int)ch.size() - 1; i >= 0; --i) {
        std::cout << " ";
        if (i < (int)ch.size() - 1 && ch[i] >= 0) std::cout << "+ ";
        std::cout << ch[i];
        if (i) std::cout << "\\lambda^{" << i << "}";
    }
    std:: cout << DOLLAR;
}