#pragma once

#include "rational.h"

#include <iostream>
#include <vector>

class Matrix {
public:
    Matrix();

    explicit Matrix(const std::vector <std::vector <Rational>>& value_);

    Matrix(size_t n, size_t m);

    size_t M() const;

    size_t N() const;

    friend Matrix& operator+=(Matrix& a, const Matrix& b);

    friend Matrix& operator*=(Matrix& a, const Matrix& b);

    const std::vector <std::vector <Rational>>& read() const;

    friend std::istream& operator>>(std::istream& cin, Matrix& a);

    Matrix GaussIteration() const;
private:
    std::vector <std::vector <Rational>> value;

    void Set(const std::vector <std::vector <Rational>>& value_);
};

std::ostream& operator<<(std::ostream& cout, const Matrix& a);

Matrix operator+(const Matrix& a, const Matrix& b);

Matrix operator*(const Matrix& a, const Matrix& b);

bool operator==(const Matrix& a, const Matrix& b);

bool operator!=(const Matrix& a, const Matrix& b);