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

    class row {
    public:
        row(size_t index, Matrix* data);

        Rational& operator[](size_t column);

    private:
        size_t index;
        Matrix* data;
    };

    class const_row {
    public:
        const_row(size_t index, const Matrix* data);

        const Rational& operator[](size_t column);

    private:
        size_t index;
        const Matrix* data;
    };

    row operator[](size_t index);

    const_row operator[](size_t index) const;

    Matrix& operator+=(const Matrix& b);

    Matrix& operator*=(const Matrix& b);

    Matrix& operator|=(const Matrix& b);

    Matrix operator-() const;

    Matrix GaussIteration() const;

    bool is_square() const;
private:
    std::vector <std::vector <Rational>> value;

    void Set(const std::vector <std::vector <Rational>>& value_);
};

Matrix unit(size_t n);

std::istream& operator>>(std::istream& cin, Matrix& a);

std::ostream& operator<<(std::ostream& cout, const Matrix& a);

Matrix operator+(const Matrix& a, const Matrix& b);

Matrix operator-(const Matrix& a, const Matrix& b);

Matrix operator*(const Matrix& a, const Matrix& b);

bool operator==(const Matrix& a, const Matrix& b);

bool operator!=(const Matrix& a, const Matrix& b);

Matrix operator|(const Matrix& a, const Matrix& b);

Matrix FindInverse(Matrix a, bool print = true);

Matrix GaussSolution(Matrix a, bool print = true);

Matrix Multiply(Matrix a, Matrix b);

Matrix Subtract(Matrix a, Matrix b);

Matrix Solve(Matrix a);

Rational det(Matrix a);

std::vector <Rational> characteristic(Matrix a);

void PrintChar(std::vector <Rational> ch);