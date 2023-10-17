#include "rational.h"
#include <sstream>
#include <numeric>
#include <string>
#include <algorithm>

int64_t MultiplyAsInt64(int a, int b) {
    return static_cast<int64_t>(a) * static_cast<int64_t>(b);
}

Rational::Rational() {
    Set(0, 1);
}

Rational::Rational(int value) {  // NOLINT
    Set(value, 1);
}

Rational::Rational(int numer, int denom) {
    Set(numer, denom);
}

int Rational::GetNumerator() const {
    return numer_;
}

int Rational::GetDenominator() const {
    return denom_;
}

void Rational::SetNumerator(int value) {
    Set(value, denom_);
}

void Rational::SetDenominator(int value) {
    Set(numer_, value);
}

Rational& operator+=(Rational& lhs, const Rational& rhs) {
    int64_t numer = MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) +
        MultiplyAsInt64(lhs.GetDenominator(), rhs.GetNumerator());
    int64_t denom = MultiplyAsInt64(lhs.GetDenominator(), rhs.GetDenominator());
    lhs.Set(numer, denom);
    return lhs;
}

Rational& operator*=(Rational& lhs, const Rational& rhs) {
    int64_t numer = MultiplyAsInt64(lhs.GetNumerator(), rhs.GetNumerator());
    int64_t denom = MultiplyAsInt64(lhs.GetDenominator(), rhs.GetDenominator());
    lhs.Set(numer, denom);
    return lhs;
}

Rational& operator++(Rational& ratio) {
    int numer = ratio.GetNumerator() + ratio.GetDenominator();
    ratio.SetNumerator(numer);
    return ratio;
}

Rational& operator--(Rational& ratio) {
    int numer = ratio.GetNumerator() - ratio.GetDenominator();
    ratio.SetNumerator(numer);
    return ratio;
}

std::istream& operator>>(std::istream& is, Rational& ratio) {
    std::string text;
    is >> text;
    std::stringstream text_to_read;
    text_to_read << text;
    int numer = 0;
    int denom = 1;
    if (std::find(text.begin(), text.end(), '/') == text.end()) {
        text_to_read >> numer;
    }
    else {
        char slash = '/';
        text_to_read >> numer >> slash >> denom;
    }
    ratio.Set(numer, denom);
    return is;
}

Rational operator+(const Rational& ratio) {
    return ratio;
}

Rational operator-(const Rational& ratio) {
    return Rational(-ratio.GetNumerator(), ratio.GetDenominator());
}

Rational& operator-=(Rational& lhs, const Rational& rhs) {
    return lhs += -rhs;
}

Rational& operator/=(Rational& lhs, const Rational& rhs) {
    Rational divisor = Rational(rhs.GetDenominator(), rhs.GetNumerator());
    return lhs *= divisor;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
    Rational ratio = lhs;
    ratio += rhs;
    return ratio;
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    Rational ratio = lhs;
    ratio -= rhs;
    return ratio;
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    Rational ratio = lhs;
    ratio *= rhs;
    return ratio;
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    Rational ratio = lhs;
    ratio /= rhs;
    return ratio;
}

Rational operator++(Rational& ratio, int) {
    const Rational saved_ratio = ratio;
    ++ratio;
    return saved_ratio;
}

Rational operator--(Rational& ratio, int) {
    const Rational saved_ratio = ratio;
    --ratio;
    return saved_ratio;
}

bool operator<(const Rational& lhs, const Rational& rhs) {
    return MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) <
        MultiplyAsInt64(rhs.GetNumerator(), lhs.GetDenominator());
}

bool operator>(const Rational& lhs, const Rational& rhs) {
    return MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) >
        MultiplyAsInt64(rhs.GetNumerator(), lhs.GetDenominator());
}

bool operator<=(const Rational& lhs, const Rational& rhs) {
    return MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) <=
        MultiplyAsInt64(rhs.GetNumerator(), lhs.GetDenominator());
}

bool operator>=(const Rational& lhs, const Rational& rhs) {
    return MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) >=
        MultiplyAsInt64(rhs.GetNumerator(), lhs.GetDenominator());
}

bool operator==(const Rational& lhs, const Rational& rhs) {
    return MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) ==
        MultiplyAsInt64(rhs.GetNumerator(), lhs.GetDenominator());
}

bool operator!=(const Rational& lhs, const Rational& rhs) {
    return MultiplyAsInt64(lhs.GetNumerator(), rhs.GetDenominator()) !=
        MultiplyAsInt64(rhs.GetNumerator(), lhs.GetDenominator());
}

std::ostream& operator<<(std::ostream& os, const Rational& ratio) {
    if (ratio.GetDenominator() == 1) {
        return os << ratio.GetNumerator();
    }
    else {
        return os << ratio.GetNumerator() << "/" << ratio.GetDenominator();
    }
}

void Rational::Set(int64_t numer, int64_t denom) {
    if (denom == 0) {
        throw RationalDivisionByZero{};
    }

    if (denom < 0) {
        numer *= -1;
        denom *= -1;
    }

    int64_t gcd_of_numer_and_denom = std::gcd(numer, denom);

    numer /= gcd_of_numer_and_denom;
    denom /= gcd_of_numer_and_denom;

    numer_ = static_cast<int>(numer);
    denom_ = static_cast<int>(denom);
}
