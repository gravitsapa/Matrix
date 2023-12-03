#include "matrix.h"

#include <iostream>
#include <vector>
#include <string>
#include <cassert>

#define PROBLEM 2

int main() {
#if PROBLEM == 2
    std::freopen("input2.txt", "r", stdin);
    std::freopen("output2.txt", "w", stdout);

    Matrix a(4, 4), b(4, 4), c(4, 4), d(4, 4);
    std::cin >> a >> b >> c >> d;

    //std::cout << a << std::endl << b << std::endl << c << std::endl << d << std::endl;

    Matrix X = Subtract(FindInverse(Multiply(Multiply(FindInverse(a), d), FindInverse(c))), b);

    Solve(X);

    //check
    assert((a * FindInverse(X + b, 0) * c) == d);
#elif PROBLEM == 3
    std::freopen("input3.txt", "r", stdin);
    std::freopen("output3.txt", "w", stdout);

    Matrix a(4, 4);
    std::cin >> a;

    auto ch = characteristic(a);
    PrintChar(ch);

    std::cout << std::endl;

    std::cout << det(a - unit(4)) << ' ' << det(a + unit(4));

    std::cout << std::endl;

#endif
	return 0;
}
