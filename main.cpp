#include "matrix.h"

#include <iostream>
#include <vector>
#include <string>

#define PROBLEM 1

void PrintGaussSolution(Matrix a) {
    Matrix changed = a;

    do {
        std::swap(changed, a);
        std::cout << a << std::endl;
        changed = a.GaussIteration();
    } while (changed != a);
}

int main() {
#if PROBLEM == 1
    std::freopen("input1.txt", "r", stdin);
    std::freopen("output1.txt", "w", stdout);

    Matrix a(5, 3), b(3, 4);
    std::cin >> a >> b;

    // std::cout << a << std::endl;
    // std::cout << b << std::endl;

    auto c = a * b;

    PrintGaussSolution(c);
#elif PROBLEM == 4
    std::freopen("input4.txt", "r", stdin);
    std::freopen("output4.txt", "w", stdout);
    Matrix a(9, 9);
    std::cin >> a;

    PrintGaussSolution(a);
#endif
	return 0;
}
