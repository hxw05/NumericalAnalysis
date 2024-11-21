#include <iostream>
#include "Algorithms.h"

using namespace std;

double f1(double x) {
    return sqrt(x);
}

double f2(double x) {
    if (x < 1e-15) return 0;
    return sin(x) / x;
}

int main() {
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, 8) == 1507906) << endl;
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, 4) == 27026) << endl;
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, 3) == 5321) << endl;
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, -8) == 1106498) << endl;

    cout << (LagrangeInterpolation({{1,9}, {2,6}, {4,0}}, 5) == -3) << endl;
    cout << (LagrangeInterpolation({{1,9}, {2,6}, {4,0}}, 2) == 6) << endl;
    cout << (LagrangeInterpolation({{1,9}, {2,6}, {4,0}}, 8) == -12) << endl;
    cout << ((LagrangeInterpolation({{2, 3}, {7, 2}, {9, 12}}, -3) - 41.142857142857146) < 0.001) << endl;
    cout << (LagrangeInterpolation({{2, 3}, {7, 2}, {9, 12}}, 9) == 12) << endl;
    cout << ((LagrangeInterpolation({{2, 3}, {7, 2}, {9, 12}}, -10) - 156.94285714285715) < 0.001) << endl;

    cout << CompositeSimpson(0, 1, 4, f2) << endl;
}
