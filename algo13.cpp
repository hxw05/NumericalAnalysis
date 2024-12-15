#include <iostream>
#include "algo13.h"

using namespace std;

double f2(double x) {
    double result;
    double x_temp = ((x == 0) ? 1e-15 : x);
    result = sin(x_temp) / x_temp;
    return result;
}

double f3(double x, double y) {
    return y - (2 * x) / y;
}

int main() {
    cout << "qinjiushao polynomial" << endl;

    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, 8) == 1507906) << endl;
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, 4) == 27026) << endl;
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, 3) == 5321) << endl;
    cout << (QinJiuShaoPolynomialVal({2, 0, 9, 8, -1, 6, 5}, -8) == 1106498) << endl;

    cout << "lagrange interpolation" << endl;

    cout << (LagrangeInterpolation({{1,9}, {2,6}, {4,0}}, 5) == -3) << endl;
    cout << (LagrangeInterpolation({{1,9}, {2,6}, {4,0}}, 2) == 6) << endl;
    cout << (LagrangeInterpolation({{1,9}, {2,6}, {4,0}}, 8) == -12) << endl;
    cout << (abs(LagrangeInterpolation({{2, 3}, {7, 2}, {9, 12}}, -3) - 41.142857142857146) < 0.001) << endl;
    cout << (LagrangeInterpolation({{2, 3}, {7, 2}, {9, 12}}, 9) == 12) << endl;
    cout << (abs(LagrangeInterpolation({{2, 3}, {7, 2}, {9, 12}}, -10) - 156.94285714285715) < 0.001) << endl;

    cout << "composite simpson" << endl;
    cout << (abs(CompositeSimpson(0, 1, 4, f2) - 0.9460869339517937) < 0.001) << endl;
    cout << (abs(CompositeSimpson(0, 2, 4, f2) - 1.6054969371344003) < 0.001) << endl;
    cout << (abs(CompositeSimpson(0, 9, 12, f2) - 1.6652433927251669) < 0.001) << endl;

    cout << "rectangle f2 a=0 b=1 step=10" << endl;
    vector<double> steps = FixedStepRectangle(0, 1, 10, f2);
    vector<double> stepRefs = {0.9397933, 0.9445135, 0.9456909, 0.9459850, 0.9460596, 0.9460769, 0.9460815, 0.9460827, 0.9460830, 0.9460831};
    for (int i = 0; i < steps.size(); i++) {
        cout << (abs(steps[i] - stepRefs[i]) < 0.001) << endl;
    }

    cout << "Romberg f2 a=0 b=1 e=1e-6" << endl;
    vector<double> Precise1E6Rectangle = VariableStepRectangle(0, 1, 0.0000001, f2);
    RombergResult Precise1E6Romberg = Romberg(0, 1, 0.0000001, f2);
    cout << "rectangle: " << Precise1E6Rectangle.back() << endl;
    cout << "romberg: " << Precise1E6Romberg.values.back() << endl;
    cout << "step_rectangle: " << Precise1E6Rectangle.size() << ", step_romberg: " << Precise1E6Romberg.values.size() << ", step_rectangle_binary: " << Precise1E6Romberg.step1 << ", error: " << Precise1E6Romberg.values.back() - Precise1E6Rectangle.back() << endl;

    cout << "ImprovedEuler: f3, y(0)=1, h=0.1, N=10" << endl;
    vector<pair<double, double>> ImprovedEulerR = ImprovedEuler(0, 1, 0.1, 10, f3);
    for (auto& p : ImprovedEulerR) {
        cout << "x: " << p.first << ", y: " << p.second << endl;
    }

    cout << "RungeKutta4th: f3, y(0)=1, h=0.2, N=5" << endl;
    vector<pair<double, double>> RK = RungeKutta4th(0, 1, 0.2, 5, f3);
    for (auto& p : RK) {
        cout << "x: " << p.first << ", y: " << p.second << endl;
    }

    cout << "AdamsPredCorrBaseRK4: f3, y(0)=1, h=0.1, N=10" << endl;
    AdamsResult AdamsPCbRK4 = AdamsPredictorCorrectorBasedOnRK4(0, 1, 0.1, 10, f3);
    cout << "==============predictors==============" << endl;
    for (auto& p : AdamsPCbRK4.predictors) {
        cout << "x: " << p.first << ", y: " << p.second << endl;
    }

    cout << "==============correctors==============" << endl;
    for (auto& p : AdamsPCbRK4.correctors) {
        cout << "x: " << p.first << ", y: " << p.second << endl;
    }
}