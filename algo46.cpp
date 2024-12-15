#include "algo46.h"
#include "iostream"
#include "iomanip"

using namespace std;

double phi1(double x) {
    return pow(M_E, -x);
}

int main() {
    cout << "---------AitkenAcc x0=0.5, e=0.0001, N=100" << endl;
    auto aitken = AitkenAcceleration(0.5, 0.0001, 100, phi1);
    if (!aitken.ok) cout << "Aitken Failed " << aitken.steps << " steps" << endl;
    else cout << "Aitken OK after " << aitken.steps << " steps" << endl;
    cout << "x_bar \t x_tilde \t x_k" << endl;
    for (auto &it: aitken.iter) {
        cout << it[0] << '\t' << it[1] << '\t' << it[2] << endl;
    }

    cout << "---------Gauss Test 1" << endl;
    vector<vector<double> > a = {{2, -1, 3},
                                 {4, 2,  5},
                                 {1, 2,  0}};
    vector<double> b = {1,
                        4,
                        7};
    for (int i = 0; i < a.size(); i++) {
        cout << a[i][0] << "x1 + " << a[i][1] << "x2 + " << a[i][2] << "x3 = " << b[i] << endl;
    }
    auto gaussRes = GaussElimination(a, b);
    cout << "RESULT: " << endl;
    for (int i = 0; i < gaussRes.size(); i++) {
        cout << "x" << i + 1 << "=" << gaussRes[i] << " ";
    }

    cout << endl;
    cout << "---------Gauss Test 2" << endl;
    cout << fixed << setprecision(5);

    a = {{1e-5, 1},
         {1,    1}};
    b = {1, 2};
    for (int i = 0; i < a.size(); i++) {
        cout << a[i][0] << "x1 + " << a[i][1] << "x2 + " << a[i][2] << "x3 = " << b[i] << endl;
    }
    gaussRes = GaussElimination(a, b);
    cout << "RESULT: " << endl;
    for (int i = 0; i < gaussRes.size(); i++) {
        cout << "x" << i + 1 << "=" << gaussRes[i] << " ";
    }
}