#include "cmath"
#include "vector"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++11-extensions"

using namespace std;

// NOT VERIFIED
pair<bool, double> IterativeMethod(double x0, double e, int N, double(*phi)(double)) {
    int k = 1;
    double x1;
    pair<bool, double> res;

    while (true) {
        x1 = phi(x0);

        if (abs(x1 - x0) < e) {
            res.first = true;
            res.second = x1;
            return res;
        }

        if (k == N) {
            res.first = false;
            return res;
        }

        k++;

        x0 = x1;
    }
}

typedef struct {
    int steps;
    bool ok;
    double res;
    vector<vector<double> > iter;
} AitkenResult;

AitkenResult AitkenAcceleration(double x0, double e, int N, double(*phi)(double)) {
    int k = 1;
    double x1, x2;
    AitkenResult result;
    vector<vector<double> > iter;
    while (true) {
        vector<double> step;
        x1 = phi(x0);
        step.push_back(x1);
        x2 = phi(x1);
        step.push_back(x2);
        x2 = x2 - pow((x2 - x1), 2) / (x2 - 2 * x1 + x0);
        step.push_back(x2);

        iter.push_back(step);

        if (abs(x2 - x0) < e) {
            return {k, true, x2, iter};
        }

        if (k == N) {
            return {k, false, 0.0, iter};
        }

        k++;
        x0 = x2;
    }
}

// NOT VERIFIED
pair<int, double> Newton(double x0, double e, int N, double(*f)(double), double(*ff)(double)) {
    int k = 1;
    double x1;
    pair<int, double> res;

    while (true) {
        if (ff(x0) == 0) {
            res.first = 2;
            return res;
        }

        x1 = x0 - f(x0) / ff(x0);

        if (abs(x1 - x0) < e) {
            res.first = 1;
            res.second = x1;
            return res;
        }

        if (k == N) {
            res.first = 0;
            return res;
        }

        k++;
        x0 = x1;
    }
}

// NOT VERIFIED
pair<bool, vector<double> > Jacobi(vector<vector<double> > a, vector<double> b, vector<double> x, double e, int N) {
    int k = 1;
    double y[1000];
    unsigned n = b.size();

    pair<bool, vector<double> > result;

    while (true) {
        for (int i = 1; i <= n; i++) {
            double sum_without_a_i = 0;
            for (int j = 1; j <= a[i].size(); j++) {
                if (j == i) continue;
                sum_without_a_i += a[i][j] * x[j];
            }

            y[i] = (b[i] - sum_without_a_i) / a[i][i];
        }

        double max = abs(x[1] - y[1]);
        for (int i = 2; i <= n; i++) {
            double val = abs(x[i] - y[i]);
            if (val > max) max = val;
        }

        if (max < e) {
            result.first = true;
            result.second = vector<double>(y, y + n);
            return result;
        }

        if (k == N) {
            result.first = false;
            return result;
        }

        k++;
        for (int i = 1; i <= n; i++) y[i] = x[i];
    }
}

// NOT VERIFIED
pair<bool, vector<double> >
GaussSeidel(vector<vector<double> > a, vector<double> b, vector<double> x0, double e, int N) {
    int k = 1;
    double y[1000], x[1000];
    unsigned n = b.size();

    pair<bool, vector<double> > result;

    for (int i = 1; i <= n; i++) {
        y[i] = x[i] = x0[i];
    }

    while (true) {
        for (int i = 1; i <= n; i++) {
            double sum_without_a_i = 0;
            for (int j = 1; j <= a[i].size(); j++) {
                if (j == i) continue;
                sum_without_a_i += a[i][j] * y[j];
            }

            y[i] = (b[i] - sum_without_a_i) / a[i][i];
        }

        double max = abs(x[1] - y[1]);
        for (int i = 2; i < n; i++) {
            double val = abs(x[i] - y[i]);
            if (val > max) max = val;
        }

        if (max < e) {
            result.first = true;
            result.second = vector<double>(y, y + n);
            return result;
        }

        if (k == N) {
            result.first = false;
            return result;
        }

        k++;
        for (int i = 1; i <= n; i++) y[i] = x[i];
    }
}

// NOT VERIFIED
vector<double> JordanElimination(vector<vector<double> > a, vector<double> b) {
    int k = 1;
    unsigned n = b.size();

    while (true) {
        for (int j = k + 1; j <= n; j++) {
            a[k][j] = a[k][j] / a[k][k];
        }

        b[k] = b[k] / a[k][k];

        int i = 1;

        while (true) {
            if (i != k) {
                for (int j = k + 1; j <= n; j++) {
                    a[i][j] = a[i][j] - a[i][k] * a[k][j];
                }
                b[i] = b[i] - a[i][k] * b[k];
            }
            if (i == n) break;
            i++;
        }

        if (k == n) {
            return b;
        }

        k++;
    }
}

vector<double> GaussElimination(vector<vector<double> > a, vector<double> b) {
    int k = 1;
    int n = static_cast<int>(b.size());

    int l, i, j;

    while (k < n + 1) {
        double d = a[k - 1][k - 1];
        l = k;
        i = k + 1;

        while (i < n + 1) {
            if (abs(a[i - 1][k - 1]) > abs(d)) {
                d = a[i - 1][k - 1];
                l = i;
            }
            i++;
        }

        if (d == 0) {
            return {};
        }

        if (l != k) {
            double t;
            for (j = k; j <= n; j++) {
                t = a[l - 1][j - 1];
                a[l - 1][j - 1] = a[k - 1][j - 1];
                a[k - 1][j - 1] = t;
            }
            t = b[k - 1];
            b[k - 1] = b[l - 1];
            b[l - 1] = t;
        }

        for (j = k + 1; j <= n; j++) {
            a[k - 1][j - 1] = a[k - 1][j - 1] / a[k - 1][k - 1];
        }
        b[k - 1] = b[k - 1] / a[k - 1][k - 1];

        for (i = k + 1; i <= n; i++) {
            for (j = k + 1; j <= n; j++) {
                a[i - 1][j - 1] = a[i - 1][j - 1] - a[i - 1][k - 1] * a[k - 1][j - 1];
            }
        }

        for (i = k + 1; i <= n; i++) {
            b[i - 1] = b[i - 1] - a[i - 1][k - 1] * b[k - 1];
        }

        k++;
    }

    for (i = n - 1; i >= 1; i--) {
        double sum = 0;
        for (j = i + 1; j <= n; j++) {
            sum += a[i - 1][j - 1] * b[j - 1];
        }
        b[i - 1] = b[i - 1] - sum;
    }

    return b;
}

// NOT VERIFIED
vector<double> Chasing(vector<double> a, vector<double> b, vector<double> c, vector<double> d) {
    c[1] = c[1] / b[1];
    d[1] = d[1] / b[1];
    int i = 2;
    unsigned n = b.size();

    double t;
    while (true) {
        t = b[i] - c[i - 1] * a[i];
        c[i] = c[i] / t;
        d[i] = (d[i] - d[i - 1] * a[i]) / t;

        if (i == n - 1) {
            d[n] = (d[n] - d[n - 1] * a[n]) / (b[n] - c[n - 1] * a[n]);
            for (int j = n - 1; j >= 1; j--) {
                d[i] = d[i] - c[i] * d[i + 1];
            }
            return d;
        }

        i++;
    }
}

#pragma clang diagnostic pop