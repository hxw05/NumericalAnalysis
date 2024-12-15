#include "vector"
#include "cmath"
#include "iostream"

using namespace std;

int QinJiuShaoPolynomialVal(const std::vector<int> &a, int x) {
    int v = a[a.size() - 1];
    int k = 1;
    while (k < a.size()) {
        v = x * v + a[a.size() - 1 - k];
        k++;
    }
    return v;
}

double LagrangeInterpolation(const std::vector<std::pair<int, int>> &data, int x) {
    double y = 0;
    int k = 0;
    while (k < data.size()) {
        double t = 1;
        for (int j = 0; j < data.size(); j++) {
            if (j == k) continue;
            t = t * (x - data[j].first) / (data[k].first - data[j].first);
        }
        y += t * data[k].second;
        k++;
    }
    return y;
}

double CompositeSimpson(double a, double b, int n, double (*f)(double)) {
    double h = (b - a) / n;
    double s = 0;
    double x = a;

    for (int k = 0; k <= n - 1; k++) {
        s += f(x);
        x += h / 2;
        s += 4 * f(x);
        x += h / 2;
        s += f(x);
    }

    s *= h / 6;
    return s;
}

std::vector<double> VariableStepRectangle(double a, double b, double e, double(*f)(double)) {
    std::vector<double> result;

    double h = b - a;
    double s, x, t2;
    double t1 = (h / 2) * (f(a) + f(b));
    while (true) {
        s = 0;
        x = a + (h / 2);
        while (x < b) {
            s += f(x);
            x += h;
        }

        t2 = (t1 / 2) + (h / 2) * s;

        result.push_back(t2);

        if (abs(t2 - t1) < e) return result;

        h /= 2;
        t1 = t2;
    }
}

std::vector<double> FixedStepRectangle(double a, double b, int step, double(*f)(double)) {
    std::vector<double> res;

    double h = b - a;
    double s, x, t2;
    double t1 = (h / 2) * (f(a) + f(b));
    for (int i = 1; i <= step; i++) {
        s = 0;
        x = a + (h / 2);
        while (x < b) {
            s += f(x);
            x += h;
        }

        t2 = (t1 / 2) + (h / 2) * s;

        res.push_back(t2);

        h /= 2;
        t1 = t2;
    }

    return res;
}

typedef struct RombergResult {
    std::vector<double> values;
    int step1 = 0; // 梯形公式的二分次数
} RombergResult;

RombergResult Romberg(double a, double b, double e, double(*f)(double)) {
    RombergResult result;

    double h = b - a;
    double t1 = (h / 2) * (f(a) + f(b));
    int k = 1;
    double s, s1, s2, x, t2, c1, c2, r1, r2;
    while (true) {
        s = 0;
        x = a + (h / 2);

        while (true) {
            s += f(x);
            x += h;
            if (x >= b) break;
        }

        t2 = (t1 / 2) + (h / 2) * s;
        result.step1++;
        s2 = t2 + (t2 - t1) / 3;

        if (k > 1) {
            c2 = s2 + (s2 - s1) / 15;
            if (k > 2) {
                r2 = c2 + (c2 - c1) / 63;
                if (k > 3) {
                    result.values.push_back(r2);
                    if (abs(r2 - r1) < e) return result;
                    else goto K3;
                } else goto K3;
            } else goto K2;
        } else goto K1;

        K3:
        r1 = r2;

        K2:
        c1 = c2;

        K1:
        k++;
        h /= 2;
        t1 = t2;
        s1 = s2;
    }
}

std::vector<std::pair<double, double>>
ImprovedEuler(double x0, double y0, double h, int N, double(*f)(double, double)) {
    std::vector<std::pair<double, double>> result;

    int n = 1;
    double x1, y1, yp, yc;
    while (n <= N) {
        x1 = x0 + h;
        yp = y0 + h * f(x0, y0);
        yc = y0 + h * f(x1, yp);
        y1 = (yp + yc) / 2;

        result.emplace_back(x1, y1);

        n++;
        x0 = x1;
        y0 = y1;
    }

    return result;
}

/**
 * 利用四阶龙格-库塔格式（经典格式）计算初值问题
 * @param x0
 * @param y0
 * @param h
 * @param N
 * @param f
 * @return
 */
std::vector<std::pair<double, double>>
RungeKutta4th(double x0, double y0, double h, int N, double(*f)(double, double)) {
    std::vector<std::pair<double, double>> result;

    int n = 1;
    double x1, y1, k1, k2, k3, k4;
    while (n <= N) {
        x1 = x0 + h;
        k1 = f(x0, y0);
        k2 = f(x0 + (h / 2), y0 + (h / 2) * k1);
        k3 = f(x0 + (h / 2), y0 + (h / 2) * k2);
        k4 = f(x1, y0 + h * k3);
        y1 = h * (k1 + 2 * k2 + 2 * k3 + k4) / 6 + y0;

        result.emplace_back(x1, y1);

        n++;
        x0 = x1;
        y0 = y1;
    }
    return result;
}

typedef struct AdamsResult {
    std::vector<std::pair<double, double>> predictors;
    std::vector<std::pair<double, double>> correctors;
} AdamsResult;

AdamsResult
AdamsPredictorCorrectorBasedOnRK4(double x0, double y0, double h, int N, double(*f)(double, double)) {
    double y[5];
    double x[5];
    double yy[5];

    AdamsResult result;

    y[0] = y0;
    x[0] = x0;
    yy[0] = 0;
    result.correctors.emplace_back(x0, y0);

    std::vector<std::pair<double, double>> rk = RungeKutta4th(x0, y0, h, 3, f);

    for (int i = 1; i <= 3; i++) {
        x[i] = rk[i - 1].first;
        y[i] = rk[i - 1].second;
        result.correctors.emplace_back(x[i], y[i]);
    }

    for (int i = 0; i <= 3; i++) {
        yy[i] = f(x[i], y[i]);
    }

    int n = 4;

    double yp, yyp;
    while (n <= N) {
        x[4] = x[3] + h;

        yp = y[3] + h * (55 * yy[3] - 59 * yy[2] + 37 * yy[1] - 9 * yy[0]) / 24;
        yyp = f(x[4], yp);

        result.predictors.emplace_back(x[4], yp);

        y[4] = y[3] + h * (9 * yyp + 19 * yy[3] - 5 * yy[2] + yy[1]) / 24;
        yy[4] = f(x[4], y[4]);

        result.correctors.emplace_back(x[4], y[4]);

        n++;
        x[3] = x[4];
        y[3] = y[4];
        for (int k = 0; k <= 3; k++) yy[k] = yy[k + 1];
    }

    return result;
}
