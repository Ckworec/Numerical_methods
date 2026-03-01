#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

const long double pi = 3.14159265358979323846;

// Параметры сетки (увеличил для демонстрации скорости)
int N1 = 10, N2 = 10;
long double x0 = 0.0, X = 1.0, Y0 = 0.0, Y = 1.0;
long double eps = 1e-8;
long double a = 0.0; // Коэффициент смешанной производной
long double h1, h2;

long double f(long double x, long double y) {
    return (y * (y - 1)) * (x * (x - 1));
}

// Оператор A (со знаком МИНУС для положительной определенности)
vector<long double> applyA(const vector<long double>& v) {
    int m = N1 - 2, n = N2 - 2;
    vector<long double> res(m * n, 0.0);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            int idx = i + m * j;
            long double vC = v[idx];
            long double vL = (i > 0) ? v[idx - 1] : 0.0;
            long double vR = (i < m - 1) ? v[idx + 1] : 0.0;
            long double vD = (j > 0) ? v[idx - m] : 0.0;
            long double vU = (j < n - 1) ? v[idx + m] : 0.0;

            // Смешанные производные (центральная разность "крест")
            long double vLU = (i > 0 && j < n - 1) ? v[idx - 1 + m] : 0.0;
            long double vRU = (i < m - 1 && j < n - 1) ? v[idx + 1 + m] : 0.0;
            long double vLD = (i > 0 && j > 0) ? v[idx - 1 - m] : 0.0;
            long double vRD = (i < m - 1 && j > 0) ? v[idx + 1 - m] : 0.0;

            long double d2x = (vL - 2 * vC + vR) / (h1 * h1);
            long double d2y = (vD - 2 * vC + vU) / (h2 * h2);
            long double d2xy = (vRU - vLU - vRD + vLD) / (4 * h1 * h2);

            res[idx] = -(d2x + d2y + 2 * a * d2xy);
        }
    }
    return res;
}

// Быстрое (относительно) обращение оператора Лапласа через DST
vector<long double> solvePoisson(const vector<long double>& r) {
    int m = N1 - 2, n = N2 - 2;
    vector<long double> res(m * n, 0.0);
    vector<long double> tmp(m * n, 0.0);

    // Прямое преобразование (Синус-базис)
    for (int k2 = 1; k2 <= n; k2++) {
        for (int k1 = 1; k1 <= m; k1++) {
            long double sum = 0;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    sum += r[i + m * j] * sin(k1 * pi * (i + 1) / (N1 - 1)) * sin(k2 * pi * (j + 1) / (N2 - 1));
                }
            }
            // Собственные числа оператора -Delta
            long double lam1 = (4.0 / (h1 * h1)) * pow(sin(k1 * pi / (2.0 * (N1 - 1))), 2);
            long double lam2 = (4.0 / (h2 * h2)) * pow(sin(k2 * pi / (2.0 * (N2 - 1))), 2);
            tmp[(k1 - 1) + m * (k2 - 1)] = sum / (lam1 + lam2);
        }
    }

    // Обратное преобразование
    long double coeff = 4.0 / ((N1 - 1) * (N2 - 1));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            long double val = 0;
            for (int k2 = 1; k2 <= n; k2++) {
                for (int k1 = 1; k1 <= m; k1++) {
                    val += tmp[(k1 - 1) + m * (k2 - 1)] * sin(k1 * pi * (i + 1) / (N1 - 1)) * sin(k2 * pi * (j + 1) / (N2 - 1));
                }
            }
            res[i + m * j] = val * coeff;
        }
    }
    return res;
}

long double dot(const vector<long double>& a, const vector<long double>& b) {
    long double res = 0;
    for (size_t i = 0; i < a.size(); i++) res += a[i] * b[i];
    return res;
}

int main() {
    h1 = (X - x0) / (N1 - 1); h2 = (Y - Y0) / (N2 - 1);
    int N = (N1 - 2) * (N2 - 2);
    vector<long double> F(N), v_n(N, 0.0);

    for (int j = 0; j < N2 - 2; j++)
        for (int i = 0; i < N1 - 2; i++)
            F[i + (N1 - 2) * j] = -f((i + 1) * h1, (j + 1) * h2); // F тоже с минусом

    int count = 0;
    while (true) {
        auto Av = applyA(v_n);
        vector<long double> r(N);
        for (int i = 0; i < N; i++) r[i] = Av[i] - F[i];

        long double r_norm = sqrt(dot(r, r));
        cout << "Iter " << count << ", norm: " << r_norm << endl;
        if (r_norm < eps || count > 1000) break;

        vector<long double> tn = solvePoisson(r);
        auto Atn = applyA(tn);

        long double alpha = dot(r, tn) / dot(Atn, tn);
        cout << "Alpha: " << fixed << setprecision(12) << alpha << endl;
        for (int i = 0; i < N; i++) v_n[i] -= alpha * tn[i];
        count++;
    }

    return 0;
}
