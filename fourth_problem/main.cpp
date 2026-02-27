#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const long double pi = 3.141592653589793;

// Параметры сетки (легко менять)
int N1 = 20;
int N2 = 20;
long double x0 = 0.0, X = 1.0;
long double Y0 = 0.0, Y = 1.0;
long double eps = 0.01;

// Коэффициент a для матрицы K
long double a = 0.0;

// Шаги сетки
long double h1, h2;

// Функция правой части
long double f(long double x, long double y) {
    // return -pi*pi*sin(pi*x)*sin(pi*y) - pi*pi*sin(pi*x)*sin(pi*y);
    return -2 * (y * (y - 1)) - 2 * (x * (x - 1));
}

// Матрица коэффициентов K
vector<vector<long double>> K(long double x, long double y) {
    return { {1.0, a}, {a, 1.0} };
}

// Норма вектора
long double norm(const vector<long double>& vec) {
    long double res = 0;
    for (auto v : vec) res += v*v;
    return sqrt(res);
}

// Предобуславливатель
long double v_ij(int i, int j, const vector<long double>& F_) {
    long double res = 0.0;
    for (int k2 = 1; k2 < N2-1; k2++) {
        long double nu = 0;
        for (int k1 = 1; k1 < N1-1; k1++) {
            long double sum = 0.0;
            for (int ii = 0; ii < N1-2; ii++) {
                for (int jj = 0; jj < N2-2; jj++) {
                    sum += F_[ii + (N1-2)*jj] * sin(k1*pi*(ii+1)/(N1-1)) * sin(k2*pi*(jj+1)/(N2-1));
                }
            }
            long double lambda = 4/(h1*h1)*sin(k1*pi*h1/(2*(X-x0)))*sin(k1*pi*h1/(2*(X-x0)))
                               + 4/(h2*h2)*sin(k2*pi*h2/(2*(Y-Y0)))*sin(k2*pi*h2/(2*(Y-Y0)));
            nu += sum * sin(k1*pi*(i+1)/(N1-1)) / lambda;
        }
        res += nu * sin(k2*pi*(j+1)/(N2-1));
    }
    return res * 4/((N1-1)*(N2-1));
}

// Оператор A*v (центральные разности)
vector<long double> applyA(const vector<long double>& v) {
    int N = (N1-2)*(N2-2);
    vector<long double> res(N, 0.0);

    for (int j = 0; j < N2-2; j++) {
        for (int i = 0; i < N1-2; i++) {
            int idx = i + (N1-2)*j;
            long double vC = v[idx];
            long double vL = (i>0) ? v[idx-1] : 0.0;
            long double vR = (i<N1-3) ? v[idx+1] : 0.0;
            long double vD = (j>0) ? v[idx-(N1-2)] : 0.0;
            long double vU = (j<N2-3) ? v[idx+(N1-2)] : 0.0;

            auto k = K(i*h1,j*h2);
            res[idx] = (k[0][0]*(vL - 2*vC + vR)/(h1*h1) + k[1][1]*(vD - 2*vC + vU)/(h2*h2)
                        + (k[0][1]+k[1][0])*(0.0)); // Смешанные производные
        }
    }
    return res;
}

int main() {
    h1 = (X - x0) / (N1-1);
    h2 = (Y - Y0) / (N2-1);
    int N = (N1-2)*(N2-2);

    // Формируем F
    vector<long double> F(N, 0.0);
    for (int j=0;j<N2-2;j++)
        for (int i=0;i<N1-2;i++)
            F[i + (N1-2)*j] = f((i+1)*h1,(j+1)*h2);

    // Инициализация решения
    vector<long double> v_n(N, 0.0);
    vector<long double> r(N, 1.0);
    vector<long double> tn(N, 0.0);

    int count = 0;
    while (norm(r) > eps) {
        // r = A*v - F
        auto Av = applyA(v_n);
        for (int i=0;i<N;i++) r[i] = Av[i] - F[i];

        // tn = B^-1 * r
        for (int j=0;j<N2-2;j++)
            for (int i=0;i<N1-2;i++)
                tn[i + (N1-2)*j] = v_ij(i,j,r);

        // alpha = (r,tn)/(A*tn,tn)
        long double numerator = 0.0;
        for (int i=0;i<N;i++) numerator += r[i]*tn[i];

        auto Atn = applyA(tn);
        long double denominator = 0.0;
        for (int i=0;i<N;i++) denominator += Atn[i]*tn[i];

        long double alpha = numerator/denominator;

        // v = v - alpha*tn
        for (int i=0;i<N;i++) v_n[i] -= alpha*tn[i];

        count++;
        cout << "Iteration " << count << ", norm(r) = " << norm(r) << endl;
    }

    cout << "Converged in " << count << " iterations." << endl;
    return 0;
}