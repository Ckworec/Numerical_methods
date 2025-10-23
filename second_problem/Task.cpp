#include "Task.h"
#include "Matrix.h"
#include "Fredholm.h"
#define N 1000000

void Task::Count(double (*K)(double x, double t), double (*f)(double x), const double alpha, const double beta, vector<double>& cn, vector<double>& Un, vector<double>& tn, const int k, Matrix& Mat_glob, vector<double>& F_glob)
{

	vector<double> c = { (beta - alpha) / 2.0, (beta - alpha) / 2.0 };
    vector<double> t = { alpha, beta };

    // первый узел сохраняем только один раз
    if (k == 0) tn[0] = t[0];

    // добавляем веса и точки
    cn[k]     += c[0];
    cn[k + 1] += c[1];

    tn[k + 1] = t[1];
}

void Task::FillMatrix(double (*K)(double x, double t), double (*f)(double x), const vector<double>& cn, const vector<double>& tn, Matrix& Mat_glob)
{
	int n = Mat_glob.a.size();
	for (int i = 0; i < n; i++)
	{
		Mat_glob.a[i][i] += 1.0;
		for (int j = 0; j < n; j++)
		{
			Mat_glob.a[i][j] -= cn[j] * K(tn[i], tn[j]);
		}
	}
}

double Task::Integral(std::function<double(double)> f, double a, double b, double n, double epsilon)
{
	double I = 0.0;
	double c, d, h;

	h = (b-a)/n;
	c = a; 
	d = c + h;

	for (int i = 0; i < n; i++)
	{
		I += (d - c)/2.0 * (f(c) + f(b));
		c = d; 
		d += h;
	}

	return I;
}

double Task::L2_norm_sqr(std::function<double(double)> g, const double a, const double b, const double epsilon)
{
	auto f = [&](double x) { return pow(g(x), 2); };
	const int p = 4;
	int n_global = 10;
	double n;
	double h = (b-a)/n_global;
	double MODIFIED_epsilon = epsilon * h / (b-a);
	double rh;
	double I = 0.0;
	double Sh, Sh2;
	double alpha, beta = a;
	while (fabs(b-beta) > epsilon)
	{
		alpha = beta; beta = std::min(beta + 2*(b-a)/n_global, b);

		n = 1;
		h = (b-a)/n_global;
		Sh = Integral(f, alpha, beta, n, epsilon);		// h
		Sh2 = Integral(f, alpha, beta, 2*n, epsilon);	// h/2
		rh = std::fabs(Sh - Sh2) / (1 - pow(2, 1 - p));
		while (rh > MODIFIED_epsilon)
		{
			Sh = Sh2;
			n *= 2;
			h /= 2;
			Sh2 = Integral(f, alpha, beta, 2*n, epsilon);
			rh = std::fabs(Sh - Sh2) / (1 - pow(2, 1 - p));
			MODIFIED_epsilon = epsilon * h / (b - a);
		}
		I += Sh;
	}
	return I;
}

void Task::main_func()
{
    Fredholm Eq = Fredholm8();
    double a = Eq.a, b = Eq.b;
    double (*u)(double x) = Eq.u;
    double (*K)(double x, double t) = Eq.K;
    double (*f)(double x) = Eq.f;
    const double epsilon = 1e-2; cout << "epsilon = " << epsilon << endl << endl;

    // m = число подотрезков; при трапециях n = m + 1 узлов
    int m = 1;
    int n = m + 1;

    vector<double> cn1, cn2;
    vector<double> Un1, Un2;
    vector<double> tn1, tn2;

    Matrix Mat_glob(n, n);
    vector<double> F_glob(n);

    // инициализация для начальной сетки
    cn2.assign(n, 0.0);
    Un2.assign(n, 0.0);
    tn2.assign(n, 0.0);
    F_glob.assign(n, 0.0);

    tn2[0] = a;
    double h = (b - a) / m;
    for (int k = 0; k < m; ++k)
    {
        double alpha = a + k * h;
        double beta  = alpha + h;
        Count(K, f, alpha, beta, cn2, Un2, tn2, k, Mat_glob, F_glob);
    }

    // Правую часть задаём напрямую как f(t_i)
    for (int i = 0; i < n; ++i) F_glob[i] = f(tn2[i]);

    FillMatrix(K, f, cn2, tn2, Mat_glob);
    Un2 = Mat_glob.Gauss(F_glob);

    std::function<double(double)> un, u2n;
    u2n = [&](double x) {
        double S = 0.0;
        for (int j = 0; j < (int)cn2.size(); ++j) S += cn2[j] * K(x, tn2[j]) * Un2[j];
        return S + f(x);
    };

    double norm_sqr;
    int n_prev;

    // цикл уточнения: удваиваем число подотрезков m -> n = m + 1
    do
    {
        cn1 = std::move(cn2);
        Un1 = std::move(Un2);
        tn1 = std::move(tn2);
        un = [&](double x) {
            double S = 0.0;
            for (int j = 0; j < (int)cn1.size(); ++j) S += cn1[j] * K(x, tn1[j]) * Un1[j];
            return S + f(x);
        };

        n_prev = n;
        m *= 2;            // удваиваем число подотрезков
        n = m + 1;         // соответственно число узлов

        Mat_glob = Matrix(n, n);
        cn2.assign(n, 0.0);
        Un2.assign(n, 0.0);
        tn2.assign(n, 0.0);
        F_glob.assign(n, 0.0);

        tn2[0] = a;
        h = (b - a) / m;
        for (int k = 0; k < m; ++k)
        {
            double alpha = a + k * h;
            double beta  = alpha + h;
            Count(K, f, alpha, beta, cn2, Un2, tn2, k, Mat_glob, F_glob);
        }

        // Правую часть снова задаём напрямую
        for (int i = 0; i < n; ++i) F_glob[i] = f(tn2[i]);

        FillMatrix(K, f, cn2, tn2, Mat_glob);
        Un2 = Mat_glob.Gauss(F_glob);

        u2n = [&](double x) {
            double S = 0.0;
            for (int j = 0; j < (int)cn2.size(); ++j) S += cn2[j] * K(x, tn2[j]) * Un2[j];
            return S + f(x);
        };

        auto diff = [&](double x) { return un(x) - u2n(x); };
        norm_sqr = L2_norm_sqr(diff, a, b, epsilon);
        cout << "n_prev = " << n_prev << "  (m_prev = " << n_prev-1 << ")" << endl;
        cout << "\t||un-u2n|| = " << sqrt(norm_sqr) << endl;

    } while (norm_sqr > epsilon * epsilon);

    auto diffFinal = [&](double x) {
        return u(x) - u2n(x);
    };
    norm_sqr = L2_norm_sqr(diffFinal, a, b, epsilon);

    cout << "\t||u-u2n|| = " << sqrt(norm_sqr) << endl;

    ofstream out("out.txt");
    for (auto& tt : tn2) out << tt << " "; out << endl;
    for (auto& tt : tn2) out << u2n(tt) << " "; out << endl;
}
