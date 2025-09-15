#include "Task1.h"


double Task1::StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk)
{
	auto phi = [&](double t) { return f(xk + hk * t); };
	double lambda = 0.5, mu = 2;
	double alpha, beta;
	double alpha_k;

	double alpha_max = 10;

	beta = 0.1;
	alpha = beta;

	if (phi(alpha) < f(xk))
	{
		do
		{
			alpha_k = alpha;
			alpha *= mu;
		} 
		while (phi(alpha) < phi(beta) && alpha <= alpha_max);
	}
	else
	{
		do
		{
			alpha_k = alpha;
			alpha *= lambda;
		} 
		while (!(phi(alpha) < f(xk)));
	}

	return alpha_k;
}

double Task1::goldenSection(double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon) {
    const double r = (3.0 - std::sqrt(5.0)) / 2.0;
	auto phi = [&](double t) { return f(xk + hk * t); };

	double a = 0;
	double b = 1;

	while (phi(b) < phi(a)) 
	{
		a = b;
		b *= 2;
	}

    double c = a + r * (b - a);
    double d = b - r * (b - a);
    double fc = phi(c);
    double fd = phi(d);
    int iter = 0;
    while ((b - a) > epsilon) {
        if (fc < fd) {
            // минимум в [a, d]
            b = d;
            d = c;
            fd = fc;
            c = a + r * (b - a);
            fc = phi(c);
        } else {
            // минимум в [c, b]
            a = c;
            c = d;
            fc = fd;
            d = b - r * (b - a);
            fd = phi(d);
        }
        ++iter;
    }

    return 0.5 * (a + b);
}

void Task1::FirstMethod(double(*f)(Vector x), const double epsilon, Vector& x0, Vector& x)
{
	size_t n = x0.GetSize();
	Vector xk(n), xk1(n);
	xk = x0;
	vector<double> df(n);
	Vector Df(n);
	Vector hk(n);
	Vector diff_iter(n);
	double alpha_k;
	double diff_f;
	double (*norm_sqr)(Vector& x);
	norm_sqr = l2_norm_square;
	double difference_norm_sqr, grad_norm_sqr;

	const double tau = 0.1 * sqrt(epsilon);

	vector<vector<double>> d(n, vector<double>(n));

	for (int i = 0; i < n; i++) 
		d[i][i] = tau;

	int iterations = 0;

	do
	{
		// Считаем градиент
		for (int i = 0; i < n; i++)
		{
			df[i] = (f(xk + d[i]) - f(xk - d[i])) / (2 * tau);
		}

		Df = df;
		hk = Df * (-1);
		alpha_k = StepSplittingMethod(f, xk, hk);
		xk1 = xk + hk * alpha_k;
		diff_iter = xk1 - xk;
		difference_norm_sqr = norm_sqr(diff_iter);
		diff_f = std::fabs(f(xk1) - f(xk));
		grad_norm_sqr = norm_sqr(Df);

		// cout << "Df = " << Df;
		// cout << "alpha_k = " << alpha_k << endl;
		// cout << "xk = " << xk;
		// cout << "xk+1 = " << xk1;
		// cout << difference_norm_sqr << "\t" << diff_f << "\t" << grad_norm_sqr << endl << endl;

		xk = xk1;
		iterations++;
	} while (difference_norm_sqr > epsilon && diff_f * diff_f > epsilon && grad_norm_sqr > epsilon);

	x = xk1;
	cout << "Iterations first method: " << iterations << endl << endl;
	cout << "Result: x = ";
	cout << x << endl;
}

void Task1::SecondMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x)
{
	size_t n = x0.GetSize();
	Vector xk(n), xk1(n);
	xk = x0;
	vector<double> grad_f(n);
	Vector Df(n);
	Vector hk(n);
	Vector diff_iter(n);
	double alpha_k;
	double diff_f;
	double (*norm_sqr)(Vector& x);
	norm_sqr = l2_norm_square;
	double difference_norm_sqr, grad_norm_sqr;

	const double tau = 0.1 * sqrt(epsilon);

	Matrix He(n, n);

	vector<vector<double>> d(n, vector<double>(n));
	for (int i = 0; i < n; i++) 
		d[i][i] = tau;

	vector<vector<double>> der2f(n, vector<double>(n));

	vector<std::function<double(Vector)>> df(n);
	vector<std::function<double(Vector)>> df4(n);
	double fij4, fi4j;
	
	int iterations = 0;
	do
	{
		// Считаем градиент
		for (int i = 0; i < n; i++)
		{
			grad_f[i] = (f(xk + d[i]) - f(xk - d[i])) / (2 * tau);
		}

		Df = grad_f;

		// Считаем матрицу Гесса
		for (int i = 0; i < n; i++)
		{
			df[i] = [&](Vector p) { return (f(p + d[i]) - f(p)) / tau; };
			df4[i] = [&](Vector p) { return (f(p) - f(p - d[i])) / tau; };
			der2f[i][i] = (f(xk + d[i]) - 2 * f(xk) + f(xk - d[i])) / (tau * tau);
			for (int j = i + 1; j < n; j++)
			{
				fij4 = (df[i](xk) - df[i](xk - d[j])) / tau;
				fi4j = (df4[i](xk + d[j]) - df4[i](xk)) / tau;
				der2f[i][j] = der2f[j][i] = 0.5 * (fij4 + fi4j);
			}
		}

		He.Move(der2f);
		// Решаем систему методом Гаусса
		hk = He.Gauss(Df) * (-1);
		alpha_k = goldenSection(f, xk, hk, epsilon);
		xk1 = xk + hk * alpha_k;
		diff_iter = xk1 - xk;
		difference_norm_sqr = norm_sqr(diff_iter);
		diff_f = std::fabs(f(xk1) - f(xk));
		grad_norm_sqr = norm_sqr(Df);

		// cout << "xk = " << xk;
		// cout << "alpha_k = " << alpha_k << endl;
		// cout << "Df = " << Df << endl;

		xk = xk1;
		iterations++;
	} while (difference_norm_sqr > epsilon*epsilon && diff_f > epsilon && grad_norm_sqr > epsilon*epsilon);

	x = xk1;
	cout << "Iterations second method: " << iterations << endl << endl;
	cout << "Answer: x = " << x << endl;
	cout << "f(x) = " << f(x) << endl;

}