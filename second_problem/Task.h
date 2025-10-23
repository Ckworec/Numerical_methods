#include "Header.h"

namespace Task
{
	double Integral(std::function<double(double)> f, double a, double b, double n, double epsilon);
	void Count(double (*K)(double x, double t), double (*f)(double x), const double a, const double b, vector<double>& cn, vector<double>& Un, vector<double>& tn, const int i, Matrix& Mat_glob, vector<double>& F_glob);
	void FillMatrix(double (*K)(double x, double t), double (*f)(double x), const vector<double>& cn, const vector<double>& tn, Matrix& Mat_glob);

	double L2_norm_sqr(std::function<double(double)> f, const double a, const double b, const double epsilon);

	void main_func();
}
