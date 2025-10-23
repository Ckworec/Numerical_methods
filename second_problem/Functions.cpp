#include "Header.h"
#include "Fredholm.h"

Fredholm Fredholm1()
{
	auto u = [](double x) { return x + exp(-x); };
	auto K = [](double x, double t) { return 0.5 * x * exp(t); };
	auto f = [](double x) { return exp(-x); };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}

Fredholm Fredholm2()
{
	auto u = [](double x) { return 1.0; };
	auto K = [](double x, double t) { return sin(x * t); };
	auto f = [](double x) { return 1.0 + (cos(x/2.0) - 1.0) / x; };
	double a = 0.0, b = 0.5;
	return Fredholm(u, K, f, a, b);
}

Fredholm Fredholm3()
{
	auto u = [](double x) { return (25.0 + 27.0 * cos(2.0*x)) / (160.0 * M_PI); };
	auto K = [](double x, double t) { 
		return 1.0 / (sin((x+t)/2.0) * sin((x+t)/2.0) + 0.25 * cos((x+t)/2.0) * cos((x+t)/2.0));
	};
	auto f = [](double x) { return (5.0 + 3.0 * cos(2.0*x)) / (16.0 * M_PI); };
	double a = 0.0, b = 2.0 * M_PI;
	return Fredholm(u, K, f, a, b);
}

Fredholm Fredholm5()
{
	auto u = [](double x) { return cos(2.0*x); };
	auto K = [](double x, double t) { return sin(x) * cos(t); };
	auto f = [](double x) { return cos(2.0*x); };
	double a = 0.0, b = 2.0 * M_PI;
	return Fredholm(u, K, f, a, b);
}

Fredholm Fredholm7()
{
	auto u = [](double x) { return 1.0 + 4.0*x/9.0; };
	auto K = [](double x, double t) { return x * t * t; };
	auto f = [](double x) { return 1.0; };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}

Fredholm Fredholm8()
{
	auto u = [](double x) { return x; };
	auto K = [](double x, double t) { return 0.5 * x * t; };
	auto f = [](double x) { return 5.0 * x / 6.0; };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}

Fredholm Fredholm9()
{
	auto u = [](double x) { return 1.0; };
	auto K = [](double x, double t) { return x*x * exp(x*t); };
	auto f = [](double x) { return 1.0 - x * (exp(x) - exp(-x)); };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}