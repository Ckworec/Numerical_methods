#include "Header.h"

struct Fredholm
{
	double (*u)(double x);
	double (*K)(double x, double t);
	double (*f)(double x);
	double a, b;
	Fredholm(double (*u)(double x), double (*K)(double x, double t), double (*f)(double x), double a, double b);
	Fredholm& operator=(const Fredholm& other);
	double GetU(double x);
	double GetK(double x, double t);
	double GetF(double x);
};