#include "Fredholm.h"

Fredholm::Fredholm(double (*u)(double x), double (*K)(double x, double t), double (*f)(double x), double a, double b)
{ 
	if (a >= b) throw - 1;
	this->u = u; this->K = K; this->f = f; this->a = a; this->b = b;
}
Fredholm& Fredholm::operator=(const Fredholm& other)
{
	if (this != &other)
	{
		u = other.u; K = other.K; f = other.f; a = other.a; b = other.b;
	}
	return *this;
}
double Fredholm::GetU(double x)
{
	return u(x);
}
double Fredholm::GetK(double x, double t)
{
	return K(x, t);
}
double Fredholm::GetF(double x)
{
	return f(x);
}