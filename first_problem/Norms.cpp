#include "Header.h"
#include "Vector.h"

double l1_norm(Vector& x)
{
	return std::accumulate(x.v.begin(), x.v.end(), 0, [](double total, double value) { return total + std::abs(value); });
}
double l2_norm(Vector& x)
{
	return sqrt(std::accumulate(x.v.begin(), x.v.end(), 0, [](double total, double value) { return total + value * value; }));
}
double l2_norm_square(Vector& x)
{
	//return std::accumulate(x.v.begin(), x.v.end(), 0, [](double total, double value) { return total + value * value; });
	double s = 0; for_each(x.v.begin(), x.v.end(), [&s](double c) {s += c*c; });
	for (int i = 0; i < x.v.size(); i++) s += x.v[i] * x.v[i];
	return s;
}
double l_inf_norm(Vector& x)
{
	return *std::max_element(x.v.begin(), x.v.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
}