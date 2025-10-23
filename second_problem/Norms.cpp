#include "Header.h"

double l2_norm_square(vector<double>& x)
{
	//return std::accumulate(x.begin(), x.end(), 0, [](double total, double value) { return total + value * value; });
	double s = 0; for_each(x.begin(), x.end(), [&s](double c) {s += c * c; });
	for (int i = 0; i < x.size(); i++) s += x[i] * x[i];
	return s;
}