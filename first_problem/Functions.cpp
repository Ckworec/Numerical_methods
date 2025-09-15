#include "Header.h"
#include "Vector.h"

double f1(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1], x3 = x.v[2];
	// return x1 * x1 + x2 * x2;
	// return std::cos(x1 * x1 + 10 * x2 * x2 + 0.1 * x3 * x3);
	// return exp(x1 * x1 + 10 * x2 * x2 + 0.1 * x3 * x3);
	// return (x1 - 1) * (x1 - 1) + 10 * (x2 - 2) * (x2 - 2) + 0.1 * (x3 - 3) * (x3 - 3);
	return -cos((x1 * x1 - 2 * x1 * x2 + 2 * x2 * x2 + 4 * x2 + 4)/100);
}