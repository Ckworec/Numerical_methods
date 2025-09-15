#include "Header.h"
#include "Matrix.h"

namespace Task1
{
	double StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk);
	double goldenSection(double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon);
	void FirstMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x); // запуск первого метода
	void SecondMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x); // запуск второго метода
}
