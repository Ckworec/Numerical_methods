#include "Header.h"
#include "Task1.h"

int main()
{
	double (*f)(Vector x);
	f = f1;
	double epsilon = 1e-8;
	vector<double> x00(3);

	cout << "Enter start point: " << endl;
	for (int i = 0; i < x00.size(); i++) 
		cin >> x00[i];

	Vector x0(x00); 
	Vector x(x0.GetSize());

	Task1::FirstMethod(f, epsilon, x0, x);
	Task1::SecondMethod(f, epsilon, x, x);

	return 0;
}