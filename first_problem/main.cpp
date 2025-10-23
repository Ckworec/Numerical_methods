#include "Header.h"
#include "Task1.h"

int main()
{
	double (*f)(Vector x);
	f = f1;
	double epsilon = 1e-6;
	int act_size;
	cout << "Enter problem size: ";
	cin >> act_size;
	vector<double> x00(act_size);

	cout << "Enter start point: ";
	for (int i = 0; i < x00.size(); i++) 
		cin >> x00[i];

	Vector x0(x00); 
	Vector x(x0.GetSize());

	Task1::FirstMethod(f, epsilon, x0, x);
	Task1::SecondMethod(f, epsilon, x, x);

	return 0;
}