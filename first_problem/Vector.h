#include "Header.h"

class Vector
{
protected:
	vector<double> v;
public:
	Vector(const int size);
	Vector(vector<double>& b);
	Vector(const Vector& b);
	Vector& operator=(const Vector& b);
	Vector& operator=(const vector<double>& vec);
	Vector operator+(const Vector& b);
	Vector operator+(const double& b);
	Vector operator-(const Vector& b);
	Vector operator-(const double& b);
	Vector operator*(const double a);
	double operator*(const Vector& b);

	size_t GetSize();

	friend class Matrix;

	friend double l1_norm(Vector& x);
	friend double l2_norm(Vector& x);
	friend double l2_norm_square(Vector& x);
	friend double l_inf_norm(Vector& x);

	friend double f1(Vector x);
	friend double f2(Vector x);
	friend double f3(Vector x);
	friend double f4(Vector x);
	friend double f5(Vector x);
	friend double f6(Vector x);
	friend double f7(Vector x);
	friend double f8(Vector x);
	friend double f9(Vector x);
	friend double f10(Vector x);
	friend double f11(Vector x);
	friend double f12(Vector x);

	friend ostream& operator<<(ostream& cout, Vector& b);
};