#include "Vector.h"

Vector::Vector(const int size)
{
	v.resize(size);
}

Vector::Vector(vector<double>& b)
{
	v = b;
}

Vector::Vector(const Vector& b)
{
	v = b.v;
}

Vector& Vector::operator=(const Vector& b)
{
	if (this != &b)
	{
		v.clear();
		v.resize(b.v.size());
		for (int i = 0; i < v.size(); i++) 
			v[i] = b.v[i];
	}
	return *this;
}

Vector& Vector::operator=(const vector<double>& vec)
{
	v.clear();
	v.resize(vec.size());
	for (int i = 0; i < v.size(); i++) 
		v[i] = vec[i];
	return *this;
}

Vector Vector::operator+(const Vector& b)
{
	Vector s(v.size());
	for (int i = 0; i < s.v.size(); i++) 
		s.v[i] = v[i] + b.v[i];
	return s;
}

Vector Vector::operator+(const double& b)
{
	Vector s(v.size());
	for (int i = 0; i < s.v.size(); i++) 
		s.v[i] = v[i] + b;
	return s;
}

Vector Vector::operator-(const Vector& b)
{
	Vector s(v.size());
	for (int i = 0; i < s.v.size(); i++) 
		s.v[i] = v[i] - b.v[i];
	return s;
}

Vector Vector::operator-(const double& b)
{
	Vector s(v.size());
	for (int i = 0; i < s.v.size(); i++) 
		s.v[i] = v[i] - b;
	return s;
}

Vector Vector::operator*(const double a)
{
	Vector s(v.size());
	for (int i = 0; i < s.v.size(); i++) 
		s.v[i] = v[i] * a;
	return s;
}

double Vector::operator*(const Vector& b)
{
	double sp = 0;
	for (int i = 0; i < v.size(); i++) 
		sp += v[i] * b.v[i];
	return sp;
}

size_t Vector::GetSize()
{
	return v.size();
}

ostream& operator<<(ostream& cout, Vector& b)
{
	int size = b.v.size();
	if (size != 0) 
		cout << "("; 
	else
		return cout;

	for (int i = 0; i < size - 1; i++)
		cout << b.v[i] << ", ";

	cout << b.v[size - 1] << ")" << endl;

	return cout;
}