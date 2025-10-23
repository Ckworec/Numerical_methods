#include "Header.h"

ostream& operator<<(ostream& cout, const vector<double>& b)
{
	if (b.size() != 0) cout << "("; else return cout;
	for (int i = 0; i < b.size() - 1; i++) cout << b[i] << ", ";
	cout << b.back() << ")" << endl;
	return cout;
}