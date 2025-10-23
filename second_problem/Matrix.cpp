#include "Matrix.h"

Matrix::Matrix() { N = M = 0; }
Matrix::Matrix(int N, int M)
{
	this->N = N; this->M = M;
	a.resize(N);
	for (auto& r : a) r.resize(M);
}
Matrix::Matrix(vector<row> rows)
{
	this->N = rows.size(); this->M = rows[0].size();
	a.resize(rows.size());
	for (int i = 0; i < N; i++) a[i].resize(M);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			a[i][j] = rows[i][j];
		}
	}
}
Matrix::Matrix(const Matrix& b)
{
	this->N = b.N;
	this->M = b.M;
	this->a.resize(N);
	for (int i = 0; i < N; i++)
	{
		this->a[i].resize(M);
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			this->a[i][j] = b.a[i][j];
		}
	}
}
size_t Matrix::GetN() { return N; }
size_t Matrix::GetM() { return M; }
void Matrix::SetValue(int i, int j, double value) { a[i][j] = value; }
Matrix& Matrix::operator=(const Matrix& b)
{
	if (this != &b)
	{
		N = b.N;
		M = b.M;
		a.resize(N);
		for (int i = 0; i < N; i++) a[i].resize(M);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				a[i][j] = b.a[i][j];
			}
		}
	}
	return *this;
}
Matrix& Matrix::operator=(Matrix&& b) noexcept
{
	if (this != &b)
	{
		a.clear();
		N = b.N;
		M = b.M;
		a.resize(N);
		for (int i = 0; i < N; i++) a[i].resize(M);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				a[i][j] = b.a[i][j];
			}
		}
		b.N = b.M = 0; b.a.clear();
	}
	return *this;
}
Matrix Matrix::operator+(const Matrix& b)
{
	if ((N != b.N) || (M != b.M))
	{
		//cout << "No" << endl;
		return Matrix();
	}
	
	Matrix S(N, M);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			S.a[i][j] = a[i][j] + b.a[i][j];
		}
	}

	return S;
}
Matrix Matrix::operator-(const Matrix& b)
{
	if ((N != b.N) || (M != b.M))
	{
		//cout << "No" << endl;
		return Matrix();
	}

	Matrix S(N, M);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			S.a[i][j] = a[i][j] - b.a[i][j];
		}
	}

	return S;
}
Matrix Matrix::operator*(double k)
{
	Matrix S(N, M);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			S.a[i][j] = k * a[i][j];
		}
	}
	return S;
}
Matrix Matrix::operator*(const Matrix& b)
{
	if (M != b.N)
	{
		//cout << "No" << endl;
		return Matrix();
	}
	Matrix S(N, b.M);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < b.M; j++)
		{
			S.a[i][j] = a[i][0] * b.a[0][j];
			for (int k = 1; k < M; k++)
			{
				S.a[i][j] += a[i][k] * b.a[k][j];
			}
		}
	}
		return S;
}
vector<double> Matrix::operator*(vector<double>& b)
{
	if (M != b.size())
	{
		//cout << "No" << endl;
		throw - 100;
	}
	vector<double> s(N);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			s[i] += a[i][j] * b[j];
		}
	}
	return s;
}
Matrix Matrix::T()
{
	Matrix S(M, N);
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			S.a[i][j] = a[j][i];
		}
	}
	return S;
}
double Matrix::Det2()
{
	return a[0][0] * a[1][1] - a[1][0] * a[0][1];
}
Matrix Matrix::Inv2()
{
	Matrix Y(2, 2);
	double detinv = 1 / (a[0][0]*a[1][1] - a[1][0]*a[0][1]);
	Y.a[0][0] = a[1][1];
	Y.a[1][1] = a[0][0];
	Y.a[1][0] = -a[1][0];
	Y.a[0][1] = -a[0][1];
	return Y * detinv;
}
void Matrix::Inversed(Matrix& Inv)
{
	Matrix Y(N, 2 * N);
	double temp;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Y.a[i][j] = a[i][j];
		}
		Y.a[i][N + i] = 1;
	}
	for (int i = N - 1; i > 0; i--)
	{
		if (Y.a[i - 1][0] < Y.a[i][0])
		{
			row temp = Y.a[i];
			Y.a[i] = Y.a[i - 1];
			Y.a[i - 1] = temp;
		}
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (j != i)
			{
				temp = Y.a[j][i] / Y.a[i][i];
				for (int k = 0; k < 2 * N; k++)
				{
					Y.a[j][k] -= Y.a[i][k] * temp;
				}
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		temp = Y.a[i][i];
		for (int j = 0; j < 2 * N; j++)
		{
			Y.a[i][j] = Y.a[i][j] / temp;
		}
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Inv.a[i][j] = Y.a[i][N + j];
		}
	}
}

void Matrix::Move(vector<row> rows)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			a[i][j] = rows[i][j];
		}
	}
}

void Matrix::SetZero()
{
	for (int i = 0; i < N; i++) std::fill(a[i].begin(), a[i].end(), 0.0);
}

void Matrix::Clean()
{
	a.clear(); N = M = 0;
}

vector<double> Matrix::Gauss(vector<double> b)
{
	if (l2_norm_square(b) < 1e-30) { cout << "<" << endl; return vector<double>(b.size()); }

	for (int i = 0; i < N - 1; i++) {
		// Поиск максимального элемента в столбце
		int maxi = i;
		double max_value = std::abs(a[i][i]);
		for (int j = i + 1; j < N; j++) {
			if (std::abs(a[j][i]) > max_value) {
				maxi = j;
				max_value = std::abs(a[j][i]);
			}
		}

		// Перестановка строк, чтобы максимальный элемент был на диагонали
		if (maxi != i) {
			std::swap(a[i], a[maxi]);
			std::swap(b[i], b[maxi]);
		}

		// Преобразование текущей строки и строк ниже
		for (int j = i + 1; j < N; j++) {
			double factor = a[j][i] / a[i][i];
			for (int k = i; k < M; k++) {
				a[j][k] -= factor * a[i][k];
			}
			b[j] -= factor * b[i];
		}
	}

	// Обратный ход
	std::vector<double> x(N);
	for (int i = N - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < M; j++) {
			sum += a[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / a[i][i];
	}

	return x;
}

vector<double> Matrix::Jacobi(vector<double>& b, const double epsilon)
{
	int max_iter = 2*N;
	int n = N;
	vector<double> x(n);
	vector<double> x_new(n);

	for (int iter = 0; iter < max_iter; ++iter) {
		for (int i = 0; i < n; ++i) {
			double sum = 0.0;
			for (int j = 0; j < n; ++j) {
				if (i != j) {
					sum += a[i][j] * x[j];
				}
			}
			x_new[i] = (b[i] - sum) / a[i][i];
		}

		// Проверка на сходимость
		double norm = 0.0;
		for (int i = 0; i < n; ++i) {
			norm += pow(x_new[i] - x[i], 2);
		}
		norm = sqrt(norm);

		if (norm < epsilon) 
		{
			break;
		}

		x = x_new;
	}
	return x;
}

ostream& operator<<(ostream& cout, const Matrix& b)
{
	for (int i = 0; i < b.N; i++)
	{
		for (int j = 0; j < b.M; j++)
		{
			cout << b.a[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	return cout;
}