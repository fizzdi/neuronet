#include "Matrix2d.h"
#include <string>

using namespace NeuroNet;

std::random_device rnd_dev;
std::mt19937 eng = std::mt19937(rnd_dev());
std::uniform_int_distribution<> rnd_gen = std::uniform_int_distribution<>(0, RAND_MAX);

Matrix2d::~Matrix2d()
{
	if (N != 0)
		delete[](matrix);
	matrix = nullptr;
}

Matrix2d::Matrix2d(const Matrix2d & rhs)
{
	if (*this == rhs) return;
	this->N = rhs.N;
	this->M = rhs.M;
	matrix = new double[N*M];
	memcpy_s(matrix, N*M * sizeof(*matrix), rhs.matrix, N*M * sizeof(*rhs.matrix));
}

Matrix2d::Matrix2d(int n, int m)
{
	this->N = n;
	this->M = m;
	matrix = new double[n*m];
}

Matrix2d::Matrix2d(const std::vector<double>& rhs)
{
	this->N = 1;
	this->M = rhs.size();
	matrix = new double[M];
	memcpy_s(matrix, N*M * sizeof(*matrix), rhs.data(), N*M * sizeof(*rhs.data()));
}

void Matrix2d::Fill(double val)
{
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			at(i, j) = val;
}

void Matrix2d::InitRandom(double minv, double maxv)
{
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			at(i, j) = getRand(minv, maxv);
}

int Matrix2d::GetHorizontalSize() const
{
	return M;
}

int Matrix2d::GetVerticalSize() const
{
	return N;
}

Matrix2d Matrix2d::operator!() const
{
	int m = GetVerticalSize();
	int n = GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			res.at(j, i) = this->at(i, j);
	return std::move(res);
}

Matrix2d Matrix2d::operator-() const
{
	Matrix2d res = *this;
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res.at(i, j) = -at(i, j);
	return std::move(res);
}

Matrix2d& Matrix2d::operator+=(const Matrix2d & rhs)
{
	if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
	{
		debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		throw std::logic_error("Wrong sizes in addition operation");
	}

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			at(i, j) += rhs.at(i, j);
	return *this;
}

Matrix2d Matrix2d::operator+(const Matrix2d & rhs) const
{
	Matrix2d res = *this;
	return std::move(res += rhs);
}

Matrix2d& Matrix2d::operator+=(const double rhs)
{
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			at(i, j) += rhs;
	return *this;
}

Matrix2d Matrix2d::operator+(const double rhs) const
{
	Matrix2d res = *this;
	return std::move(res += rhs);
}

Matrix2d& Matrix2d::operator-=(const Matrix2d & rhs)
{
	if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
	{
		debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		throw std::logic_error("Wrong sizes in subtraction operation");
	}

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			at(i, j) -= rhs.at(i, j);
	return *this;
}

Matrix2d Matrix2d::operator-(const Matrix2d & rhs) const
{
	Matrix2d res = *this;
	return std::move(res -= rhs);
}

Matrix2d& Matrix2d::operator-=(const double rhs)
{
	return *this += -rhs;
}

Matrix2d Matrix2d::operator-(const double rhs) const
{
	Matrix2d res = *this;
	return std::move(res += -rhs);
}

Matrix2d & NeuroNet::Matrix2d::operator*=(const double rhs)
{
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			at(i, j) *= rhs;
	return *this;
}

Matrix2d Matrix2d::operator*(const Matrix2d & rhs)
{
	if (GetHorizontalSize() != rhs.GetVerticalSize())
	{
		debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		throw std::logic_error("Wrong sizes in multiplication operation");
	}

	int n = GetVerticalSize();
	int m = rhs.GetHorizontalSize();
	int nm = GetHorizontalSize();
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
		{
			res.at(i, j) = 0.0;
			for (int k = 0; k < nm; ++k)
				res.at(i, j) += this->at(i, k) * rhs.at(k, j);
		}
	return std::move(res);
}

Matrix2d Matrix2d::operator*(const double & rhs)
{
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res.at(i, j) = at(i, j) * rhs;
	return std::move(res);
}

Matrix2d& Matrix2d::operator=(const std::vector<std::vector<double>>& rhs)
{
	//if (matrix != nullptr)
	delete[](matrix);
	this->N = rhs.size();
	this->M = N > 0 ? rhs[0].size() : 0;
	matrix = new double[N*M];

	for (int i = 0; i < (int)rhs.size(); ++i)
		memcpy_s(matrix + i*M, M * sizeof(*matrix), rhs[i].data(), M * sizeof(*rhs.data()));
	//for (int j = 0; j < m; ++j)
	//	at(i,j) = rhs[i][j];
	return *this;
}

Matrix2d& Matrix2d::operator=(const Matrix2d & rhs)
{
	if (*this == rhs)
	{
		return *this;
	}

	//if (matrix != nullptr)
	delete[](matrix);


	this->N = rhs.N;
	this->M = rhs.M;
	matrix = new double[N*M];

	memcpy_s(matrix, N*M * sizeof(*matrix), rhs.matrix, N*M * sizeof(*rhs.matrix));
	return *this;
}

Matrix2d& Matrix2d::operator=(const std::vector<double>& rhs)
{
	//if (matrix != nullptr)
	delete[](matrix);
	this->N = 1;
	this->M = rhs.size();
	matrix = new double[N*M];
	memcpy_s(matrix, N*M * sizeof(*matrix), rhs.data(), N*M * sizeof(*rhs.data()));
	return *this;
}

bool Matrix2d::operator==(const Matrix2d & rhs)
{
	return N == rhs.N && M == rhs.M && matrix == rhs.matrix;
}

Matrix2d Matrix2d::abs() const
{
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	Matrix2d old(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			old.at(i, i) = std::abs(at(i, j));
	return std::move(old);
}

Matrix2d Matrix2d::multiplication(const Matrix2d & rhs) const
{
	if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
	{
		debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		throw std::logic_error("Wrong sizes in multiplication operation");
	}

	Matrix2d res(rhs.N, rhs.M);
	for (int i = 0; i < rhs.N; ++i)
		for (int j = 0; j < rhs.M; ++j)
			res.at(i, j) = this->at(i, j) * rhs.at(i, j);
	return  std::move(res);
}

const double Matrix2d::sum() const
{
	double ans = 0;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
			ans += at(i, j);
	return ans;
}

void NeuroNet::Matrix2d::copy(int pos_source, int pos_dest, int len, const Matrix2d & source)
{
	int sz = len * sizeof(*matrix);
	memcpy_s(matrix + pos_dest, sz, source.matrix + pos_source, sz);
}

double & Matrix2d::at(const int &i, const int &j)
{
	return matrix[i*M + j];
}

double Matrix2d::at(const int& i, const int& j) const
{
	return matrix[i*M + j];
}

double NeuroNet::getRand(double vmin, double vmax, bool integer)
{
	if (integer)
		return (int)(vmin + myrand() * (vmax - vmin) / RAND_MAX);
	return vmin + myrand() * (vmax - vmin) / RAND_MAX;//
}

int NeuroNet::myrand()
{
	return rnd_gen(eng);
}

Matrix2d NeuroNet::operator*(const double lhs, const Matrix2d& rhs)
{
	Matrix2d res(rhs.N, rhs.M);
	for (int i = 0; i < rhs.N; ++i)
		for (int j = 0; j < rhs.M; ++j)
			res.at(i, j) = lhs * rhs.at(i, j);
	return  std::move(res);
}

Matrix2d NeuroNet::operator/(const double lhs, const Matrix2d & rhs)
{
	Matrix2d res(rhs.N, rhs.M);
	for (int i = 0; i < rhs.N; ++i)
		for (int j = 0; j < rhs.M; ++j)
			res.at(i, j) = lhs / rhs.at(i, j);
	return std::move(res);
}

std::ostream & NeuroNet::operator<<(std::ostream & os, const Matrix2d & m)
{
	for (int j = 0; j < m.GetVerticalSize(); ++j)
	{
		for (int k = 0; k < m.GetHorizontalSize(); ++k)
			os << m.at(j, k) << " ";
		os << std::endl;
	}

	return os;
}

std::istream & NeuroNet::operator >> (std::istream & os, Matrix2d & m)
{
	for (int j = 0; j < m.GetVerticalSize(); ++j)
		for (int k = 0; k < m.GetHorizontalSize(); ++k)
			os >> m.at(j, k);

	return os;
}

Matrix2d NeuroNet::sqrt(const Matrix2d & rhs)
{
	int n = rhs.GetVerticalSize();
	int m = rhs.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res.at(i, j) = std::sqrt(rhs.at(i, j));
	return std::move(res);
}

