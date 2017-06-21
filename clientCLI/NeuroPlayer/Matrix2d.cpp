#include "Matrix2d.h"
#include <string>

using namespace NeuroNet;
Matrix2d::~Matrix2d()
{
	if (n != 0)
		delete[](_m);
	_m = nullptr;
}

Matrix2d::Matrix2d(const Matrix2d & rhs)
{
	if (*this == rhs) return;
	this->n = rhs.n;
	this->m = rhs.m;
	_m = new double[n*m];
	memcpy_s(_m, n*m * sizeof(*_m), rhs._m, n*m * sizeof(*rhs._m));
}

Matrix2d::Matrix2d(Matrix2d && rhs)
{
/*	if (n > 0 && n < 1500000)
	{
		delete[] _m;
	}*/
	_m = std::move(rhs._m);
	n = rhs.n;
	m = rhs.m;
	rhs._m = nullptr;
}

Matrix2d::Matrix2d(int n, int m)
{
	this->n = n;
	this->m = m;
	_m = new double[n*m];
}

Matrix2d::Matrix2d(const std::vector<double>& rhs)
{
	this->n = 1;
	this->m = rhs.size();
	_m = new double[m];
	memcpy_s(_m, n*m * sizeof(*_m), rhs.data(), n*m * sizeof(*rhs.data()));
}

void Matrix2d::Fill(double val)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			at(i, j) = val;
}

void Matrix2d::InitRandom(double minv, double maxv)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			at(i, j) = getRand(minv, maxv);
}

int Matrix2d::GetHorizontalSize() const
{
	return m;
}

int Matrix2d::GetVerticalSize() const
{
	return n;
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

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
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
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
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

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
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

Matrix2d Matrix2d::operator*(const Matrix2d & rhs) const
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
	res.Fill(0.0);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			for (int k = 0; k < nm; ++k)
				res.at(i, j) += this->at(i, k) * rhs.at(k, j);
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
	//if (_m != nullptr)
	delete[](_m);
	this->n = rhs.size();
	this->m = n > 0 ? rhs[0].size() : 0;
	_m = new double[n*m];

	for (int i = 0; i < (int)rhs.size(); ++i)
		memcpy_s(_m + (i*m) * sizeof(*_m), m * sizeof(*_m), rhs.data(), m * sizeof(*rhs.data()));
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

	//if (_m != nullptr)
	delete[](_m);


	this->n = rhs.n;
	this->m = rhs.m;
	_m = new double[n*m];

	memcpy_s(_m, n*m * sizeof(*_m), rhs._m, n*m * sizeof(*rhs._m));
	return *this;
}

Matrix2d & Matrix2d::operator=(Matrix2d && rhs)
{
	if (*this == rhs) return *this;
	delete[](_m);
	_m = std::move(rhs._m);
	n = rhs.n;
	m = rhs.m;

	rhs._m = nullptr;
	return *this;
}

Matrix2d& Matrix2d::operator=(const std::vector<double>& rhs)
{
	//if (_m != nullptr)
	delete[](_m);
	this->n = 1;
	this->m = rhs.size();
	_m = new double[n*m];
	memcpy_s(_m, n*m * sizeof(*_m), rhs.data(), n*m * sizeof(*rhs.data()));
	return *this;
}

bool Matrix2d::operator==(const Matrix2d & rhs)
{
	return n == rhs.n && m == rhs.m && _m == rhs._m;
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

	Matrix2d res(rhs.n, rhs.m);
	for (int i = 0; i < rhs.n; ++i)
		for (int j = 0; j < rhs.m; ++j)
			res.at(i, j) = this->at(i, j) * rhs.at(i, j);
	return  std::move(res);
}

const double Matrix2d::sum() const
{
	double ans = 0;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			ans += at(i, j);
	return ans;
}

double & Matrix2d::at(const int i, const int j)
{
	if (i >= n || j >= m)
	{
		debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		throw std::logic_error("Wrong sizes in at operation");
	}
	return _m[i*m + j];
}

double Matrix2d::at(const int i, const int j) const
{
	if (i >= n || j >= m)
	{
		debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		throw std::logic_error("Wrong sizes in at operation");
	}
	return _m[i*m + j];
}

double NeuroNet::getRand(double vmin, double vmax, bool integer)
{
	if (integer)
		return (int)(vmin + rand() * (vmax - vmin) / RAND_MAX);
	return vmin + rand() * (vmax - vmin) / RAND_MAX;//
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

