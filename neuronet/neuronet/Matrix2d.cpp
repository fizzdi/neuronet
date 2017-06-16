#include "Matrix2d.h"
#include <exception>
#include <algorithm>
#include <iostream>
#include "Common.h"


using namespace NeuroNet;

Matrix2d::Matrix2d(int n, int m, double val)
{
	_m.resize(0);
	_m.resize(n, std::vector<double>(m, val));
}

void Matrix2d::Init(int n, int m, double val)
{
	_m.resize(0);
	_m.resize(n, std::vector<double>(m, val));
}

void Matrix2d::InitRandom(int n, int m, double minv, double maxv)
{
	Init(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			_m[i][j] = Common::getRand(minv, maxv);
		}
	}
}

void Matrix2d::Clear()
{
	Init(GetVerticalSize(), GetHorizontalSize());
}

int Matrix2d::GetHorizontalSize() const
{
	return (GetVerticalSize() > 0 ? (int)_m[0].size() : 0);
}

int Matrix2d::GetVerticalSize() const
{
	return (int)_m.size();
}

Matrix2d Matrix2d::operator!() const
{
	Matrix2d res;
	int m = GetVerticalSize();
	int n = GetHorizontalSize();
	res.Init(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			res._m[j][i] = _m[i][j];
	return res;
}

Matrix2d Matrix2d::operator-() const
{
	Matrix2d res = *this;
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res._m[i][j] = -_m[i][j];
	return res;
}

Matrix2d Matrix2d::operator+=(const Matrix2d & rhs)
{
	if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
		throw std::logic_error("Wrong sizes in addition operation");

	for (int i = 0; i < this->_m.size(); ++i)
		for (int j = 0; j < this->_m[i].size(); ++j)
			_m[i][j] += rhs._m[i][j];
	return *this;
}

Matrix2d Matrix2d::operator+(const Matrix2d & rhs) const
{
	Matrix2d res = *this;
	return res += rhs;
}

Matrix2d NeuroNet::Matrix2d::operator+=(const double rhs)
{
	for (int i = 0; i < this->_m.size(); ++i)
		for (int j = 0; j < this->_m[i].size(); ++j)
			_m[i][j] += rhs;
	return *this;
}

Matrix2d Matrix2d::operator+(const double rhs) const
{
	Matrix2d res = *this;
	return res += rhs;
}

Matrix2d Matrix2d::operator-=(const Matrix2d & rhs)
{
	if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
		throw std::logic_error("Wrong sizes in subtraction operation");

	for (int i = 0; i < this->_m.size(); ++i)
		for (int j = 0; j < this->_m[i].size(); ++j)
			_m[i][j] -= rhs._m[i][j];
	return *this;
}

Matrix2d Matrix2d::operator-(const Matrix2d & rhs) const
{
	Matrix2d res = *this;
	return res -= rhs;
}

Matrix2d NeuroNet::Matrix2d::operator-=(const double rhs)
{
	return *this += -rhs;
}

Matrix2d NeuroNet::Matrix2d::operator-(const double rhs) const
{
	Matrix2d res = *this;
	return res += -rhs;
}

Matrix2d Matrix2d::operator*(const Matrix2d & rhs) const
{
	if (GetHorizontalSize() != rhs.GetVerticalSize())
		throw std::logic_error("Wrong sizes in multiplication operation");

	int n = GetVerticalSize();
	int m = rhs.GetHorizontalSize();
	int nm = GetHorizontalSize();
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			for (int k = 0; k < nm; ++k)
				res._m[i][j] += _m[i][k] * rhs._m[k][j];
	return res;
}

Matrix2d Matrix2d::operator*(const double & rhs)
{
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res._m[i][j] = _m[i][j] * rhs;
	return res;
}

Matrix2d Matrix2d::operator=(const std::vector<std::vector<double>>& rhs)
{
	_m.resize(rhs.size());
	for (int i = 0; i < rhs.size(); ++i)
	{
		_m[i].resize(rhs[i].size());
		for (int j = 0; j < rhs[i].size(); ++j)
			_m[i][j] = rhs[i][j];
	}
	return *this;
}

Matrix2d Matrix2d::operator=(const Matrix2d & rhs)
{
	_m.resize(rhs._m.size());
	for (int i = 0; i < this->_m.size(); ++i)
	{
		_m[i].resize(rhs._m[i].size());
		for (int j = 0; j < this->_m[i].size(); ++j)
		{
			this->_m[i][j] = rhs._m[i][j];
		}
	}
	return *this;
}

Matrix2d Matrix2d::operator=(const std::vector<double>& rhs)
{
	this->_m.resize(1);
	_m[0].resize(rhs.size());
	std::copy(rhs.begin(), rhs.end(), _m[0].begin());
	return *this;
}

Matrix2d Matrix2d::abs() const
{
	Matrix2d old;
	int n = GetHorizontalSize();
	int m = GetVerticalSize();
	old.Init(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			old._m[j][i] = std::abs(_m[i][j]);
	return old;
}


Matrix2d Matrix2d::multiplication(const Matrix2d & rhs) const
{
	if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
		throw std::logic_error("Wrong sizes in multiplication operation");

	Matrix2d res((int)rhs._m.size(), (int)rhs._m[0].size());
	for (int i = 0; i < rhs._m.size(); ++i)
		for (int j = 0; j < rhs._m[i].size(); ++j)
			res._m[i][j] = this->_m[i][j] * rhs._m[i][j];
	return res;
}

const double Matrix2d::sum() const
{
	double ans = 0;
	for (int i = 0; i < _m.size(); ++i)
		for (int j = 0; j < _m[i].size(); ++j)
			ans += _m[i][j];
	return ans;
}

double & Matrix2d::operator()(const int i, const int j)
{
	return _m[i][j];
}


std::ostream & NeuroNet::operator<<(std::ostream & os, const Matrix2d & m)
{
	for (int j = 0; j < m.GetVerticalSize(); ++j)
	{
		for (int k = 0; k < m.GetHorizontalSize(); ++k)
			os << m._m[j][k] << " ";
		os << std::endl;
	}

	return os;
}

Matrix2d NeuroNet::sqrt(const Matrix2d & rhs)
{
	Matrix2d res;
	int n = rhs.GetHorizontalSize();
	int m = rhs.GetVerticalSize();
	res.Init(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			res._m[j][i] = std::abs(rhs._m[i][j]);
	return res;
}
