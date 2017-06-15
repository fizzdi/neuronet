#include "Matrix2d.h"
#include <exception>
#include <algorithm>
#include <iostream>
#include "Common.h"


NeuroNet::Matrix2d::Matrix2d(int n, int m, double val)
{
	_matrix.resize(0);
	_matrix.resize(n, std::vector<double>(m, val));
}

void NeuroNet::Matrix2d::Init(int n, int m, double val)
{
	_matrix.resize(0);
	_matrix.resize(n, std::vector<double>(m, val));
}

void NeuroNet::Matrix2d::InitRandom(int n, int m)
{
	Init(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			_matrix[i][j] = NeuroNet::Common::getRand(-0.5, 0.5);
		}
	}
}

void NeuroNet::Matrix2d::InitRandom(int n, int m, double minv, double maxv)
{
	Init(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			_matrix[i][j] = NeuroNet::Common::getRand(minv, maxv);
		}
	}
}

int NeuroNet::Matrix2d::GetHorizontalSize() const
{
	return (GetVerticalSize() > 0 ? (int)_matrix[0].size() : 0);
}

int NeuroNet::Matrix2d::GetVerticalSize() const
{
	return (int)_matrix.size();
}

void NeuroNet::Matrix2d::Clear()
{
	Init(GetVerticalSize(), GetHorizontalSize());
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator*(const Matrix2d & rhs)
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
				res._matrix[i][j] += _matrix[i][k] * rhs._matrix[k][j];
	return res;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator*(const double & rhs)
{
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res._matrix[i][j] = _matrix[i][j] * rhs;
	return res;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator+(const Matrix2d & rhs)
{
	int n = (int)this->_matrix.size();
	int m = rhs._matrix.size() > 0 ? (int)rhs._matrix[0].size() : 0;
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			res._matrix[i][j] = _matrix[i][j] + rhs._matrix[i][j];
		}
	}
	return res;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator+(const double & rhs)
{
	int n = (int)this->_matrix.size();
	int m = this->_matrix.size() > 0 ? (int)this->_matrix[0].size() : 0;
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			res._matrix[i][j] = _matrix[i][j] + rhs;
		}
	}
	return res;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator-(const Matrix2d & rhs)
{
	int n = (int)this->_matrix.size();
	int m = rhs._matrix.size() > 0 ? (int)rhs._matrix[0].size() : 0;
	Matrix2d res(n, m);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			res._matrix[i][j] = _matrix[i][j] - rhs._matrix[i][j];
		}
	}
	return res;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator+=(const Matrix2d & rhs)
{
	for (int i = 0; i < this->_matrix.size(); ++i)
	{
		for (int j = 0; j < this->_matrix[i].size(); ++j)
		{
			_matrix[i][j] += rhs._matrix[i][j];
		}
	}
	return *this;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator=(const std::vector<std::vector<double>>& rhs)
{
	_matrix.resize(rhs.size());
	for (int i = 0; i < rhs.size(); ++i)
	{
		_matrix[i].resize(rhs[i].size());
		for (int j = 0; j < rhs[i].size(); ++j)
			_matrix[i][j] = rhs[i][j];
	}
	return *this;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator=(const Matrix2d & rhs)
{
	_matrix.resize(rhs._matrix.size());
	for (int i = 0; i < this->_matrix.size(); ++i)
	{
		_matrix[i].resize(rhs._matrix[i].size());
		for (int j = 0; j < this->_matrix[i].size(); ++j)
		{
			this->_matrix[i][j] = rhs._matrix[i][j];
		}
	}
	return *this;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator!() const
{
	Matrix2d old;
	int n = GetHorizontalSize();
	int m = GetVerticalSize();
	old.Init(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			old._matrix[j][i] = _matrix[i][j];
	return old;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::sqrt() const
{
	Matrix2d old;
	int n = GetHorizontalSize();
	int m = GetVerticalSize();
	old.Init(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			old._matrix[j][i] = std::sqrt(_matrix[i][j]);
	return old;
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::abs() const
{
	Matrix2d old;
	int n = GetHorizontalSize();
	int m = GetVerticalSize();
	old.Init(n, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			old._matrix[j][i] = std::abs(_matrix[i][j]);
	return old;
}

std::vector<double>& NeuroNet::Matrix2d::operator[](const int i)
{
	return _matrix[i];
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::multiplication(const Matrix2d & rhs)
{
	Matrix2d res((int)rhs._matrix.size(), (int)rhs._matrix[0].size());
	for (int i = 0; i < rhs._matrix.size(); ++i)
		for (int j = 0; j < rhs._matrix[i].size(); ++j)
			res[i][j] = this->_matrix[i][j] * rhs._matrix[i][j];
	return res;
}

const double NeuroNet::Matrix2d::sum()
{
	double ans = 0;
	for (int i = 0; i < _matrix.size(); ++i)
		for (int j = 0; j < _matrix[i].size(); ++j)
			ans += _matrix[i][j];
	return ans;
}

std::vector<std::vector<double>>::iterator NeuroNet::Matrix2d::rowbegin()
{
	return _matrix.begin();
}

std::vector<std::vector<double>>::iterator NeuroNet::Matrix2d::rowend()
{
	return _matrix.end();
}

NeuroNet::Matrix2d NeuroNet::Matrix2d::operator=(const std::vector<double>& rhs)
{
	this->_matrix.resize(1);
	_matrix[0].resize(rhs.size());
	std::copy(rhs.begin(), rhs.end(), _matrix[0].begin());
	return *this;
}

std::ostream & NeuroNet::operator<<(std::ostream & os, const Matrix2d & m)
{
	for (int j = 0; j < m.GetVerticalSize(); ++j)
	{
		for (int k = 0; k < m.GetHorizontalSize(); ++k)
			os << m._matrix[j][k] << " ";
		os << std::endl;
	}

	return os;
}
