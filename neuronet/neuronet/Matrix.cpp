#include "Matrix.h"

NeuroNet::Matrix::Matrix()
{
}

NeuroNet::Matrix::Matrix(int n, int m)
{
	_matrix.resize(0);
	_matrix.resize(n, std::vector<double>(m));
}

void NeuroNet::Matrix::Init(int n, int m)
{
	_matrix.resize(0);
	_matrix.resize(n, std::vector<double>(m));
}

void NeuroNet::Matrix::InitRandom(int n, int m)
{
	Init(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			_matrix[i][j] = (rand() % 101)*1.0 / 100.0;
		}
	}
}

int NeuroNet::Matrix::size()
{
	return _matrix.size();
}

void NeuroNet::Matrix::Clear()
{
	int n = this->size();
	int m = n > 0 ? (*this)[0].size() : 0;
	Init(n, m);
}

NeuroNet::Matrix NeuroNet::Matrix::operator*(const Matrix & rhs)
{
	int n = this->_matrix.size();
	int m = rhs._matrix.size() > 0 ? rhs._matrix[0].size() : 0;
	Matrix res(n, m);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			for (int k = 0; k < rhs._matrix.size(); ++k)
			{
				res._matrix[i][j] += _matrix[i][k] * rhs._matrix[k][j];
			}
		}
	}
	return res;
}

NeuroNet::Matrix NeuroNet::Matrix::operator+(const Matrix & rhs)
{
	int n = this->_matrix.size();
	int m = rhs._matrix.size() > 0 ? rhs._matrix[0].size() : 0;
	Matrix res(n, m);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			res._matrix[i][j] = _matrix[i][j] + rhs._matrix[i][j];
		}
	}
	return res;
}

NeuroNet::Matrix NeuroNet::Matrix::operator+=(const Matrix & rhs)
{
	for (int i = 0; i < this->_matrix.size(); ++i)
	{
		for (int j = 0; j < this->_matrix[i].size(); ++j)
		{
			this->_matrix[i][j] = _matrix[i][j] + rhs._matrix[i][j];
		}
	}
	return *this;
}

NeuroNet::Matrix NeuroNet::Matrix::operator=(const std::vector<std::vector<double>>& rhs)
{
	_matrix = rhs;
	return *this;
}

std::vector<double>& NeuroNet::Matrix::operator[](const int i)
{
	return _matrix[i];
}
