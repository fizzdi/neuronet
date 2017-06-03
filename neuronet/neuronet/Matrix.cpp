#include "Matrix.h"

NeuroNet::Matrix::Matrix()
{
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
			_matrix[i][j] = (rand()%101)*1.0/100.0;
		}
	}
}

double NeuroNet::Matrix::get(int i, int j)
{
	return _matrix[i][j];
}

double NeuroNet::Matrix::set(int i, int j, double val)
{
	return _matrix[i][j] = val;
}

double NeuroNet::Matrix::add(int i, int j, double val)
{
	return _matrix[i][j] += val;
}

int NeuroNet::Matrix::GetHorizontalSize()
{
	return _matrix.size() > 0 ? _matrix[0].size() : 0;
}

int NeuroNet::Matrix::GetVerticalSize()
{
	return _matrix.size();
}

void NeuroNet::Matrix::Clear()
{
	int n = GetVerticalSize();
	int m = GetHorizontalSize();
	Init(n, m);
}
