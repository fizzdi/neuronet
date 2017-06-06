#include <vector>
#pragma once
namespace NeuroNet
{
	class Matrix
	{
	private:
		std::vector<std::vector<double>> _matrix;
	public:
		Matrix();
		Matrix(int n, int m);
		void Init(int n, int m);
		void InitRandom(int n, int m);
		const int size();
		void Clear();
		Matrix operator* (const Matrix &rhs);
		Matrix operator+ (const Matrix &rhs);
		Matrix operator- (const Matrix &rhs);
		Matrix operator+= (const Matrix &rhs);
		Matrix operator= (const std::vector<std::vector<double>> &rhs);
		Matrix operator= (const Matrix &rhs);
		Matrix operator= (const std::vector<double> &rhs);
		Matrix operator! ();
		std::vector<double>& operator[] (const int i);
		Matrix multiplication(const Matrix &rhs);
		const double sum();
	};
}