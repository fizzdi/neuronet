#pragma once
#include <vector>
namespace NeuroNet
{
	class Matrix2d
	{
	private:
		std::vector<std::vector<double>> _matrix;
	public:
		Matrix2d() {};
		Matrix2d(int n, int m, double val = 0.0);
		void Init(int n, int m, double val = 0.0);
		void InitRandom(int n, int m);
		int GetHorizontalSize() const;
		int GetVerticalSize() const;
		void Clear();
		Matrix2d operator* (const Matrix2d &rhs);
		Matrix2d operator* (const double &rhs);
		Matrix2d operator+ (const Matrix2d &rhs);
		Matrix2d operator- (const Matrix2d &rhs);
		Matrix2d operator+= (const Matrix2d &rhs);
		Matrix2d operator= (const std::vector<std::vector<double>> &rhs);
		Matrix2d operator= (const Matrix2d &rhs);
		Matrix2d operator= (const std::vector<double> &rhs);
		Matrix2d operator! ();
		std::vector<double>& operator[] (const int i);
		Matrix2d multiplication(const Matrix2d &rhs);
		const double sum();
		std::vector<std::vector<double>>::iterator rowbegin();
		std::vector<std::vector<double>>::iterator rowend();
	};
}