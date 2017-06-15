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
		void InitRandom(int n, int m, double minv, double maxv);
		int GetHorizontalSize() const;
		int GetVerticalSize() const;
		void Clear();
		Matrix2d operator* (const Matrix2d &rhs);
		Matrix2d operator* (const double &rhs);
		Matrix2d operator+ (const Matrix2d &rhs);
		Matrix2d operator+ (const double &rhs);
		Matrix2d operator- (const Matrix2d &rhs);
		Matrix2d operator+= (const Matrix2d &rhs);
		Matrix2d operator= (const std::vector<std::vector<double>> &rhs);
		Matrix2d operator= (const Matrix2d &rhs);
		Matrix2d operator= (const std::vector<double> &rhs);
		Matrix2d operator! () const;
		Matrix2d sqrt() const;
		Matrix2d abs() const;
		std::vector<double>& operator[] (const int i);
		Matrix2d multiplication(const Matrix2d &rhs);
		const double sum();
		std::vector<std::vector<double>>::iterator rowbegin();
		std::vector<std::vector<double>>::iterator rowend();

		friend std::ostream& operator<< (std::ostream &os, const Matrix2d &m);
	};
}