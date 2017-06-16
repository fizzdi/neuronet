#pragma once
#include <vector>
namespace NeuroNet
{
	class Matrix2d
	{
	private:
		std::vector<std::vector<double>> _m;
	public:
		Matrix2d() {};
		Matrix2d(int n, int m, double val = 0.0);
		void Init(int n, int m, double val = 0.0);
		void InitRandom(int n, int m, double minv = -1.0, double maxv = 1.0);

		void Clear();
		int GetHorizontalSize() const;
		int GetVerticalSize() const;
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
		Matrix2d operator- () const;
		Matrix2d abs() const;
		Matrix2d multiplication(const Matrix2d &rhs);
		const double sum();
		std::vector<std::vector<double>>::iterator rowbegin();
		std::vector<std::vector<double>>::iterator rowend();
		double& operator() (const int i, const int j);

		friend std::ostream& operator<< (std::ostream &os, const Matrix2d &m);
		friend Matrix2d sqrt(const Matrix2d&rhs);
	};
}