#pragma once
#include <vector>
#include <fstream>
#include <random>
#define EPS (1e-4)
#define ERRORDEF "###ERROR"
extern std::ofstream debug;

extern std::random_device rnd_dev;
extern std::mt19937 eng;
extern std::uniform_int_distribution<> rnd_gen;

namespace NeuroNet
{
	double getRand(double vmin, double vmax, bool integer = false);
	int myrand();

	class Matrix2d
	{
	private:
		int N, M;
	public:
		double* matrix;
		Matrix2d() { N = M = 0; matrix = nullptr; };
		~Matrix2d();
		Matrix2d(const Matrix2d& rhs);
		Matrix2d(int n, int m);
		Matrix2d(const std::vector<double> &rhs);
		void Fill(double val = 0.0);
		void InitRandom(double minv = -1.0, double maxv = 1.0);

		int GetHorizontalSize() const;
		int GetVerticalSize() const;

		Matrix2d operator! () const;
		Matrix2d operator- () const;

		Matrix2d& operator+= (const Matrix2d &rhs);
		Matrix2d operator+ (const Matrix2d &rhs) const;
		Matrix2d& operator+= (const double rhs);
		Matrix2d operator+ (const double rhs) const;

		Matrix2d& operator-= (const Matrix2d &rhs);
		Matrix2d operator- (const Matrix2d &rhs) const;
		Matrix2d& operator-= (const double rhs);
		Matrix2d operator- (const double rhs) const;

		Matrix2d& operator*= (const double rhs);
		Matrix2d operator* (const Matrix2d &rhs);
		Matrix2d operator* (const double &rhs);

		Matrix2d& operator= (const std::vector<std::vector<double>> &rhs);
		Matrix2d& operator= (const Matrix2d &rhs);
		Matrix2d& operator= (const std::vector<double> &rhs);

		bool operator== (const Matrix2d &rhs);
		friend Matrix2d operator*(const double lhs, const Matrix2d& rhs);
		friend Matrix2d operator/(const double lhs, const Matrix2d& rhs);


		Matrix2d abs() const;
		Matrix2d multiplication(const Matrix2d &rhs) const;
		const double sum() const;

		double& at(const int& i, const int& j);
		//double& at(const int i, const int j);
		double at(const int& i, const int& j) const;
		//double at(const int i, const int j) const;

		friend std::ostream& operator<< (std::ostream &os, const Matrix2d &m);
		friend Matrix2d sqrt(const Matrix2d&rhs);
	};

	struct Problem
	{
		Matrix2d inputs, outputs;
		Problem(std::vector<double> input, std::vector<double> output)
		{
			inputs = input;
			outputs = output;
		}
		Problem()
		{};
	};
}