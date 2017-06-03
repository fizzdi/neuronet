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
		void Init(int n, int m);
		void InitRandom(int n, int m);
		double get(int i, int j);
		double set(int i, int j, double val);
		double add(int i, int j, double val);
		int GetHorizontalSize();
		int GetVerticalSize();
		void Clear();
	};
}