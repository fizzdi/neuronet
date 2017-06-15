#include "Layer.h"
#include <sstream>
#include <iostream>
#include "Common.h"

NeuroNet::Matrix2d NeuroNet::Layer::sigm_function(Matrix2d x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res[i][j] = x[i][j] <= -35 ? x[i][j] = 10e-15 : 1.0 / (1.0 + exp(-x[i][j]));
	return res;
}

NeuroNet::Matrix2d NeuroNet::Layer::tanh_function(Matrix2d x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
		{
			//res[i][j] = (exp(2 * x[i][j]) - 1.0) / (exp(2 * x[i][j]) + 1.0);
			res[i][j] = std::tanh(x[i][j]);
			auto ex = exp(2 * x[i][j]);
			auto ttt = std::tanh(x[i][j]);
			std::ostringstream ss;
			ss << res[i][j];
			auto a = ss.str();
			ss.str("");
			//ss.clear();
			ss << -nan("");
			auto b = ss.str();
			if (a == b)
			{
				int y = 0;
				std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << "!!!!!!!!!!!!ALARM" << std::endl;
				system("pause");
			}
		}
	return res;
}

NeuroNet::Matrix2d NeuroNet::Layer::diff_tanh_function(Matrix2d x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res[i][j] = 1.0 - x[i][j] * x[i][j];
	return res;
}

NeuroNet::Matrix2d NeuroNet::Layer::diff_sigm_function(Matrix2d x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res[i][j] = (1.0 - x[i][j]) * x[i][j];
	return res;
}

NeuroNet::Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction, bool bias)
{
	_aftype = activationFunction;
	Weights.InitRandom(neuronCount, prevNeuronCount, -1, 1);
	States.Init(1, neuronCount);
	Axons.Init(1, neuronCount);
	Delta.Init(1, neuronCount);
	LastDelta.Init(1, neuronCount);
	Grad.Init(neuronCount, prevNeuronCount);
	LastGrad.Init(neuronCount, prevNeuronCount);
	CorrectVal.InitRandom(neuronCount, prevNeuronCount,0,1);
	if (bias)
		Bias.InitRandom(1, neuronCount);
	else Bias.Init(1, neuronCount);
	//BiasCorrectVal.Init(1, neuronCount, 0.1);
	BiasCorrectVal.InitRandom(1, neuronCount);
}

//TODO add const operation
void NeuroNet::Layer::CalculateStates(Layer & prevLayer)
{
	States = prevLayer.Axons * !Weights + Bias;
}

void NeuroNet::Layer::CalculateAxons()
{
	switch (_aftype)
	{
	case SIGM:
		Axons = sigm_function(States);
		break;
	case LINE:
		Axons = States;
		break;
	case TANH:
		Axons = tanh_function(States);
		break;
	default:
		break;
	}
}

NeuroNet::Matrix2d NeuroNet::Layer::GetDiff()
{
	switch (_aftype)
	{
	case SIGM:
		return diff_sigm_function(Axons);
	case LINE:
		return Matrix2d(Axons.GetVerticalSize(), Axons.GetHorizontalSize(), 1.0);
	case TANH:
		return diff_tanh_function(Axons);
	default:
		break;
	}
}

double NeuroNet::Layer::GetDiff(double val)
{
	switch (_aftype)
	{
	case SIGM:
		return (1.0 - val) * val;
	case LINE:
		return val;
	case TANH:
		return 1.0 - val * val;
	default:
		break;
	}
}

void NeuroNet::Layer::NguenWidrow(double Xmin, double Xmax, double Ymin, double Ymax)
{
	double beta = 0.7*std::pow(Weights.GetVerticalSize(), 1.0 / Weights.GetHorizontalSize());
	for (int i = 0; i < Weights.GetVerticalSize(); ++i)
	{
		double mvij = 0.0;
		for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
			mvij += Weights[i][j] * Weights[i][j];
		mvij = std::sqrt(mvij);
		if (mvij == 0) mvij = 1;

		for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
		{
			Weights[i][j] *= beta / mvij;
		}
		Bias[0][i] = NeuroNet::Common::getRand(-beta, beta);
	}

	double x = 0.5*(Xmax - Xmin), y = 0.5*(Xmax + Xmin);
	Weights = Weights*x; Bias = Bias*x + y;
	Matrix2d a(1, Weights.GetHorizontalSize()), c(1, Weights.GetHorizontalSize());
	for (int j = 0; j < Weights.GetHorizontalSize(); j++) {
		a[0][j] = 2.0 / (Ymax - Ymin);
		c[0][j] = 1.0 - Ymax*a[0][j];
	}
	Bias = !(Weights*!c + !Bias);

	for (int j = 0; j < Weights.GetVerticalSize(); j++)
		for (int k = 0; k < Weights.GetHorizontalSize(); k++)
			Weights[j][k]*=a[0][k];
}


double operator* (const std::vector<double> &lhs, const std::vector<double> &rhs)
{
	double res = 0.0;
	for (int i = 0; i < lhs.size(); ++i)
	{
		res += lhs[i] * rhs[i];
	}
	return res;
}


