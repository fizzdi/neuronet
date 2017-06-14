#include "Layer.h"

NeuroNet::Matrix2d NeuroNet::Layer::sigm_function(Matrix2d x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n,m);
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
			res[i][j] = (exp(2 *x[i][j]) - 1.0) / (exp(2 * x[i][j]) + 1.0);
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

NeuroNet::Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction)
{
	_aftype = activationFunction;
	Weights.InitRandom(neuronCount, prevNeuronCount);
	Correct.Init(neuronCount, prevNeuronCount);
	States.Init(neuronCount,1);
	Axons.Init (neuronCount,1);
	Delta.Init (neuronCount,1);
	LastDelta.Init(neuronCount,1);
}

//TODO add const operation
void NeuroNet::Layer::CalculateStates(Layer & prevLayer)
{
	States =  prevLayer.Axons * !Weights;
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
		return Axons;
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


double operator* (const std::vector<double> &lhs, const std::vector<double> &rhs)
{
	double res = 0.0;
	for (int i = 0; i < lhs.size(); ++i)
	{
		res += lhs[i] * rhs[i];
	}
	return res;
}


