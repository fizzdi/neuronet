#include "Layer.h"

NeuroNet::Matrix NeuroNet::Layer::sigm_function(Matrix m)
{
	Matrix res(m.size(), m[0].size());
	for (int i = 0; i < m.size(); ++i)
		for (int j = 0; j < m[0].size(); ++j)
			res[i][j] = m[i][j] <= -35 ? m[i][j] = 10e-15 : 1.0 / (1.0 + exp(-m[i][j]));
	return res;
}

NeuroNet::Matrix NeuroNet::Layer::tanh_function(Matrix m)
{
	//TODO get vertical and horizontal sizes
	Matrix res(m.size(), m[0].size());
	for (int i = 0; i < m.size(); ++i)
		for (int j = 0; j < m[0].size(); ++j)
			res[i][j] = (exp(2 * m[i][j]) - 1.0) / (exp(2 * m[i][j]) + 1.0);
	return res;
}

NeuroNet::Matrix NeuroNet::Layer::diff_tanh_function(Matrix m)
{
	Matrix res(m.size(), m[0].size());
	for (int i = 0; i < m.size(); ++i)
		for (int j = 0; j < m[0].size(); ++j)
			res[i][j] = 1.0 - m[i][j] * m[i][j];
	return res;
}

NeuroNet::Matrix NeuroNet::Layer::diff_sigm_function(Matrix m)
{
	Matrix res(m.size(), m[0].size());
	for (int i = 0; i < m.size(); ++i)
		for (int j = 0; j < m[0].size(); ++j)
			res[i][j] = (1.0 - m[i][j]) * m[i][j];
	return res;
}

NeuroNet::Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction)
{
	_aftype = activationFunction;
	Weights.InitRandom(prevNeuronCount, neuronCount);
	Correct.Init(prevNeuronCount, neuronCount);
	States.Init(neuronCount, 1);
	Axons.Init(neuronCount, 1);
	Delta.Init(neuronCount, 1);
}

//TODO add const operation
void NeuroNet::Layer::CalculateStates(Layer & prevLayer)
{
	States = !(!prevLayer.Axons * Weights);
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

NeuroNet::Matrix NeuroNet::Layer::GetDiff()
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


