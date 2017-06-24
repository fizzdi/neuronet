#include "Layer.h"
#include <string>


using namespace NeuroNet;
Matrix2d Layer::sigm_function(Matrix2d& x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res.at(i, j) = x.at(i, j) <= -35 ? x.at(i, j) = 10e-15 : 1.0 / (1.0 + exp(-x.at(i, j)));
	return std::move(res);
}

Matrix2d Layer::tanh_function(Matrix2d& x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
		{
			//res(i,j) = (exp(2 * x(i,j)) - 1.0) / (exp(2 * x(i,j)) + 1.0);
			res.at(i, j) = std::tanh(x.at(i, j));
		}
	return std::move(res);
}

Matrix2d Layer::diff_tanh_function(Matrix2d& x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res.at(i, j) = 1.0 - x.at(i, j) * x.at(i, j);
	return std::move(res);
}

Matrix2d Layer::diff_sigm_function(Matrix2d& x)
{
	int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
	Matrix2d res(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			res.at(i, j) = (1.0 - x.at(i, j)) * x.at(i, j);
	return std::move(res);
}

Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction)
{
	_aftype = activationFunction;
	Weights = Matrix2d(neuronCount, prevNeuronCount);
	States = Matrix2d(1, neuronCount);
	Axons = Matrix2d(1, neuronCount);
	Delta = Matrix2d(1, neuronCount);
	LastDelta = Matrix2d(1, neuronCount);
	Grad = Matrix2d(neuronCount, prevNeuronCount);
	LastGrad = Matrix2d(neuronCount, prevNeuronCount);
	CorrectVal = Matrix2d(neuronCount, prevNeuronCount);
	DeltaSum = Matrix2d(1, neuronCount);
	GradSum = Matrix2d(neuronCount, prevNeuronCount);
	LastDeltaSum = Matrix2d(1, neuronCount);
	LastGradSum = Matrix2d(neuronCount, prevNeuronCount);
	Bias = Matrix2d(1, neuronCount);
	BiasCorrectVal = Matrix2d(1, neuronCount);
	RMS = Matrix2d(neuronCount, prevNeuronCount);
	RMSBias = Matrix2d(1, neuronCount);
	RMSN = Matrix2d(neuronCount, prevNeuronCount);
	RMSNBias = Matrix2d(1, neuronCount);

	Weights.InitRandom(0, 1);
	CorrectVal.InitRandom(0, 1);
	DeltaSum.Fill(0.0);
	GradSum.Fill(0.0);
	LastDeltaSum.Fill(0.0);
	LastGradSum.Fill(0.0);
	BiasCorrectVal.InitRandom(-1.0, 1.0);
	Bias.InitRandom(-1.0, 1.0);
	States.Fill(0.0);
	RMS.Fill(0.0);
	RMSBias.Fill(0.0);
	RMSN.Fill(0.0);
	RMSNBias.Fill(0.0);
}

void Layer::CalculateStates(Layer & prevLayer)
{
	States = prevLayer.Axons * !Weights + Bias;
}

void Layer::CalculateAxons()
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

Matrix2d Layer::GetDiff()
{
	switch (_aftype)
	{
	case SIGM:
		return diff_sigm_function(Axons);
	case LINE:
	{
		Matrix2d res(Axons.GetVerticalSize(), Axons.GetHorizontalSize());
		res.Fill(1.0);
		return res;
	}
	case TANH:
		return diff_tanh_function(Axons);
	default:
		Matrix2d res(Axons.GetVerticalSize(), Axons.GetHorizontalSize());
		res.Fill(1.0);
		return res;
		break;
	}
}

void Layer::NguenWidrow(double Xmin, double Xmax, double Ymin, double Ymax)
{
	double beta = 0.7*std::pow(Weights.GetVerticalSize(), 1.0 / Weights.GetHorizontalSize());
	for (int i = 0; i < Weights.GetVerticalSize(); ++i)
	{
		double mvij = 0.0;
		for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
			mvij += Weights.at(i, j) * Weights.at(i, j);
		mvij = std::sqrt(mvij);
		if (mvij == 0) mvij = 1;

		for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
		{
			Weights.at(i, j) *= beta / mvij;
		}
		Bias.at(0, i) = getRand(-beta, beta);
	}

	double x = 0.5*(Xmax - Xmin), y = 0.5*(Xmax + Xmin);
	Weights = Weights*x; Bias = Bias*x + y;
	Matrix2d a(1, Weights.GetHorizontalSize()), c(1, Weights.GetHorizontalSize());
	for (int j = 0; j < Weights.GetHorizontalSize(); j++) {
		a.at(0, j) = 2.0 / (Ymax - Ymin);
		c.at(0, j) = 1.0 - Ymax*a.at(0, j);
	}
	Bias = !(Weights * !c + !Bias);

	for (int j = 0; j < Weights.GetVerticalSize(); j++)
		for (int k = 0; k < Weights.GetHorizontalSize(); k++)
			Weights.at(j, k) *= a.at(0, k);
}