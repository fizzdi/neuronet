#include "NeuralNetwork.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Common.h"

void NeuroNet::NeuralNetwork::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount, HiddenLayerFunction, true));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	_layers.push_back(Layer(OutputCount, NeuronCount, LINE, true));
	//_layers.back().NguenWidrow(-1, 1, -1, 1);

	_learn = Learn;

	LastDeltaSum.resize(3);
	LastGradSum.resize(3);
	for (int i = 0; i < 3; ++i)
	{
		LastGradSum[i].Init(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
		LastDeltaSum[i].Init(_layers[i].Delta.GetVerticalSize(), _layers[i].Delta.GetHorizontalSize());
	}
}

//TODO global var debug
double NeuroNet::NeuralNetwork::RunTrainingSet(bool print)
{
	double error = 0.0;
	for each (Problem test in TrainingSet)
	{
		Run(test);
		if (print) PrintProblemResult(test);
		error = CalculateError(test);
		CalcCorrectWeights(test);
	}
	return error;
}

double NeuroNet::NeuralNetwork::RunTrainingSetOffline(bool print)
{
	int countLayers = (int)_layers.size();
	std::vector<Matrix2d> GradSum(3);
	std::vector<Matrix2d> DeltaSum(3);
	
	for (int i = 0; i < 3; ++i)
	{
		GradSum[i] = Matrix2d(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
		DeltaSum[i] = Matrix2d(_layers[i].Delta.GetVerticalSize(), _layers[i].Delta.GetHorizontalSize());
	}
	int y = 0;
	for each (Problem test in TrainingSet)
	{
		Run(test);

		for (int i = 1; i < 3; ++i)
		{
			GradSum[i] += _layers[i].Grad;
			DeltaSum[i] += _layers[i].Delta;
		}
		if (print) PrintProblemResult(test);
	}
	if (print) std::cout << std::endl << "=======CORRECT==========" << std::endl;
	/*std::cout << "GradSum" << std::endl;
	std::cout << GradSum[1] << std::endl;
	std::cout << GradSum[2] << std::endl;
	std::cout << "DeltaSum" << std::endl;
	std::cout << DeltaSum[1] << std::endl;
	std::cout << DeltaSum[2] << std::endl;
	std::cout << "CorrectVal" << std::endl;
	std::cout << _layers[1].CorrectVal << std::endl;
	std::cout << _layers[2].CorrectVal << std::endl;
	std::cout << "BiasCorrectVal" << std::endl;
	std::cout << _layers[1].BiasCorrectVal << std::endl;
	std::cout << _layers[2].BiasCorrectVal << std::endl;*/

	double normGrad = 0.0;
	for (int i = 0; i < countLayers; ++i)
		normGrad += sqrt(GradSum[i] * !GradSum[i]).sum();

	if (normGrad < 1e-6)
		return 0.0;
	ResilientPropagation(LastGradSum, GradSum, LastDeltaSum, DeltaSum);
	LastGradSum = GradSum;
	LastDeltaSum = DeltaSum;


	double error = 0.0;
	for each (Problem test in TrainingSet)
	{
		Run(test);
		if (print) PrintProblemResult(test);
		error += CalculateError(test, print);
	}
	return error;
}

void NeuroNet::NeuralNetwork::PrintProblemResult(Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.GetHorizontalSize(); ++i)
		std::cout << test.outputs(0,i) << " ";
	std::cout << std::endl;
}

double NeuroNet::NeuralNetwork::CalculateError(Problem & test, bool print)
{
	//MSE
	double error = (test.outputs - _layers.back().Axons).multiplication(test.outputs - _layers.back().Axons).sum() / 2.0;
	//RootMSE
	//double error = (test.outputs - _layers.back().Axons).abs().sqrt().sum() / 2.0;

	if (print)
	{
		std::cout << "MSE: " << error << std::endl;
		std::cout << "---------------------------------------------------------------------------------" << std::endl;
	}
	return error;
}

void NeuroNet::NeuralNetwork::CalcCorrectWeights(Problem& test)
{
	switch (_learn)
	{
	case BACKPROP:
		BackPropagation(test);
		break;
	case RPROP:
		ResilientPropagation(test);
		break;
	default:
		break;
	}
}

void NeuroNet::NeuralNetwork::ResilientPropagation(Problem & test)
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	int countLayers = (int)_layers.size();
	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < _layers[i].Grad.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < _layers[i].Grad.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double cur_mult = _layers[i].Grad(j,k) * _layers[i].LastGrad(j,k);
				if (cur_mult == 0.0)
					cur_correct = Common::getRand(0, 1);
				else if (cur_mult > 0.0)
					cur_correct = std::min(EttaPlus * _layers[i].CorrectVal(j,k), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = std::max(EttaMinus * _layers[i].CorrectVal(j,k), 1e-6);

				_layers[i].CorrectVal(j,k) = cur_correct;

				if (_layers[i].Grad(j,k) == 0.0) continue;

				if (_layers[i].Grad(j,k) > 0)
					_layers[i].Weights(j,k) += -cur_correct;
				else
					_layers[i].Weights(j,k) += cur_correct;
			}
		}

		for (int j = 0; j < _layers[i].Delta.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double cur_mult = _layers[i].Delta(0,j) * _layers[i].LastDelta(0,j);
			if (cur_mult == 0.0)
				cur_correct = Common::getRand(0, 1);
			else if (cur_mult > 0.0)
				cur_correct = std::min(EttaPlus * _layers[i].BiasCorrectVal(0,j), 50.0);
			else if (cur_mult < 0.0)
				cur_correct = std::max(EttaMinus * _layers[i].BiasCorrectVal(0,j), 1e-6);

			_layers[i].BiasCorrectVal(0,j) = cur_correct;

			if (_layers[i].Delta(0,j) == 0.0) continue;

			if (_layers[i].Delta(0,j) > 0)
				_layers[i].Bias(0,j) += -cur_correct;
			else
				_layers[i].Bias(0,j) += cur_correct;
		}
	}
}

void NeuroNet::NeuralNetwork::ResilientPropagation(std::vector<NeuroNet::Matrix2d> PrevSumGrad, std::vector<NeuroNet::Matrix2d> SumGrad, std::vector<NeuroNet::Matrix2d> PrevSumDelta, std::vector<NeuroNet::Matrix2d> SumDelta)
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	int countLayers = (int)_layers.size();
	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < SumGrad[i].GetVerticalSize(); ++j)
		{
			for (int k = 0; k < SumGrad[i].GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double cur_mult = SumGrad[i](j,k) * PrevSumGrad[i](j, k);
				if (cur_mult == 0.0)
					cur_correct = Common::getRand(0, 1);
				else if (cur_mult > 0.0)
					cur_correct = std::min(EttaPlus * _layers[i].CorrectVal(j,k), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = std::max(EttaMinus * _layers[i].CorrectVal(j,k), 1e-6);

				_layers[i].CorrectVal(j,k) = cur_correct;

				if (SumGrad[i](j, k) == 0.0) continue;

				if (SumGrad[i](j, k) > 0)
					_layers[i].Weights(j,k) += -cur_correct;
				else
					_layers[i].Weights(j,k) += cur_correct;
			}
		}

		for (int j = 0; j < SumDelta[i].GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double cur_mult = SumDelta[i](0,j) * PrevSumDelta[i](0,j);
			if (cur_mult == 0.0)
				cur_correct = Common::getRand(0, 1);
			else if (cur_mult > 0.0)
				cur_correct = std::min(EttaPlus * _layers[i].BiasCorrectVal(0,j), 50.0);
			else if (cur_mult < 0.0)
				cur_correct = std::max(EttaMinus * _layers[i].BiasCorrectVal(0,j), 1e-6);

			_layers[i].BiasCorrectVal(0,j) = cur_correct;

			if (SumDelta[i](0,j) == 0.0) continue;

			if (SumDelta[i](0,j) > 0)
				_layers[i].Bias(0,j) += -cur_correct;
			else
				_layers[i].Bias(0,j) += cur_correct;
		}
	}
}

void NeuroNet::NeuralNetwork::BackPropagation(Problem & test)
{
	int countLayers = (int)_layers.size();

	for (int i = 1; i < countLayers; ++i)
	{
		_layers[i].Weights += -_layers[i].Grad * EducationalSpeed;
		_layers[i].Bias += -_layers[i].Delta * EducationalSpeed;
	}
}

void NeuroNet::NeuralNetwork::Run(NeuroNet::Problem test)
{
	int countLayers = (int)_layers.size();
	//init input layer
	for (int i = 0; i < test.inputs.GetHorizontalSize(); ++i)
		_layers[0].States(0,i) = test.inputs(0,i);
	_layers[0].CalculateAxons();

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}

	for (int i = 0; i < countLayers; ++i)
	{
		_layers[i].LastDelta = _layers[i].Delta;
		_layers[i].LastGrad = _layers[i].Grad;
	}

	_layers.back().Delta = (_layers.back().Axons - test.outputs).multiplication(_layers.back().GetDiff());

	for (int i = countLayers - 2; i >= 0; --i)
		_layers[i].Delta = (_layers[i + 1].Delta * _layers[i + 1].Weights).multiplication(_layers[i].GetDiff());

	for (int i = 1; i < countLayers; ++i)
		_layers[i].Grad = !_layers[i].Delta * _layers[i - 1].Axons;
}

NeuroNet::Matrix2d NeuroNet::NeuralNetwork::GetOut() const
{
	return _layers.back().Axons;
}

std::ostream & NeuroNet::operator<<(std::ostream & os, NeuralNetwork & net)
{
	for (int i = 1; i < net._layers.size(); ++i)
	{
		os << "Weights " << i << " -> " << i - 1 << ": " << std::endl;
		os << net._layers[i].Weights << std::endl;
	}

	for (int i = 1; i < net._layers.size(); ++i)
	{
		os << "Bias " << i << " -> " << i - 1 << ": " << std::endl;
		os << net._layers[i].Bias << std::endl;
	}

	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Axons.GetHorizontalSize(); ++i)
		os << net._layers[0].Axons(0,i) << " ";

	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Axons.GetHorizontalSize(); ++i)
		os << net._layers.back().Axons(0,i) << " ";
	os << std::endl;
	os << "===============" << std::endl;
	return os;
}