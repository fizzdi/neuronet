#include "NeuralNetwork.h"
#include <iostream>
#include <algorithm>

void NeuroNet::NeuralNetwork::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount, HiddenLayerFunction));
	_layers.push_back(Layer(OutputCount, NeuronCount, LINE));
}

//TODO global var debug
double NeuroNet::NeuralNetwork::RunTrainingSet(bool print)
{
	double maxError = 0.0;
	for (int i = 0; i < _layers.size(); ++i)
		_layers[i].Correct.Clear();

	for each (Problem test in TrainingSet)
	{
		//init input layer
		std::copy(
			test.inputs.rowbegin(),
			test.inputs.rowend(),
			_layers[0].States.rowbegin()
			);
		_layers[0].CalculateAxons();

		//init hidden and output layers
		for (int i = 1; i < _layers.size(); ++i)
		{
			_layers[i].CalculateStates(_layers[i - 1]);
			_layers[i].CalculateAxons();
		}
		if (print)
			PrintProblemResult(test);
		maxError = std::max(maxError, CalculateError(test, print));
		CalcCorrectWeights(test);
	}

	return maxError;
}

void NeuroNet::NeuralNetwork::PrintProblemResult(Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.GetHorizontalSize(); ++i)
		std::cout << test.outputs[0][i] << " ";
	std::cout << std::endl;
}

double NeuroNet::NeuralNetwork::CalculateError(Problem & test, bool print)
{
	//MSE
	double error = (test.outputs - _layers.back().Axons).multiplication(test.outputs - _layers.back().Axons).sum() / 2;
	if (print)
	{
		std::cout << "Error: " << error << std::endl;
		std::cout << "---------------------------------------------------------------------------------" << std::endl;
	}
	return error;
}

void NeuroNet::NeuralNetwork::CorrectWeights()
{
	for (int i = 0; i < _layers.size(); ++i)
	{
		_layers[i].Weights += _layers[i].Correct;
		_layers[i].CorrectVal = _layers[i].Correct;
		_layers[i].Correct.Clear();

	}
}

void NeuroNet::NeuralNetwork::CalcCorrectWeights(Problem& test)
{
	//RPROP
	int countLayers = _layers.size();

	//calc delta
	_layers.back().Delta = (test.outputs - _layers.back().Axons).multiplication(_layers.back().GetDiff());

	for (int i = countLayers - 2; i >= 0; --i)
		_layers[i].Delta = (_layers[i + 1].Delta * _layers[i + 1].Weights).multiplication(_layers[i].GetDiff());

	for (int i = 1; i < countLayers; ++i)
	{
		//_layers[i].Grad = !_layers[i].Delta * _layers[i - 1].Axons;
		//TODO STL
		Matrix2d Grad = !_layers[i].Delta * _layers[i - 1].Axons;
		for (int j = 0; j < Grad.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < Grad.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double cur_mult = Grad[j][k] * _layers[i].LastGrad[j][k];
				if (cur_mult >= 0.0)
					cur_correct = 1.2 * _layers[i].CorrectVal[j][k];
				else if (cur_mult < 0.0)
					cur_correct = 0.5 * _layers[i].CorrectVal[j][k];

				cur_correct = std::min(cur_correct, 50.0);
				cur_correct = std::max(cur_correct, 1e-6);
				_layers[i].CorrectVal[j][k] = cur_correct;

				if (Grad[j][k] == 0)
					_layers[i].Correct[j][k] = 0.0;
				else if (Grad[j][k] > 0)
					_layers[i].Weights[j][k] += -cur_correct;
				else if (Grad[j][k] < 0)
					_layers[i].Weights[j][k] += cur_correct;
				int y = 0;
			}
		}
		_layers[i].LastGrad = Grad;
		//_layers[i].Weights += _layers[i].Correct;
		//_layers[i].Correct.Clear();
	}
}

NeuroNet::Matrix2d NeuroNet::NeuralNetwork::Run(NeuroNet::Matrix2d& inputs)
{
	//init input layer
	std::copy(
		inputs.rowbegin(),
		inputs.rowend(),
		_layers[0].States.rowbegin()
		);
	_layers[0].CalculateAxons();

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}
	return _layers.back().Axons;
}

std::ostream & NeuroNet::operator<<(std::ostream & os, NeuralNetwork & net)
{
	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Axons.GetHorizontalSize(); ++i)
		os << net._layers[0].Axons[0][i] << " ";
	os << std::endl;

	os << std::endl << "Hidden axons:" << std::endl;
	for (int i = 0; i < net._layers[1].Axons.GetHorizontalSize(); ++i)
		os << net._layers[1].Axons[0][i] << " ";
	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Axons.GetHorizontalSize(); ++i)
		os << net._layers.back().Axons[0][i] << " ";
	os << std::endl;

	os << std::endl;
	return os;
}