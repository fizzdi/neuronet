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
	for (int i = 0; i < test.outputs.GetVerticalSize(); ++i)
		std::cout << "\toutput[" << i << "] = " << test.outputs[i][0] << std::endl;
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
	}
}

void NeuroNet::NeuralNetwork::CalcCorrectWeights(Problem& test)
{
	//PROP
	int countLayers = _layers.size();

	//calc delta
	_layers.back().Delta = (test.outputs - _layers.back().Axons).multiplication(_layers.back().GetDiff());

	for (int i = countLayers - 2; i >= 0; --i)
		_layers[i].Delta = (_layers[i + 1].Delta * _layers[i + 1].Weights).multiplication(_layers[i].GetDiff());

	for (int i = 1; i < countLayers; ++i)
		_layers[i].Correct += !_layers[i].Delta * _layers[i - 1].Axons * EducationalSpeed;
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
	for (int i = 0; i < net._layers[0].Axons.GetVerticalSize(); ++i)
		os << "\tinput[" << i << "] = " << net._layers[0].Axons[i][0] << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Axons.GetVerticalSize(); ++i)
		os << "\toutput[" << i << "] = " << net._layers.back().Axons[i][0] << std::endl;
	os << std::endl;
	return os;
}