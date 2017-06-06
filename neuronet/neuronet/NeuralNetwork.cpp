#include "NeuralNetwork.h"
#include <iostream>
#include <algorithm>

void NeuroNet::NeuralNetwork::Init(int InputCount, int OutputCount, int NeuronCount)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount, TANH));
	_layers.push_back(Layer(OutputCount, NeuronCount, LINE));
}

double NeuroNet::NeuralNetwork::RunTrainingSet(bool print)
{
	double maxError = 0.0;
	for (int i = 0; i < _layers.size(); ++i)
	{
		_layers[i].OldCorrect = _layers[i].Correct;
		_layers[i].Correct.Clear();
	}

	for each (Problem test in TrainingSet)
	{
		//init input layer
		_layers[0].States = test.inputs;
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
	for (int i = 0; i < test.outputs.size(); ++i)
		std::cout << "\toutput[" << i << "] = " << test.outputs[i][0] << std::endl;
}

double NeuroNet::NeuralNetwork::CalculateError(Problem & test, bool print)
{
	//MSE
	double error = ((test.outputs - _layers.back().Axons)*(test.outputs - _layers.back().Axons)).sum() / 2;
	if (print)
	{
		std::cout << "Error: " << error * 100 << "% (" << error << ")" << std::endl;
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
	int countLayers = _layers.size();
	for (int i = 0; i < _layers.back().Axons.size(); ++i)
	{
		_layers.back().Delta[i][0] = ((test.outputs[i][0] - _layers.back().Axons[i][0]) * _layers.back().GetDiff(_layers.back().Axons[i][0]));
	}
	//int countLayers = _layers.size();
	//calc delta
	_layers.back().Delta = (test.outputs - _layers.back().Axons) * _layers.back().GetDiff();

	for (int i = countLayers - 2; i >= 0; --i)
	{
		for (int j = 0; j < _layers[i].Axons.size(); ++j)
		{
			//calc sum(Wi * deltai)
			double sum = 0.0;
			for (int k = 0; k < _layers[i + 1].Axons.size(); ++k)
			{
				sum += _layers[i + 1].Delta[k][0] * _layers[i + 1].Weights[j][k];
			}
			_layers[i].Delta[j][0] = _layers[i].GetDiff(_layers[i].Axons[j][0]) * sum;
		}
	}

	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < _layers[i].Axons.size(); ++j)
		{
			for (int k = 0; k < _layers[i - 1].Axons.size(); ++k)
			{
				double grad = _layers[i].Delta[j][0] * _layers[i - 1].Axons[k][0];
				//double grad = _layers[i-1][k].Delta() * _layers[i][j].GetAxon();
				//TODO Moment
				_layers[i].Correct[k][j] += EducationalSpeed * grad + Alpha*_layers[i].OldCorrect[k][j];
			}
		}
	}
}

std::vector<double> NeuroNet::NeuralNetwork::Run(NeuroNet::Matrix& inputs)
{
	//init input layer
	_layers[0].States = inputs;
	_layers[0].CalculateAxons();

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}
	std::vector<double> ans;
	for (int i = 0; i < _layers.back().Axons.size(); ++i)
	{
		ans.push_back(_layers.back().Axons[i][0]);
	}
	return ans;
}

std::ostream & NeuroNet::operator<<(std::ostream & os, NeuralNetwork & net)
{
	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Axons.size(); ++i)
		os << "\tinput[" << i << "] = " << net._layers[0].Axons[i][0] << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Axons.size(); ++i)
		os << "\toutput[" << i << "] = " << net._layers.back().Axons[i][0] << std::endl;
	os << std::endl;
	return os;
}