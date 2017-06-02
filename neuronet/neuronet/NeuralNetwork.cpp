#include "NeuralNetwork.h"
#include <iostream>
#include <algorithm>

void NeuroNet::NeuralNetwork::Init(int InputCount, int OutputCount, int NeuronCount)
{
	_layers.resize(0, Layer(0, 0));
	_layers.push_back(Layer(InputCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount));
	_layers.push_back(Layer(OutputCount, NeuronCount));
}

double NeuroNet::NeuralNetwork::RunProblemSet()
{
	double maxError = 0.0;
	for each (Problem test in ProblemSet)
	{
		//init input layer
		for (int i = 0; i < test.inputs.size(); ++i)
			_layers[0][i].SetState(test.inputs[i]);
		_layers[0].CalculateAxons();

		//init hidden and output layers
		for (int i = 1; i < _layers.size(); ++i)
		{
			correct[i].Clear();
			_layers[i].CalculateStates(_layers[i - 1]);
			_layers[i].CalculateAxons();
		}

		PrintProblemResult(test);
		maxError = std::max(maxError, CalculateError(test));
		CalcCorrectWeights(test);
	}
	return maxError;
}

void NeuroNet::NeuralNetwork::PrintProblemResult(const Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.size(); ++i)
		std::cout << "\toutput[" << i << "] = " << test.outputs[i] << std::endl;
	std::cout << "---------------------------" << std::endl;
}

double NeuroNet::NeuralNetwork::CalculateError(const Problem & test)
{
	//MSE
	double error = 0.0;
	for (int i = 0; i < test.outputs.size(); ++i)
	{
		error += (test.outputs[i] - _layers.back()[i].GetAxon())*(test.outputs[i] - _layers.back()[i].GetAxon());
	}
	error /= test.outputs.size();
	std::cout << "Error: " << error * 100 << "% (" << error << ")" << std::endl;
	return error;
}

void NeuroNet::NeuralNetwork::CorrectWeights()
{
	for (int i = 0; i < correct.size(); ++i)
	{
		for (int j = 0; j < correct[i].GetVerticalSize(); ++j)
		{
			for (int k = 0; k < correct[i].GetHorizontalSize(); ++k)
			{
				_layers[i].Weights.add(j, k, correct[i].get(j, k));
			}
		}
	}
}

void NeuroNet::NeuralNetwork::CalcCorrectWeights(const Problem& test)
{
	int countLayers = _layers.size();
	//calc delta
	for (int i = 0; i < _layers.back().Count(); ++i)
	{
		_layers.back()[i].Delta((test.outputs[i] - _layers.back()[i].GetAxon()) * _layers.back()[i].GetDiff(_layers.back()[i].GetAxon()));
	}

	for (int i = countLayers - 2; i >= 0; --i)
	{
		for (int j = 0; j < _layers[i].Count(); ++j)
		{
			//calc sum(Wi * deltai)
			double sum = 0.0;
			for (int k = 0; k < _layers[i + 1].Count(); ++k)
			{
				sum += _layers[i + 1][k].Delta() * _layers[i + 1].Weights.get(j, k);
			}
			double newDelta = _layers[i][j].GetDiff(_layers[i][j].GetAxon()) * sum;
			_layers[i][j].Delta(_layers[i][j].GetDiff(_layers[i][j].GetAxon()) * sum);
		}
	}

	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < _layers[i].Count(); ++j)
		{
			for (int k = 0; k < _layers[i - 1].Count(); ++k)
			{
				double grad = _layers[i][j].Delta() * _layers[i - 1][k].GetAxon();
				//TODO Moment
				correct[i].add(k, j, EducationalSpeed * grad);
			}
		}
	}
}

std::ostream & NeuroNet::operator<<(std::ostream & os, NeuralNetwork & net)
{
	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Count(); ++i)
		os << "\tinput[" << i << "] = " << net._layers[0][i].GetAxon() << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Count(); ++i)
		os << "\toutput[" << i << "] = " << net._layers.back()[i].GetAxon() << std::endl;
	os << std::endl;
	return os;
}
