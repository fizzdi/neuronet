#include "NeuralNetwork.h"
#include <iostream>
#include <algorithm>

void NeuroNet::NeuralNetwork::Init(int InputCount, int OutputCount, int NeuronCount)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount, TANH));
	_layers.push_back(Layer(OutputCount, NeuronCount, TANH));
}

double NeuroNet::NeuralNetwork::RunProblemSet()
{
	double maxError = 0.0;
	for (int i = 0; i < _layers.size(); ++i)
	{
		_layers[i].OldCorrect = _layers[i].Correct;
		_layers[i].Correct.Clear();
	}

	for each (Problem test in ProblemSet)
	{
		//init input layer
		for (int i = 0; i < test.inputs.size(); ++i)
			_layers[0][i].SetState(test.inputs[i]);
		_layers[0].CalculateAxons();

		//init hidden and output layers
		for (int i = 1; i < _layers.size(); ++i)
		{
			_layers[i].CalculateStates(_layers[i - 1]);
			_layers[i].CalculateAxons();
		}

		PrintProblemResult(test);
		maxError = std::max(maxError, CalculateError(test));
		CalcCorrectWeights(test);
	}

	/*for (int i = 0; i < _layers.size(); ++i)
	{
		if (_layers[i].Correct.GetVerticalSize() != _layers[i].Weights.GetVerticalSize() ||
			_layers[i].Correct.GetHorizontalSize() != _layers[i].Weights.GetHorizontalSize())
			int y = 0;
		for (int j = 0; j < _layers[i].Correct.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < _layers[i].Correct.GetHorizontalSize(); ++k)
			{
				std::cout << _layers[i].Weights.get(j,k) << " " << _layers[i].Correct.get(j, k) << std::endl;
			}
		}
	}
	std::cout << std::endl;*/

	return maxError;
}

void NeuroNet::NeuralNetwork::PrintProblemResult(const Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.size(); ++i)
		std::cout << "\toutput[" << i << "] = " << test.outputs[i] << std::endl;
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
	std::cout << "---------------------------------------------------------------------------------" << std::endl;
	return error;
}

void NeuroNet::NeuralNetwork::CorrectWeights()
{
	for (int i = 0; i < _layers.size(); ++i)
	{
		for (int j = 0; j < _layers[i].Correct.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < _layers[i].Correct.GetHorizontalSize(); ++k)
			{
				_layers[i].Weights.add(j, k, _layers[i].Correct.get(j, k));
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
		_layers.back()[i].Delta((test.outputs[i] - _layers.back()[i].GetAxon()) * _layers.back()[i].GetDiff());
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
			double newDelta = _layers[i][j].GetDiff() * sum;
			_layers[i][j].Delta(newDelta);
		}
	}

	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < _layers[i].Count(); ++j)
		{
			for (int k = 0; k < _layers[i - 1].Count(); ++k)
			{
				double grad = _layers[i][j].Delta() * _layers[i - 1][k].GetAxon();
				//double grad = _layers[i-1][k].Delta() * _layers[i][j].GetAxon();
				//TODO Moment
				_layers[i].Correct.add(k, j, EducationalSpeed * grad + Alpha*_layers[i].OldCorrect.get(k,j));
			}
		}
	}
}

std::vector<double> NeuroNet::NeuralNetwork::Run(std::vector<double>& inputs)
{
	//init input layer
	for (int i = 0; i < inputs.size(); ++i)
		_layers[0][i].SetState(inputs[i]);
	_layers[0].CalculateAxons();

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}
	std::vector<double> ans;
	for (int i = 0; i < _layers.back().Count(); ++i)
	{
		ans.push_back(_layers.back()[i].GetAxon());
	}
	return ans;
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
