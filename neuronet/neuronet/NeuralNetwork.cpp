#include "NeuralNetwork.h"
#include <iostream>
#include <sstream>
#include <algorithm>

void NeuroNet::NeuralNetwork::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount, HiddenLayerFunction, true));
	_layers.back().NguenWidrow();
	_layers.push_back(Layer(OutputCount, NeuronCount, LINE, true));
	_layers.back().NguenWidrow();

	_learn = Learn;
}

//TODO global var debug
double NeuroNet::NeuralNetwork::RunTrainingSet(bool print)
{
	double error = 0.0;
	for each (Problem test in TrainingSet)
	{
		Run(test.inputs);
		error = CalculateError(test);
		CalcCorrectWeights(test);		
	}
	return error;
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
	double error = (test.outputs - _layers.back().Axons).multiplication(test.outputs - _layers.back().Axons).sum()/2.0;

	//RootMSE
	//double error = (test.outputs - _layers.back().Axons).abs().sqrt().sum() / 2.0;

	if (print)
	{
		std::cout << "Error: " << error << std::endl;
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
	int countLayers = _layers.size();

	//calc delta
	_layers.back().Delta = (test.outputs - _layers.back().Axons).multiplication(_layers.back().GetDiff());

	for (int i = countLayers - 2; i >= 0; --i)
		_layers[i].Delta = (_layers[i + 1].Delta * _layers[i + 1].Weights).multiplication(_layers[i].GetDiff());

	for (int i = 1; i < countLayers; ++i)
	{
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
				cur_correct = std::max(cur_correct, 1e-7);
				_layers[i].CorrectVal[j][k] = cur_correct;

				if (Grad[j][k] == 0.0) continue;

				if (Grad[j][k] > 0)
				{
					_layers[i].Weights[j][k] += -cur_correct;
					_layers[i].Bias[0][j] += -cur_correct;
				}
				else
				{
					_layers[i].Weights[j][k] += cur_correct;
					_layers[i].Bias[0][j] += cur_correct;
				}
				int y = 0;
			}

		}
		_layers[i].LastGrad = Grad;
	}
}

void NeuroNet::NeuralNetwork::BackPropagation(Problem & test)
{
	int countLayers = _layers.size();

	//calc delta
	_layers.back().Delta = (test.outputs - _layers.back().Axons).multiplication(_layers.back().GetDiff());

	for (int i = countLayers - 2; i >= 0; --i)
		_layers[i].Delta = (_layers[i + 1].Delta * _layers[i + 1].Weights).multiplication(_layers[i].GetDiff());

	for (int i = 1; i < countLayers; ++i)
	{
		_layers[i].Weights += !_layers[i].Delta * _layers[i - 1].Axons * EducationalSpeed;
		_layers[i].Bias += _layers[i].Delta * EducationalSpeed;
	}
}

void NeuroNet::NeuralNetwork::Run(NeuroNet::Matrix2d& inputs)
{
	//init input layer
	for (int i = 0; i < inputs.GetHorizontalSize(); ++i)
		_layers[0].States[0][i] = inputs[0][i];
	_layers[0].CalculateAxons();

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}

	//DEBUG----------------------------------
	for (int i = 0; i < _layers.size(); ++i)
	{
		for (int j = 0; j < _layers[i].Axons.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < _layers[i].Axons.GetHorizontalSize(); ++k)
			{
				std::ostringstream ss;
				ss << _layers[i].Axons[j][k];
				auto a = ss.str();
				ss.str("");
				//ss.clear();
				ss << -nan("");
				auto b = ss.str();
				if (a == b)
				{
					int y = 0;
				}
			}
		}
	}
}

NeuroNet::Matrix2d NeuroNet::NeuralNetwork::GetOut() const
{
	return _layers.back().Axons;
}

std::ostream & NeuroNet::operator<<(std::ostream & os, NeuralNetwork & net)
{
	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Axons.GetHorizontalSize(); ++i)
		os << net._layers[0].Axons[0][i] << " ";

	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Axons.GetHorizontalSize(); ++i)
		os << net._layers.back().Axons[0][i] << " ";
	os << std::endl;
	os << std::endl;
	return os;
}