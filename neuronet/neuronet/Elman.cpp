#include "Elman.h"
#include <algorithm>

void NeuroNet::Elman::Init(int InputCount, int OutputCount, int HiddenNeuronCount, AFType HiddenLayerFunction, LearningType Learn)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount, 0, LINE));
	for (int i = InputCount; i < _layers[0].States.GetHorizontalSize(); ++i)
		_layers[0].States[0][i] = 0.5;
	_layers.push_back(Layer(HiddenNeuronCount, InputCount, HiddenLayerFunction));
	_layers.push_back(Layer(OutputCount, HiddenNeuronCount, SIGM));
	_layers.push_back(Layer(HiddenNeuronCount, HiddenNeuronCount, SIGM));
	_learn = Learn;
}

double NeuroNet::Elman::RunTrainingSet(bool print)
{
	double maxError = 0.0;

	for each (Problem test in TrainingSet)
	{
		//init input layer
		_layers[0].States = test.inputs;
		_layers[0].CalculateAxons();

		_layers[1].States = _layers[0].Axons * !_layers[1].Weights + _layers[3].Axons * !_layers[3].Weights;
		_layers[1].CalculateAxons();

		_layers[2].CalculateStates(_layers[1]);
		_layers[2].CalculateAxons();

		//copy last hiddent layer in input layer
		_layers[3].Axons = _layers[1].Axons;

		if (print)
			PrintProblemResult(test);
		maxError = std::max(maxError, CalculateError(test, print));
		CalcCorrectWeights(test);
	}

	return maxError;
}

NeuroNet::Matrix2d NeuroNet::Elman::Run(Matrix2d & inputs)
{
	//init input layer
	for (int i = 0; i < inputs.GetHorizontalSize(); ++i)
		_layers[0].Axons[0][i] = inputs[0][i];
	_layers[0].CalculateAxons();

	_layers[1].States = _layers[0].Axons * !_layers[1].Weights;
	_layers[1].States +=_layers[3].Axons * !_layers[3].Weights;
	_layers[1].CalculateAxons();

	_layers[2].CalculateStates(_layers[1]);
	_layers[2].CalculateAxons();

	//copy last hiddent layer in input layer
	_layers[3].Axons = _layers[1].Axons;

	return _layers[2].Axons;
}


std::ostream & NeuroNet::operator<<(std::ostream & os, NeuroNet::Elman & net)
{
	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Axons.GetHorizontalSize(); ++i)
		os << net._layers[0].Axons[0][i] << " ";
	os << std::endl;

	os << std::endl << "Context axons:" << std::endl;
	for (int i = 0; i < net._layers[3].Axons.GetHorizontalSize(); ++i)
		os << net._layers[3].Axons[0][i] << " ";
	os << std::endl;

	os << std::endl << "Hidden axons:" << std::endl;
	for (int i = 0; i < net._layers[1].Axons.GetHorizontalSize(); ++i)
		os << net._layers[1].Axons[0][i] << " ";
	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers[2].Axons.GetHorizontalSize(); ++i)
		os << net._layers[2].Axons[0][i] << " ";
	os << std::endl;

	os << std::endl;
	return os;
}

void NeuroNet::Elman::PrintProblemResult(Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.GetHorizontalSize(); ++i)
		std::cout << test.outputs[0][i] << " ";
	std::cout << std::endl;
}