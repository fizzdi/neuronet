#include "Elman.h"
#include <algorithm>

void NeuroNet::Elman::Init(int InputCount, int OutputCount, int HiddenNeuronCount, AFType HiddenLayerFunction)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount + HiddenNeuronCount, 0, LINE));
	for (int i = InputCount; i < _layers[0].States.GetHorizontalSize(); ++i)
		_layers[0].States[0][i] = 0.5;
	_layers.push_back(Layer(HiddenNeuronCount, InputCount + HiddenNeuronCount, HiddenLayerFunction));
	_layers.push_back(Layer(OutputCount, HiddenNeuronCount, SIGM));
}

double NeuroNet::Elman::RunTrainingSet(bool print)
{
	double maxError = 0.0;
	for (int i = 0; i < _layers.size(); ++i)
		_layers[i].Correct.Clear();

	for each (Problem test in TrainingSet)
	{
		//init input layer
		for (int i = 0; i < test.inputs.GetHorizontalSize(); ++i)
			_layers[0].Axons[0][i] = test.inputs[0][i];
		_layers[0].CalculateAxons();

		//init hidden and output layers
		for (int i = 1; i < _layers.size(); ++i)
		{
			_layers[i].CalculateStates(_layers[i - 1]);
			_layers[i].CalculateAxons();
		}

		//copy last hiddent layer in input layer
		for (int i = 0; i < _layers[_layers.size() - 2].Axons.GetHorizontalSize(); ++i)
			_layers[0].States[0][_layers[0].States.GetHorizontalSize() - i - 1] = _layers[_layers.size() - 2].Axons[0][_layers[_layers.size() - 2].Axons.GetHorizontalSize() - i - 1];

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

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}

	//copy last hiddent layer in input layer
	for (int i = 0; i < _layers[_layers.size() - 2].Axons.GetHorizontalSize(); ++i)
		_layers[0].States[0][_layers[0].States.GetHorizontalSize() - i - 1] = _layers[_layers.size() - 2].Axons[0][_layers[_layers.size() - 2].Axons.GetHorizontalSize() - i - 1];

	return _layers.back().Axons;
}
