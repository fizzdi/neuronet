#include "Elman.h"
#include <algorithm>

void NeuroNet::Elman::Init(int InputCount, int OutputCount, int HiddenNeuronCount, AFType HiddenLayerFunction)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount + HiddenNeuronCount, 0, LINE));
	_layers.push_back(Layer(HiddenNeuronCount, InputCount + HiddenNeuronCount, HiddenLayerFunction));
	_layers.push_back(Layer(OutputCount, HiddenNeuronCount, LINE));
}

double NeuroNet::Elman::RunTrainingSet(bool print)
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

		//copy last hiddent layer in input layer
		std::copy(
			_layers[_layers.size() - 2].Axons.rowbegin(),
			_layers[_layers.size() - 2].Axons.rowend(),
			_layers[0].States.rowbegin() + (_layers[0].States.GetVerticalSize() - _layers[_layers.size() - 2].States.GetVerticalSize())
		);

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
