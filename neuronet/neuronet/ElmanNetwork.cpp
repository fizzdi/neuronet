#include "ElmanNetwork.h"

void NeuroNet::ElmanNetwork::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn)
{
	_layers.clear();
	_layers.push_back(Layer(InputCount+ NeuronCount, 0, LINE));
	_layers.push_back(Layer(NeuronCount, InputCount + NeuronCount, HiddenLayerFunction, true));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	_layers.push_back(Layer(OutputCount, NeuronCount, LINE, true));
	//_layers.back().NguenWidrow(-1, 1, -1, 1);
	_countlayers = _layers.size();

	_learn = Learn;

	LastDeltaSum.resize(3);
	LastGradSum.resize(3);
	for (int i = 0; i < 3; ++i)
	{
		LastGradSum[i].Init(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
		LastDeltaSum[i].Init(_layers[i].Delta.GetVerticalSize(), _layers[i].Delta.GetHorizontalSize());
	}
}

void NeuroNet::ElmanNetwork::Run(Problem test)
{
	//init input layer
	for (int i = 0; i < test.inputs.GetHorizontalSize(); ++i)
		_layers[0].States(0, i) = test.inputs(0, i);
	_layers[0].CalculateAxons();

	//init hidden and output layers
	for (int i = 1; i < _layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}

	//copy hidden into input (context)
	for (int i = 1; i <= _layers[1].Axons.GetHorizontalSize(); ++i)
		_layers[0].States(0, _layers[0].States.GetHorizontalSize() - i) = _layers[1].Axons(0, _layers[1].Axons.GetHorizontalSize() - i);

	CalcGradDelta(test);
}
