#include "ElmanNetwork.h"
#include <algorithm>

using namespace NeuroNet;

ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	_layers.clear();
	_layers.emplace_back(Layer(InputCount + NeuronCount, 0, LINE));
	_layers.emplace_back(Layer(NeuronCount, InputCount + NeuronCount, HiddenLayerFunction));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	_layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	_layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	_layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	_countlayers = _layers.size();

	eligibility.resize(_countlayers);
	for (int i = 0; i < _countlayers; ++i)
	{
		eligibility[i] = Matrix2d(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
		eligibility[i].Fill(0.0);
	}
}

void ElmanNetwork::Run()
{
	//init hidden and output layers
	for (int i = 1; i < (int)_layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}

	//copy hidden into input (context)
	for (int i = 1; i <= _layers[1].Axons.GetHorizontalSize(); ++i)
		_layers[0].States.at(0, _layers[0].States.GetHorizontalSize() - i) = _layers[1].Axons.at(0, _layers[1].Axons.GetHorizontalSize() - i);
}