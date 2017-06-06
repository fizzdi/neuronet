#include "Layer.h"

NeuroNet::Layer::Layer(int neuronCount, int prevNeuronCount)
{
	Weights.InitRandom(prevNeuronCount, neuronCount);
	Correct.Init(prevNeuronCount, neuronCount);
	_layer.resize(neuronCount);
}

NeuroNet::Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction)
{
	Weights.InitRandom(prevNeuronCount, neuronCount);
	Correct.Init(prevNeuronCount, neuronCount);
	_layer.resize(neuronCount, Neuron(activationFunction));
}

void NeuroNet::Layer::CalculateStates(Layer & prevLayer)
{
	for (int i = 0; i < this->Count(); ++i)
	{
		double newState = 0;
		for (int j = 0; j < prevLayer.Count(); ++j)
			newState += prevLayer[j].GetAxon() * this->Weights[j][i];
		this->_layer[i].SetState(newState);
	}
}

void NeuroNet::Layer::CalculateAxons()
{
	for (int i = 0; i < this->Count(); ++i)
		_layer[i].CalculateAxon();
}

const int NeuroNet::Layer::Count()
{
	return _layer.size();
}

NeuroNet::Neuron& NeuroNet::Layer::operator[](int n)
{
	return _layer[n];
}
