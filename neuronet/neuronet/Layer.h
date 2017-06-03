#pragma once
#include "Neuron.h"
#include "Matrix.h"
#include <vector>
namespace NeuroNet
{
	class Layer
	{
	private:
		std::vector<NeuroNet::Neuron> _layer;
	public:
		Layer(int neuronCount, int prevNeuronCount);
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		const int Count();
		NeuroNet::Neuron& operator[] (int n);

		Matrix Weights;
		Matrix Correct;
		Matrix OldCorrect;
	};
}