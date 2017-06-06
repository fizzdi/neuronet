#pragma once
#include "Neuron.h"
#include "Matrix.h"
#include <vector>
namespace NeuroNet
{
	//activation function type
	enum AFType { LINE, SIGM, TANH };

	class Layer
	{
	private:
		AFType _aftype;

		Matrix sigm_function(Matrix m);
		Matrix tanh_function(Matrix m);
		Matrix diff_tanh_function(Matrix m);
		Matrix diff_sigm_function(Matrix m);
	public:
		NeuroNet::Matrix Axons;
		NeuroNet::Matrix States;
		NeuroNet::Matrix Delta;
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		Matrix GetDiff();
		double GetDiff(double val);
		Matrix Weights;
		Matrix Correct;
		Matrix OldCorrect;
	};
}