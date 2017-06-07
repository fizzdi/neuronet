#pragma once
#include "Neuron.h"
#include "Matrix2d.h"
#include <vector>
namespace NeuroNet
{
	//activation function type
	enum AFType { LINE, SIGM, TANH };

	class Layer
	{
	private:
		AFType _aftype;

		Matrix2d sigm_function(Matrix2d m);
		Matrix2d tanh_function(Matrix2d m);
		Matrix2d diff_tanh_function(Matrix2d m);
		Matrix2d diff_sigm_function(Matrix2d m);
	public:
		NeuroNet::Matrix2d Axons;
		NeuroNet::Matrix2d States;
		NeuroNet::Matrix2d Delta;
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		Matrix2d GetDiff();
		double GetDiff(double val);
		Matrix2d Weights;
		Matrix2d Correct;
	};
}