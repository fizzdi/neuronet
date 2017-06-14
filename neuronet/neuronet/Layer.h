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
		Matrix2d Axons;
		Matrix2d States;
		Matrix2d Delta;
		Matrix2d LastDelta;
		Matrix2d LastGrad;
		Matrix2d CorrectVal;
		Matrix2d GetDiff();
		Matrix2d Weights;

		//----------------------------------
		//Functions
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		double GetDiff(double val);		
	};
}