#pragma once
#include "Matrix2d.h"

namespace NeuroNet
{
	enum AFType { LINE, SIGM, TANH };

	class Layer
	{
	private:
		AFType _aftype;

		Matrix2d sigm_function(Matrix2d& m);
		Matrix2d tanh_function(Matrix2d& m);
		Matrix2d diff_tanh_function(Matrix2d& m);
		Matrix2d diff_sigm_function(Matrix2d& m);
	public:
		Matrix2d Axons;
		Matrix2d States;
		Matrix2d Delta;
		Matrix2d LastDelta;
		Matrix2d Grad;
		Matrix2d GradSum;
		Matrix2d DeltaSum;
		Matrix2d LastGrad;
		Matrix2d CorrectVal;
		Matrix2d BiasCorrectVal;
		Matrix2d Weights;
		Matrix2d Bias;
		Matrix2d LastDeltaSum;
		Matrix2d LastGradSum;

		//----------------------------------
		//Functions
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		void NguenWidrow(double Xmin, double Xmax, double Ymin, double Ymax);
		Matrix2d GetDiff();
	};
}