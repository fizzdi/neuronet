#pragma once
#include "NeuralNetwork.h"

namespace NeuroNet
{
	class Elman : public NeuralNetwork
	{
	public:
		void Init(int InputCount, int OutputCount, int HiddenNeuronCount, AFType HiddenLayerFunction, LearningType Learn);
		double RunTrainingSet(bool print = false);
		Matrix2d Run(Matrix2d &inputs);
		friend std::ostream& operator<< (std::ostream &os, Elman &net);
		void PrintProblemResult(Problem& test);
	};



}