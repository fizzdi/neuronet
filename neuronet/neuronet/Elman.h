#pragma once
#include "NeuralNetwork.h"

namespace NeuroNet
{
	class Elman : public NeuralNetwork
	{
	public:
		void Init(int InputCount, int OutputCount, int HiddenNeuronCount, AFType HiddenLayerFunction);
		double RunTrainingSet(bool print = false);
		std::vector<double> Run(Matrix2d &inputs);
	};



}