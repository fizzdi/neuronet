#pragma once
#include "NeuralNetwork.h"

namespace NeuroNet
{
	class Elman : public NeuralNetwork
	{
	public:
		void Init(int InputCount, int OutputCount, int HiddenNeuronCount);
		double RunTrainingSet(bool print = false);

	};



}