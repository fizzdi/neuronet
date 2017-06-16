#pragma once
#include "NeuralNetwork.h"

namespace NeuroNet 
{
	class ElmanNetwork : public NeuralNetwork
	{
	public:
		void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn);
		void Run(Problem test);
	};
}
