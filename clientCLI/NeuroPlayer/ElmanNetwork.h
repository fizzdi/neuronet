#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>
#include "NeuralNetwork.h"

namespace NeuroNet
{
	class ElmanNetwork : public NeuralNetwork
	{
	protected:
		int context_neuron;
		void Run();
	public:
		using NeuralNetwork::Run;
		ElmanNetwork() {};
		ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		Matrix2d GetContext();
		void SetContext(const Matrix2d& context);
		double RMSTraining(training_set &TrainingSet);
	};		
}
