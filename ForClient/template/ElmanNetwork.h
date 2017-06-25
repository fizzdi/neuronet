#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>
#include "MLP.h"

namespace NeuroNet
{
	class ElmanNetwork : public MLP
	{
	protected:
		int context_neuron;
		virtual void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		virtual void Run();
	public:
		using MLP::Run;
		//construct 
		ElmanNetwork() {};
		ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		ElmanNetwork(const std::string& FileName);

		//trainings
		double RMSTraining(training_set &TrainingSet);

		//Gets/Sets
		Matrix2d GetContext();
		void SetContext(const Matrix2d& context);
	};		
}
