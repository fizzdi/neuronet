#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>


//Neural config
const int SENSOR_COUNT = 8;
const int INPUT_NEURON_COUNT = SENSOR_COUNT * 2;
const int HIDDEN_NEURON_COUNT = 40;
const int OUTPUT_NEURON_COUNT = 1;
const int TEST_COUNT = 300;
const int TRAIN_EPOCH = 10;
const int TRAIN_PERIOD = 5;
const double TRAIN_EPS = 1e-1;

namespace NeuroNet
{
	class ElmanNetwork
	{
	protected:
		int _countlayers;
	public:
		std::vector<Layer> _layers;
		ElmanNetwork() {};
		ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		std::vector<Matrix2d> eligibility;
		void Run();
		void Run(const std::vector<double>& input);
		void Run(Matrix2d &input);
		void AddTest(std::deque<Problem> &TrainingSet, const std::vector<double>& ideal) const;
		void AddTest(std::deque<Problem> &TrainingSet, Matrix2d& ideal) const;
		double RunTrainingSetOffline(std::deque<Problem> &TrainingSet, bool print = false);
		double CalculateError(Problem& test, bool print = false);
		void PrintProblemResult(Problem& test);
		void ResilientPropagation();
		void ResilientPropagationOffline();
		void CalcGradDelta(const double output);
		void CalcGradDelta(const std::vector<double>& output);
		void CalcGradDelta(Matrix2d &output);
		Matrix2d& GetOutputGrad();
		void MCQLCorrect();
		void debuginfo(std::ostream &debug) const
		{
			for (int i = 1; i < _countlayers; ++i)
				debug << _layers[i].Weights << std::endl;
			debug << std::endl;
		}
		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, ElmanNetwork &net);
	};
}
