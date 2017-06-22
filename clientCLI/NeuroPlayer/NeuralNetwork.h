#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>

//Neural config
const int SENSOR_COUNT = 8;
const int INPUT_NEURON_COUNT = SENSOR_COUNT * 2;
const int HIDDEN_NEURON_COUNT = INPUT_NEURON_COUNT + 10;
const int OUTPUT_NEURON_COUNT = 1;
const int TEST_COUNT = 500;
const int TRAIN_EPOCH = 30;
const int TRAIN_PERIOD = 10;
const double TRAIN_EPS = 1e-3;

namespace NeuroNet
{
	class NeuralNetwork
	{
	protected:
		int _countlayers;
		std::vector<Layer> _layers;
		virtual void Run();
	public:
		std::vector<Matrix2d> eligibility;

		NeuralNetwork() {};
		NeuralNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
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
		void debuginfo(std::ostream &debug) const;
		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);
	};
}
