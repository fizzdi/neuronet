#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>

//Neural config
const int SENSOR_COUNT = 8;
const int INPUT_NEURON_COUNT = SENSOR_COUNT * 2;
const int HIDDEN_NEURON_COUNT = INPUT_NEURON_COUNT * 3;
const int OUTPUT_NEURON_COUNT = 1;
const int TEST_COUNT = 50;
const int TRAIN_EPOCH = 30;
const int TRAIN_PERIOD = 5;
const double TRAIN_EPS = 1e-3;
const int RANDOM_ACTION_PERIOD = 5;

//RMS config
const double RMS_GAMMA = 0.95;
const double RMS_LEARNRATE = 1e-3;
const double RMS_EPSILON = 1e-2;

namespace NeuroNet
{
	class NeuralNetwork
	{
	protected:
		int _countlayers;
		virtual void Run();
	public:
		std::vector<Layer> Layers;
		std::vector<Matrix2d> eligibility;

		NeuralNetwork() {};
		NeuralNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		void Run(const std::vector<double>& input);
		void Run(Matrix2d &input);
		void AddTest(std::deque<Problem> &TrainingSet, const std::vector<double> &input, const std::vector<double> &ideal) const;
		void AddTest(std::deque<Problem> &TrainingSet, const Matrix2d& input, const Matrix2d& ideal) const;
		double RunTrainingSetOffline(std::deque<Problem> &TrainingSet, bool print = false);
		double CalculateError(Problem& test, bool print = false);
		void PrintProblemResult(Problem& test);
		void ResilientPropagation();
		void ResilientPropagationOffline();
		double RMSTraining(std::deque<Problem> &TrainingSet);
		void RMSPropagation();
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
