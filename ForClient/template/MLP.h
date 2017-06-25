#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>
#include <string>
#include "NeuralNetwork.h"

//Neural config
static int SENSOR_COUNT = {0};
static int INPUT_NEURON_COUNT = SENSOR_COUNT * 2;
static int HIDDEN_NEURON_COUNT = {1};
static int HIDDEN_LAYER_COUNT = {2};
static int OUTPUT_NEURON_COUNT = SENSOR_COUNT;
static int TEST_COUNT = {3};
static int MAX_TEST_COUNT = {4};
static int TRAIN_EPOCH = {5};
static int TRAIN_PERIOD = {6};
static double TRAIN_EPS = {7};
static int END_TRAIN_TICK = {8};
static NeuroNet::AFType FUN_ACT = (NeuroNet::AFType){9};

//RMS config
static double RMS_GAMMA = {10};
static double RMS_LEARNRATE = {11};
static double RMS_EPSILON = {12};

namespace NeuroNet
{
	class MLP : public NeuralNetwork
	{
	protected:
		//Methods
		virtual void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		virtual void Run();
		//Training
		void ResilientPropagation();
		void RMSPropagation();
		//Calc
		virtual double CalculateError(Problem& test, bool print = false);
		virtual void CalcGradDelta(const double output);
		virtual void CalcGradDelta(const std::vector<double>& output);
		virtual void CalcGradDelta(const Matrix2d &output); 
	
	public:
		//Constructors
		MLP() {};
		MLP(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		MLP(const std::string& FileName);

		//Methods
		virtual void Run(const std::vector<double>& input);
		virtual void Run(const Matrix2d &input);

		//Training
		double RPROPTraining(training_set &TrainingSet);
		double RMSTraining(training_set &TrainingSet);

		//Gets
		virtual void PrintFullInfo(std::ostream &debug) const;
		virtual Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, MLP &net);

		//Sets
		virtual void SetInput(const Matrix2d& inputs);
	};
}
