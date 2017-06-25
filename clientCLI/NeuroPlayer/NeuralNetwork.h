#pragma once
#include <vector>
#include "Layer.h"
#include <iostream>
#include <fstream>
#include <deque>
#include <string>

//World config
const int PARAMS_COUNT = 3; //HEALTH, FULLNESS, ANGLE
const int FOOD_COUNT = 100;
const int POISON_COUNT = 0;
const int TRAP_COUNT = 0;
const int CORNUCOPIA_COUNT = 0;
const int BLOCK_COUNT = 0;
const int PLAYER_COUNT = 1;

//Neural config
static int SENSOR_COUNT = 16;
static int INPUT_NEURON_COUNT = SENSOR_COUNT * 2;
static int HIDDEN_NEURON_COUNT = INPUT_NEURON_COUNT*4;
static int HIDDEN_LAYER_COUNT = 2;
static int OUTPUT_NEURON_COUNT = SENSOR_COUNT;
static int TEST_COUNT = 15;
static int MAX_TEST_COUNT = 4000;
static int TRAIN_EPOCH = 1;
static int TRAIN_PERIOD = 3;
static double TRAIN_EPS = 1e-3;
static int END_TRAIN_TICK = 10000;

//RMS config
static double RMS_GAMMA = 0.9;
static double RMS_LEARNRATE = 1e-3;
static double RMS_EPSILON = 1e-2;

namespace NeuroNet
{
	typedef std::deque<Problem> training_set;
	void AddTest(training_set &TrainingSet, const std::vector<double> &input, const std::vector<double> &ideal);
	void AddTest(training_set &TrainingSet, const Matrix2d& input, const Matrix2d& ideal);

	class NeuralNetwork
	{
	protected:
		int countLayers;
		std::vector<Layer> layers;
		int inputNeuron;
		int outputNeuron;
		
		//Methods
		virtual void Run();
	public:
		//Constructors
		NeuralNetwork() {};
		NeuralNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		NeuralNetwork(std::string FileName);
		
		//Methods
		void Run(const std::vector<double>& input);
		void Run(const Matrix2d &input);

		//Training
		double RunTrainingSetOffline(training_set &TrainingSet);
		void ResilientPropagation();
		void ResilientPropagationOffline();
		double RMSTraining(training_set &TrainingSet);
		void RMSPropagation();

		//Calc
		double CalculateError(Problem& test, bool print = false);
		void CalcGradDelta(const double output);
		void CalcGradDelta(const std::vector<double>& output);
		void CalcGradDelta(const Matrix2d &output);

		//Gets
		void PrintFullInfo(std::ostream &debug) const;
		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);

		//Sets
		void SetInput(const Matrix2d& inputs);
	};
}
