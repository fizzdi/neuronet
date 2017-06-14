#pragma once

#include "Layer.h"
#include <vector>
#include <iostream>

namespace NeuroNet
{
	struct Problem
	{
		Matrix2d inputs, outputs;
		Problem(std::vector<double> input, std::vector<double> output)
		{
			inputs = input;
			outputs = output;
		}
	};

	enum LearningType { RPROP, PROP };


	class NeuralNetwork
	{
	protected:
		std::vector<Layer> _layers;
		LearningType _learn;
	public:
		void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn);
		double RunTrainingSet(bool print = false);
		std::vector<Problem> TrainingSet;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);
		void PrintProblemResult(Problem& test);
		double CalculateError(Problem& test, bool print = false);
		void CorrectWeights();
		void CalcCorrectWeights(Problem& test);
		void ResilientPropagation(Problem& test);
		void BackPropagation(Problem& test);

		Matrix2d Run(Matrix2d &inputs);
	};
}