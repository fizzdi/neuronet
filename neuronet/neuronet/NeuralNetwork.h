#pragma once

#include "Layer.h"
#include <vector>
#include <iostream>
const double EducationalSpeed = 0.001;

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

	enum LearningType { RPROP, BACKPROP };

	class NeuralNetwork
	{
	protected:
		std::vector<Layer> _layers;
		LearningType _learn;
	public:
		std::vector<Problem> TrainingSet;

		void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn);
		double RunTrainingSet(bool print = false);
		double CalculateError(Problem& test, bool print = false);
		void Run(const Matrix2d &inputs);
		void PrintProblemResult(Problem& test);
		void CalcCorrectWeights(Problem& test);
		void ResilientPropagation(Problem& test);
		void BackPropagation(Problem& test);


		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);
	};
}