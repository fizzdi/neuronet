#pragma once

#include "Layer.h"
#include <vector>
#include <iostream>
const double EducationalSpeed = 0.0001;

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

	class NeuralNetwork
	{
	protected:
		std::vector<Layer> _layers;
	public:
		void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		double RunTrainingSet(bool print = false);
		std::vector<Problem> TrainingSet;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);
		void PrintProblemResult(Problem& test);
		double CalculateError(Problem& test, bool print = false);
		void CorrectWeights();
		void CalcCorrectWeights(Problem& test);

		Matrix2d Run(Matrix2d &inputs);
	};
}