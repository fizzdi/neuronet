#pragma once

#include "Layer.h"
#include <vector>
#include <iostream>
const double EducationalSpeed = 0.3;
const double Alpha = 0.01;

namespace NeuroNet
{
	struct Problem
	{
		std::vector<double> inputs, outputs;
		Problem(std::vector<double> input, std::vector<double> output)
		{
			inputs = input;
			outputs = output;
		}
	};

	class NeuralNetwork
	{
	private:
		std::vector<Layer> _layers;
	public:
		void Init(int InputCount, int OutputCount, int NeuronCount);
		double RunProblemSet();
		std::vector<Problem> ProblemSet;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);
		void PrintProblemResult(const Problem& test);
		double CalculateError(const Problem& test);
		void CorrectWeights();
		void CalcCorrectWeights(const Problem& test);

		std::vector<double> Run(std::vector<double> &inputs);
	};
}