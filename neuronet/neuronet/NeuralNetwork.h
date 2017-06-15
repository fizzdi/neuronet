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
		std::vector<Matrix2d> LastDeltaSum;
		std::vector<Matrix2d> LastGradSum;
	public:
		std::vector<Problem> TrainingSet;
		double Xmin, Xmax, Ymin, Ymax;

		void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction, LearningType Learn);
		void Run(Problem test);
		double RunTrainingSet(bool print = false);
		double RunTrainingSetOffline(bool print = false);
		double CalculateError(Problem& test, bool print = false);
		void PrintProblemResult(Problem& test);
		void CalcCorrectWeights(Problem& test);
		void ResilientPropagation(Problem& test);
		void ResilientPropagation(std::vector<Matrix2d> PrevSumGrad, std::vector<Matrix2d> SumGrad, std::vector<Matrix2d> PrevSumDelta, std::vector<Matrix2d> SumDelta);
		void BackPropagation(Problem& test);


		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, NeuralNetwork &net);
	};
}