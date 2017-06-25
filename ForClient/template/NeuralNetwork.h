#pragma once
#include "Layer.h"

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
		virtual void Run() = 0;
	public:
		//Methods
		virtual void Run(const std::vector<double>& input) = 0;
		virtual void Run(const Matrix2d &input) = 0;

		//Calc
		virtual double CalculateError(Problem& test, bool print = false) = 0;
		virtual void CalcGradDelta(const double output) = 0;
		virtual void CalcGradDelta(const std::vector<double>& output) = 0;
		virtual void CalcGradDelta(const Matrix2d &output) = 0;

		//Gets
		virtual	void PrintFullInfo(std::ostream &debug) const = 0;
		virtual Matrix2d GetOut() const = 0;

		//Sets
		virtual void SetInput(const Matrix2d& inputs) = 0;
	};
}
